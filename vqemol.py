import sys
#print(sys.path)
from openfermion.chem import MolecularData
from openfermionpyscf import prepare_pyscf_molecule
from openfermionpyscf import run_pyscf
from qiskit import Aer,aqua,execute,QuantumCircuit
from qiskit.visualization.utils import _validate_input_state
from math import pi,degrees,acos,sqrt,cos
from joblib import Parallel,delayed
import scipy.optimize
import numpy
import cmath
from ANSATZ import (uccsd_circuit,he_circuit,combination,
                    mp2_uccsd_init,zero_uccsd_init, random_uccsd_init,
                    zero_he_init, random_he_init, hf_state,
                    zero_he_type2_init, random_he_type2_init,he_circuit_type2)
from MEASUREMENT import statevector_expect
from HAMILTONIAN import compute_spin_integral, compute_scf, ph_hamiltonian, hamiltonian_grouping 
import copy 

class VQEMOL:

    def __init__(self,molecule):

        self.opf_mol = molecule
        self.scf_mol = prepare_pyscf_molecule(molecule)
        self.ele_num = self.opf_mol.get_n_alpha_electrons()+self.opf_mol.get_n_beta_electrons()
        self.aodip = self.scf_mol.intor_symmetric('int1e_r')
        self.aoint2 = self.scf_mol.intor('int2e')
        self.hf_mol = compute_scf(self.scf_mol)
        self.hcore = self.hf_mol.get_hcore()
        pyscf_mol = run_pyscf(molecule)
        self.mo_a = pyscf_mol.canonical_orbitals
        self.one_integral,self.two_integral = pyscf_mol.get_integrals()
        self.qubit_num = 2*self.one_integral.shape[0]
        self.g_jw_hamil = compute_spin_integral(self.mo_a,self.mo_a,self.hcore,self.aoint2)
        self.g_jw_opr, self.g_jw_num = hamiltonian_grouping(list(self.g_jw_hamil.terms.keys()),
                                                            list(self.g_jw_hamil.terms.values()))
        self.s_comb,self.d_comb = combination(self.ele_num,self.qubit_num)

        self.machine = Aer.get_backend('statevector_simulator')


    def get_HF_result(self):
        hf_keys =['replusion','hf_energy','orbital energy']
        hf_values = [self.scf_mol.energy_nuc(),self.scf_mol.e_tot,self.scf_mol.mo_energy]
        return  dict(zip(hf_keys, hf_values))

     
            
    def runVQE(self,ansatz,init_guess='ZERO',use_core=None,vqe_param=None,t_num=1,d_num=1,
               opt_detail={'Method':'BFGS','gtol':1e-4,'eps':1e-5,'maxiter':None},
               EPR=False,opt_log=False):


        if use_core is None:
            ncpu = 1
        else:
            ncpu = use_core
        
        if opt_log:
            energy_log = []
            ele_num_log = []
            param_log = []

        if EPR:
            opt_detail['Method']='trust-constr'
        else:
            opt_detail['Method']='BFGS'
        self.t_num = t_num
        self.d_num = d_num
        self.ansatz=ansatz
        def cons_elenum(phi):
            if self.ansatz == 'UCCSD':
                vqe_state=0
                vector = 0
                vqe_state = hf_state(self.ele_num,self.qubit_num)
                vector = uccsd_circuit(self.ele_num,self.qubit_num,vqe_state,phi,trotter=self.t_num)
            elif self.ansatz == 'HE_TYPE1':
                vector = 0
                vector = he_circuit(self.qubit_num,phi,depth=self.d_num)
            elif self.ansatz == 'HE_TYPE2':
                vector = 0
                vector = he_circuit_type2(self.qubit_num,phi,depth=self.d_num)
            
            alpha_num,beta_num = self.count_elenum(vector=vector)
            return self.ele_num -(sum(alpha_num)-sum(beta_num))

        if ansatz=='UCCSD':
            if vqe_param is None:
                if init_guess=='ZERO':
                    self.init_param = zero_uccsd_init(self.ele_num,self.qubit_num)
                elif init_guess=='RANDOM':
                    self.init_param = random_uccsd_init(self.ele_num,self.qubit_num)
                elif init_guess=='MP2':
                    self.init_param = mp2_uccsd_init(self.ele_num,self.qubit_num,self.mo_a,
                                                     self.mo_a,self.hcore,self.aoint2)
            else:
                self.init_param = vqe_param

            def uccsd_cost(phi):
                vqe_state = 0
                vector = 0
                vqe_state = hf_state(self.ele_num,self.qubit_num)
                vector = uccsd_circuit(self.ele_num,self.qubit_num,vqe_state,phi,trotter=t_num)
                def uccsd_hamil_cost(hamil_num,ucc_wave):
                    return statevector_expect(ucc_wave,self.g_jw_opr[hamil_num],self.g_jw_num[hamil_num],
                                              self.qubit_num,hamil_num,len(self.g_jw_opr))
                hamil_processed = Parallel(n_jobs=ncpu,backend="threading",verbose=0)([delayed(uccsd_hamil_cost)(i,vector) 
                                                                                       for i in range(len(self.g_jw_opr))])
                if opt_log:
                    alpha,beta = self.count_elenum(vector=vector)
                    energy_log.append(sum(hamil_processed))
                    ele_num_log.append(dict(zip(['alpha','beta','total'],[alpha,beta,sum(alpha)+sum(beta)])))
                    param_log.append(phi)
                #print(sum(hamil_processed))
                return sum(hamil_processed)
                        
            if EPR:                
                cons_ele = ({'type':'eq','fun':cons_elenum})
                res = scipy.optimize.minimize(uccsd_cost,self.init_param,constraints=cons_ele,
                                              method = opt_detail['Method'],
                                              options = {'xtol': 1e-04, 'gtol': 1e-04, 'barrier_tol': 1e-04})
           
            else:
                res = scipy.optimize.minimize(uccsd_cost,self.init_param,method = opt_detail['Method'],
                                              options = {"gtol": opt_detail['gtol'],
                                                         'maxiter': opt_detail['maxiter'],
                                                         'eps': opt_detail['eps']})
            
        elif ansatz=='HE_TYPE1':
            if vqe_param is None:
                if init_guess=='ZERO':
                    self.init_param = zero_he_init(self.qubit_num,d_num)
                elif init_guess=='RANDOM':
                    self.init_param = random_he_init(self.qubit_num,d_num)
            else:
                self.init_param = vqe_param

            def he_cost(phi):
                vqe_state = 0
                vqe_state = hf_state(self.ele_num,self.qubit_num)
                vector = 0
                vector = he_circuit(self.qubit_num,phi,depth=d_num)
                def he_hamil_cost(hamil_num,he_wave):
                    return statevector_expect(he_wave,self.g_jw_opr[hamil_num],self.g_jw_num[hamil_num],
                                              self.qubit_num,hamil_num,len(self.g_jw_opr))
                hamil_processed = Parallel(n_jobs=ncpu,backend="threading",verbose=0)([delayed(he_hamil_cost)(i,vector)
                                                                                       for i in range(len(self.g_jw_opr))])
                #print(sum(hamil_processed))
                
                if opt_log:
                    alpha,beta = self.count_elenum(vector=vector)
                    energy_log.append(sum(hamil_processed))
                    ele_num_log.append(dict(zip(['alpha','beta','total'],[alpha,beta,sum(alpha)+sum(beta)])))
                    param_log.append(phi)
                return sum(hamil_processed)
            
            if EPR:                
                cons_ele = ({'type':'eq','fun':cons_elenum})
                res = scipy.optimize.minimize(he_cost,self.init_param,constraints=cons_ele,
                                              method = opt_detail['Method'],
                                              options = {'xtol': 1e-04, 'gtol': 1e-04, 'barrier_tol': 1e-04})
                
            else:
                res = scipy.optimize.minimize(he_cost,self.init_param,method = opt_detail['Method'],
                                              options = {"gtol": opt_detail['gtol'],
                                                         'maxiter': opt_detail['maxiter'],
                                                         'eps': opt_detail['eps']})
            
        elif ansatz=='HE_TYPE2':
            if vqe_param is None:
                if init_guess=='ZERO':
                    self.init_param = zero_he_type2_init(self.qubit_num,d_num)
                elif init_guess=='RANDOM':
                    self.init_param = random_he_type2_init(self.qubit_num,d_num)
            else:
                self.init_param = vqe_param

            def he_cost_2(phi):
                vqe_state = 0
                vqe_state = hf_state(self.ele_num,self.qubit_num)
                vector = 0
                vector = he_circuit_type2(self.qubit_num,phi,depth=d_num)
                def he_hamil_cost_2(hamil_num,he_wave):
                    return statevector_expect(he_wave,self.g_jw_opr[hamil_num],self.g_jw_num[hamil_num],
                                              self.qubit_num,hamil_num,len(self.g_jw_opr))
                hamil_processed = Parallel(n_jobs=ncpu,backend="threading",verbose=0)([delayed(he_hamil_cost_2)(i,vector)
                                                                                       for i in range(len(self.g_jw_opr))])
                if opt_log:
                    alpha,beta = self.count_elenum(vector=vector)
                    energy_log.append(sum(hamil_processed))
                    ele_num_log.append(dict(zip(['alpha','beta','total'],[alpha,beta,sum(alpha)+sum(beta)])))
                    param_log.append(phi)
                return sum(hamil_processed)

                        
            if EPR:                
                cons_ele = ({'type':'eq','fun':cons_elenum})
                res = scipy.optimize.minimize(he_cost_2,self.init_param,constraints=cons_ele,
                                              method = opt_detail['Method'],
                                              options = {'xtol': 1e-04, 'gtol': 1e-04, 'barrier_tol': 1e-04})
                
            else:
                res = scipy.optimize.minimize(he_cost_2,self.init_param,method = opt_detail['Method'],
                                              options = {"gtol": opt_detail['gtol'],
                                                         'maxiter': opt_detail['maxiter'],
                                                         'eps': opt_detail['eps']})
        self.opt_param = res.x
        vqe_keys = ['success','opt_param','energy','nfev','nit']
        vqe_values = [res.success,res.x,res.fun+self.scf_mol.energy_nuc(),res.nfev,res.nit]
        if opt_log:
            log_keys =['energy_log','param_log','elenum_log']
            log_values = [energy_log,param_log,ele_num_log]
            return dict(zip(vqe_keys,vqe_values)), dict(zip(log_keys,log_values))
        else:
            return dict(zip(vqe_keys,vqe_values))
            
    def cal_NaturalOrb(self):
        if self.ansatz == 'UCCSD':
            vqe_state=0
            vector = 0
            vqe_state = hf_state(self.ele_num,self.qubit_num)
            vector = uccsd_circuit(self.ele_num,self.qubit_num,vqe_state,self.opt_param,trotter=self.t_num)
            self.opt_vector = vector
        elif self.ansatz == 'HE_TYPE1':
            vector = 0
            vector = he_circuit(self.qubit_num,self.opt_param,depth=self.d_num)
            self.opt_vector = vector
        elif self.ansatz == 'HE_TYPE2':
            vector = 0
            vector = he_circuit_type2(self.qubit_num,self.opt_param,depth=self.d_num)
            self.opt_vector = vector

        rho = _validate_input_state(self.opt_vector)
        num_rho = int(numpy.log2(len(rho)))
        state = [bin(i)[2:].zfill(num_rho) for i in range(2**num_rho)]
        state_keys = []
        for ns in state:
            state_keys.append(str(ns))
            
        return dict(zip(state_keys,self.opt_vector))

    def count_elenum(self,vector=None):
        if vector is None:
            count = copy.copy(self.opt_vector)
        else:
            count = copy.copy(vector)

        for r in range(len(count)):
            count[r] = (abs(count[r]) ** 2).real

        rho = _validate_input_state(count)
        num_rho = int(numpy.log2(len(rho)))
        state = [bin(i)[2:].zfill(num_rho) for i in range(2**num_rho)]
        alpha_elec = numpy.zeros(self.qubit_num//2)
        beta_elec = numpy.zeros(self.qubit_num//2)

        for r in range(len(state)):
            temp_alpha= []
            temp_beta = []
            temp_list = list(state[r])
            for j in range(len(temp_list)):
                num = int(temp_list[j])
                if j % 4 == 0 or j % 4 == 2:
                    temp_alpha.append((num * count[r]).real)
                    #alpha_elec = alpha_elec + temp_alpha
                else:
                    temp_beta.append((num * count[r]).real)
                    #beta_elec = beta_elec + temp_beta
            if len(temp_alpha) != 0:
                alpha_elec = alpha_elec + temp_alpha
            if len(temp_beta) != 0:
                beta_elec = beta_elec + temp_beta

        #print(alpha_elec)
        r_alpha_elec = numpy.flip(alpha_elec)
        r_beta_elec = numpy.flip(beta_elec)

        alpha = sum(r_alpha_elec)
        beta = sum(r_beta_elec)

        num = alpha + beta
        return r_alpha_elec,r_beta_elec

"""
basis = "sto-3g" #basis set
multiplicity = 1 #spin multiplicity
charge = 0   #total charge for the molecule
geometry = [("H",(0.37,0,0)),("H",(-0.37,0,0))]
molecule = MolecularData(geometry, basis, multiplicity, charge)
mol = VQEMOL(molecule)

ans = mol.runVQE('UCCSD',init_guess='MP2'))
opt_param = ans['opt_param']
"""

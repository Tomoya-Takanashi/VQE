import numpy
from openfermionpyscf import prepare_pyscf_molecule
from openfermion.transforms import jordan_wigner,get_fermion_operator
from openfermion.ops import InteractionOperator
from openfermion.ops import QubitOperator
from pyscf import scf

def compute_scf(pyscf_molecule,UHF=False):
    """
    Perform a Hartree-Fock calculation.
    Args:
    pyscf_molecule: A pyscf molecule instance.
    Returns:
    pyscf_scf: A PySCF "SCF" calculation object.
    """
    if UHF:
        pyscf_scf = scf.UHF(pyscf_molecule)
    else:
        if pyscf_molecule.spin:
            pyscf_scf = scf.ROHF(pyscf_molecule)
        else:
            pyscf_scf = scf.RHF(pyscf_molecule)

    return pyscf_scf

def compute_spin_integral(mo_a,mo_b,hcore,ao2int,central=None,EDA=False,mp2=False):
    
    repulsion = 0
    ao1int = hcore
    aolen = len(ao2int)
    n_spin = aolen * 2    
    moint = 0
    moint1 = numpy.zeros((n_spin,n_spin))
    moint2 = numpy.zeros((n_spin,n_spin,n_spin,n_spin))
    moint_fock = numpy.zeros((n_spin,n_spin))
    st = 2

    if central is None:
        pass 
    else:
        central_s = central[0]
        central_e = central[len(central)-1]+1

    for i in range(n_spin):
        for a in range(n_spin):
            for a1 in range(aolen):
                if EDA:
                    for a2 in range(central_s,central_e,1):
                        p = i // 2
                        q = a // 2
                        
                        if i % 2 == 0 and a % 2 == 0 :
                            moint = moint + mo_a[a1][p]*mo_a[a2][q]*ao1int[a1][a2]
                        elif i % 2 == 0 and a % 2 != 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 == 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 != 0:
                            moint = moint + mo_b[a1][p]*mo_b[a2][q]*ao1int[a1][a2]

                else:
                    for a2 in range(aolen):
                        p = i // 2
                        q = a // 2
                        if i % 2 == 0 and a % 2 == 0 :
                            moint = moint + mo_a[a1][p]*mo_a[a2][q]*ao1int[a1][a2]
                        elif i % 2 == 0 and a % 2 != 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 == 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 != 0:
                            moint = moint + mo_b[a1][p]*mo_b[a2][q]*ao1int[a1][a2]
            
            moint1[i][a] = moint
            moint = 0

    moint = 0
    #print(moint1)
    for p in range(n_spin):
        for q in range(n_spin):
            for r in range(n_spin):
                for s in range(n_spin):
                    for a1 in range(aolen):
                        for a2 in range(aolen):
                            for a3 in range(aolen):
                                if EDA:
                                    for a4 in range(central_s,central_e,1):
                                        if p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 ==0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 == 0 and q % 2 == 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 == 0 and s % 2 == 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        else:
                                            moint = 0

                                else:
                                    for a4 in range(aolen):
                                        if p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 ==0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 == 0 and q % 2 == 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 == 0 and s % 2 == 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        else:
                                            moint = 0
                    moint2[p][q][r][s] = moint/2
                    moint = 0

    EQ_TOLERANCE = 0.000000001
    moint1[numpy.absolute(moint1) < EQ_TOLERANCE] = 0.
    moint2[numpy.absolute(moint2) < EQ_TOLERANCE] = 0.
    if mp2:
        moint1_mp2 = moint1
        moint2_mp2 = moint2
    moint2 =numpy.asarray(moint2.transpose(0, 2, 3, 1), order='C')
    
    molecular_hamiltonian = InteractionOperator(repulsion, moint1, moint2)
    #print(N_hamiltonian(repulsion,moint1,moint2,2))

#    print(molecular_hamiltonian)
    #fermion_opr=get_fermion_operator(molecular_hamiltonian)
    #print('test')
    #print(jordan_wigner_fermion_operator(fermion_opr))
    jw_hamiltonian = jordan_wigner(get_fermion_operator(molecular_hamiltonian))
    if mp2:
        return moint1_mp2, moint2_mp2
    else:
        return jw_hamiltonian

def ph_hamiltonian(mo_a,mo_b,hcore,ao2int,hf_ene,fock_a,fock_b,central=None,sub=None,EDA=False):
    aolen = len(hcore)
    n_spin = aolen * 2    
    moint = 0
    moint_hcore = 0 
    moint1_hcore = numpy.zeros((n_spin,n_spin))
    moint1 = numpy.zeros((n_spin,n_spin))
    moint2 = numpy.zeros((n_spin,n_spin,n_spin,n_spin))
    moint1_ver = numpy.zeros((n_spin,n_spin))
    moint2_ver = numpy.zeros((n_spin,n_spin,n_spin,n_spin))
    moint2_ver_2 = numpy.zeros((n_spin,n_spin,n_spin,n_spin))
    st = 2
    if central is not None:
        central_s = min(central) - min(sub)
        central_e = max(central) - min(sub)
        
    for i in range(n_spin):
        for a in range(n_spin):
            for a1 in range(aolen):
                if EDA:
                    for a2 in range(central_s,central_e,1):
                        p = i // 2
                        q = a // 2
                        if i % 2 == 0 and a % 2 == 0 :
                            moint = moint + mo_a[a1][p]*mo_a[a2][q]*fock_a[a1][a2]
                            moint_hcore = moint_hcore + mo_a[a1][p]*mo_a[a2][q]*hcore[a1][a2]
                        elif i % 2 == 0 and a % 2 != 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 == 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 != 0:
                            moint = moint + mo_b[a1][p]*mo_b[a2][q]*fock_b[a1][a2]
                            moint_hcore = moint_hcore + mo_b[a1][p]*mo_b[a2][q]*hcore[a1][a2]
                else:
                    for a2 in range(aolen):
                        p = i // 2
                        q = a // 2
                        if i % 2 == 0 and a % 2 == 0 :
                            moint = moint + mo_a[a1][p]*mo_a[a2][q]*fock_a[a1][a2]
                            moint_hcore = moint_hcore + mo_a[a1][p]*mo_a[a2][q]*hcore[a1][a2]
                        elif i % 2 == 0 and a % 2 != 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 == 0:
                            moint = 0
                        elif i % 2 != 0 and a % 2 != 0:
                            moint = moint + mo_b[a1][p]*mo_b[a2][q]*fock_b[a1][a2]
                            moint_hcore = moint_hcore + mo_b[a1][p]*mo_b[a2][q]*hcore[a1][a2]

            moint1[i][a]=moint
            moint1_hcore[i][a] = moint_hcore
            moint = 0
            moint_hcore = 0

    moint = 0
                    
    for p in range(n_spin):
        for q in range(n_spin):
            for r in range(n_spin):
                for s in range(n_spin):
                    for a1 in range(aolen):
                        for a2 in range(aolen):
                            for a3 in range(aolen):
                                if EDA:
                                    for a4 in range(central_s,central_e,1):
                                        if p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 ==0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 == 0 and q % 2 == 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 == 0 and s % 2 == 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        else:
                                            moint = 0
                                else:
                                    for a4 in range(aolen):
                                        if p % 2 == 0 and q % 2 == 0 and r % 2 == 0 and s % 2 ==0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 == 0 and q % 2 == 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_a[a1][p//2]*mo_a[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 == 0 and s % 2 == 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_a[a3][r//2]*mo_a[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        elif p % 2 != 0 and q % 2 != 0 and r % 2 != 0 and s % 2 != 0:
                                            moint = moint + mo_b[a1][p//2]*mo_b[a2][q//2]*mo_b[a3][r//2]*mo_b[a4][s//2]*ao2int[a1][a2][a3][a4]
                                        else:
                                            moint = 0
                    moint2[p][q][r][s] = moint
                    moint = 0
    moint = 0
    temp = 0
    moint2_ver = moint2/2
    
    one_const = 0
    two_const = 0
    temp = 0
    total_const = 0
    
    for i in range(aolen):
        one_const = one_const + (moint1[i][i] -moint1_hcore[i][i])

    for i in range(aolen):
        for j in range(aolen):
            two_const = two_const + 1*moint2[i][i][j][j] + (-1)*moint2[i][j][i][j]

    total_const = one_const/2 - two_const/2
    #print(total_const)
    EQ_TOLERANCE = 0.000000001
    moint1[numpy.absolute(moint1) < EQ_TOLERANCE] = 0.
    moint2[numpy.absolute(moint2) < EQ_TOLERANCE] = 0.
    moint1_hcore[numpy.absolute(moint1_hcore) < EQ_TOLERANCE] = 0.
    moint1_ver[numpy.absolute(moint1_ver) < EQ_TOLERANCE] = 0.
    moint2_ver[numpy.absolute(moint2_ver) < EQ_TOLERANCE] = 0.
    moint2_ver =numpy.asarray(moint2_ver.transpose(0, 2, 3, 1), order='C')
    #    moint2 =numpy.asarray(moint2.transpose(0, 2, 3, 1), order='C')
    molecular_hamiltonian = InteractionOperator(total_const, moint1_hcore, moint2_ver)
    jw_hamiltonian = jordan_wigner(get_fermion_operator(molecular_hamiltonian))
    return jw_hamiltonian

def hamiltonian_grouping(hamiltonianListOpr,hamiltonianListNum):

    hamiltonian_z_Opr = []
    hamiltonian_z_Num = []
    
    for i in range(len(hamiltonianListOpr)):
        exist_Z = []
        for j in range(len(hamiltonianListOpr[i])):
            exist_Z.append('Z' in hamiltonianListOpr[i][j])

        if(False in exist_Z):
            continue
        else:
            hamiltonian_z_Opr.append(hamiltonianListOpr[i])
            hamiltonian_z_Num.append(hamiltonianListNum[i])
            hamiltonianListOpr[i] = 'dummy'
            hamiltonianListNum[i] = 'dummy'

    hamiltonianOpr = []
    hamiltonianNum = []
    for i in range(len(hamiltonianListOpr)):
        if(hamiltonianListOpr[i] != 'dummy'):
            hamiltonianOpr.append(hamiltonianListOpr[i])
            hamiltonianNum.append(hamiltonianListNum[i])

    for i in range(len(hamiltonianOpr)):
        for j in range(i+1):
            if j <= len(hamiltonianOpr):
                if len(hamiltonianOpr[i]) > len(hamiltonianOpr[j]):
                    oprtmp = hamiltonianOpr[i]
                    numtmp = hamiltonianNum[i]
                    hamiltonianOpr[i] = hamiltonianOpr[j]
                    hamiltonianNum[i] = hamiltonianNum[j]
                    hamiltonianOpr[j] = oprtmp
                    hamiltonianNum[j] = numtmp
    hamiltonian_group_Opr = []
    hamiltonian_group_Num = []
    Opr_hamiltonian = []
    Num_hamiltonian = []

    for i in range(len(hamiltonianOpr)):

        ref = hamiltonianOpr[i]
        ref_Num = hamiltonianNum[i]
            
        if ref == 'dummy':
            continue
        else:
            hamiltonian_group_Opr = []
            hamiltonian_group_Num = []
            hamiltonian_group_Opr.append(ref)
            hamiltonian_group_Num.append(ref_Num)
            for j in range(i+1,len(hamiltonianOpr),1):
                if hamiltonianOpr[j] == 'dummy':
                    continue
                else:
                    exist = []
                    if j <= len(hamiltonianOpr):
                        for k in range(len(hamiltonianOpr[j])):
                            exist.append(hamiltonianOpr[j][k] in ref)
                        if False in exist:
                            continue
                        else:
                            hamiltonian_group_Opr.append(hamiltonianOpr[j])
                            hamiltonian_group_Num.append(hamiltonianNum[j])
                            hamiltonianOpr[j] = 'dummy'
                            hamiltonianNum[j] = 'dummy'
            Opr_hamiltonian.append(hamiltonian_group_Opr)
            Num_hamiltonian.append(hamiltonian_group_Num)

    Opr_hamiltonian.append(hamiltonian_z_Opr)
    Num_hamiltonian.append(hamiltonian_z_Num)

    return Opr_hamiltonian,Num_hamiltonian

    

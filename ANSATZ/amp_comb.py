import copy

def combination(ele_num,qubit_num,sort=True):
    Comb = ele_num * (qubit_num - ele_num)
    occ = ele_num
    vir = qubit_num - occ
    occComb = []
    virComb = []
    occNum_s = []
    virNum_s = []
    single_comb = []


    for i in range(ele_num):
        occNum_s.append(i)

    for i in range(ele_num,qubit_num):
        virNum_s.append(i)
    nvir = len(virNum_s)

    occ_a=[]
    occ_b=[]
    vir_a=[]
    vir_b=[]

    for i in occNum_s:
        if i % 4 == 0 or i % 4 == 2:
            occ_a.append(i)
        else:
            occ_b.append(i)

    for i in virNum_s:
        if i % 4 == 0 or i % 4 == 2:
            vir_a.append(i)
        else:
            vir_b.append(i)
    single_comb_a = []
    single_comb_b = []
    for i in occ_a:
        for j in vir_a:
            single_comb_a.append([i,j])

    for i in occ_b:
        for j in vir_b:
            single_comb_b.append([i,j])

    single_comb = single_comb_a + single_comb_b
    occ_aa = []
    occ_ab = []
    occ_bb = []

    vir_aa = []
    vir_ab = []
    vir_bb = []

    for i in occ_a:
        for j in occ_a:
            if(i < j):
                occ_aa.append([i,j])

    for i in occ_b:
        for j in occ_b:
            if(i < j):
                occ_bb.append([i,j])

    for i in occ_a:
        for j in occ_b:
            occ_ab.append([i,j])

    for i in vir_a:
        for j in vir_a :
            if(i < j):
                vir_aa.append([i,j])

    for i in vir_b:
        for j in vir_b:
            if(i < j):
                vir_bb.append([i,j])

    for i in vir_a:
        for j in vir_b:
            vir_ab.append([i,j])


    Combnum = len(occComb) * len(virComb)
    Combnum_s = len(occNum_s) * len(virNum_s)

    double_comb_aa = []
    double_comb_ab = []
    double_comb_bb = []

    for i in range(len(occ_aa)):
        for j in range(len(vir_aa)):
            double_comb_aa.append([occ_aa[i],vir_aa[j]])

    for i in range(len(occ_bb)):
        for j in range(len(vir_bb)):
            double_comb_bb.append([occ_bb[i],vir_bb[j]])

    for i in range(len(occ_ab)):
        for j in range(len(vir_ab)):
            double_comb_ab.append([occ_ab[i],vir_ab[j]])

    double_comb = double_comb_aa + double_comb_bb + double_comb_ab
    
    if sort:
        for i in range(len(double_comb)):
            if double_comb[i][0][0] > double_comb[i][0][1]:
                temp = double_comb[i][0][0]
                double_comb[i][0][0] = double_comb[i][0][1]
                double_comb[i][0][1] = temp
            if double_comb[i][1][0] > double_comb[i][1][1]:
                temp = double_comb[i][1][0]
                double_comb[i][1][0] = double_comb[i][1][1]
                double_comb[i][1][1] = temp
        return single_comb,double_comb

    else:
        return single_comb,double_comb

import glob
import copy
import os
import shutil 
#from Bio.PDB import *
import numpy as np
from itertools import combinations,product
def get_coordinates_pdb(filename):
    V = list()

    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                for j in range(len(line)):
                    if '.' == line[j] and '.' == line[j + 8]:
                        x_index = j
                        break
                x = line[x_index - 4:x_index + 4]
                y = line[x_index + 4:x_index + 12]
                z = line[x_index + 12:x_index + 20]
                V.append(np.asarray([x, y, z], dtype=float))

    V = np.asarray(V)

    return V

def centroid(X):
    C = X.mean(axis=0)  # 计算每一列的均值
    return C


def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def rmsd(P, Q):
    diff = P - Q
    return np.sqrt((diff * diff).sum() / P.shape[0])

def rmsdcal(m, n):
    m -= centroid(m)
    n -= centroid(n)
    U = kabsch(m, n)
    m = np.dot(m, U)
    r = rmsd(m, n)
    
    return r
def assemble_frag(frags, fragsplit, ftype, num, bpn, tn): # ['0', '5_12', '12']
    fn = len(frags)
    kb_num = 10
    j_fraglist = []
    j_lines = []
    
    for f in frags:
        with open(f) as f1:
            line = f1.readlines()
        j_lines = j_lines + line
        f1.close()
    
    for i in range(fn):
        j_num = fragsplit[i].split('_')
        j_way = len(j_num)
        j_num = [int(x) for x in j_num]
        j_list = []
        fc =get_coordinates_pdb(frags[i])
        #print(len(fc), j_num)
            
        for j in range(j_way):
            #j_list.append(fc[(12*bpn*j + sum(j_num[0:j])*3): (12*bpn*(j+1) + sum(j_num[0:(j+1)])*3)])
            j_list.append(fc[((4*bpn*j + sum(j_num[0:j]))*3): ((4*bpn*(j+1) + sum(j_num[0:(j+1)]))*3)])
            #print(((4*bpn*j + sum(j_num[0:j]))*3), ((4*bpn*(j+1) + sum(j_num[0:(j+1)]))*3))
        j_fraglist.append(j_list)
        
    #print(len(j_fraglist[1][1]))
    
    # 第一次拼装
    for i in range(fn-1):
        i_all = j_fraglist[i]
        j_all = j_fraglist[i+1]
        
        i_insert = np.array(i_all[-1][-6*bpn:])
        j_insert = np.array(j_all[0][0:6*bpn])
                
        #print(i_insert.shape,j_insert.shape)
        
        for j in range(len(j_all)):
            j_ac, r =  kb_frag(j_all[j], i_insert, j_insert)
            j_fraglist[i+1][j] = j_ac
        
    #循环拼装
    j_alllist = []
    for m in range(len(j_fraglist)):
            j_alllist = j_alllist + j_fraglist[m]
            
    #print(len(j_alllist))
    frag_num = len(j_alllist)
    
    r_all = 0
    
    for p in range(kb_num):
        #print('第%s次循环：'%(p))
        for q in range(frag_num):
            i_all = j_alllist[0:q] + j_alllist[(q+1):]
            j_all = j_alllist[q]
            if (q == 0) or (q == frag_num - 1):
                i_insert = np.append(i_all[0][0:2*bpn*3], i_all[-1][-(2*bpn*3):], axis=0)
                j_insert = np.append(j_all[-(2*bpn*3):], j_all[0:2*bpn*3], axis=0)
            else:
                i_insert = np.append(j_alllist[q+1][0:2*bpn*3], j_alllist[q-1][-(2*bpn*3):], axis=0)
                j_insert = np.append(j_all[-(2*bpn*3):], j_all[0:2*bpn*3], axis=0)
            
            j_ac, r = kb_frag(j_all, i_insert, j_insert)
            j_alllist[q] = j_ac
            
            if p == kb_num - 1:
                r_all += r
    
    #print(len(j_alllist[0]),len(j_alllist[1]), len(j_alllist[2]),len(j_lines))
    lines = change_line(j_alllist, j_lines)
    #print(j_alllist, j_lines)
    #print(lines)
    
    loop_l = []
    line_n = lines[0:(2*bpn*3)]
    
    for i in  fragsplit:
        loop_l = loop_l + [int(x) for x in i.split('_')]
        
    for i in range(len(loop_l)):
        #line_n = line_n + lines[(4*i+ 1 +sum(loop_l[0:i]))*3:(4*i+3+sum(loop_l[0:(i+1)]))*3]
        line_n = line_n + lines[((4*i+2)*bpn +sum(loop_l[0:i]))*3:((4*i+4)*bpn +sum(loop_l[0:(i+1)]))*3]
        #print(((4*i+2)*bpn +sum(loop_l[0:i]))*3,((4*i+4)*bpn +sum(loop_l[0:(i+1)]))*3)
        #print((4*i+ 1 +sum(loop_l[0:i]))*3,(4*i+3+sum(loop_l[0:(i+1)]))*3)
        
    del line_n[-(2*bpn)*3:]    
    line_n = line_n + line_n[0:(bpn*3)]
    del line_n[0:(bpn*3)]
    
    if not os.path.exists('./frag_e_%s'%(tn)):
        os.mkdir('./frag_e_%s'%(tn))
 
    with open('./frag_e_%s/%s.pdb'%(tn,num), 'w+') as f:
        for l in line_n:
            f.write(l)
    f.close()
        
    return r_all
            

def change_line(j_list, line_list, decimals = 3):
    fmt = ("{:8." + str(decimals) + "f}") * 3
    all_list = []
    i = 0
    for j in j_list:
        for c in j:
            line_list[i] = line_list[i][:30] + fmt.format(c[0], c[1], c[2]) + '  1.00\n'
            i += 1
            
    
    return line_list
    
        

def get_coordinates_pdb(filename):
    V = list()

    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                for j in range(len(line)):
                    if '.' == line[j] and '.' == line[j + 8]:
                        x_index = j
                        break
                x = line[x_index - 4:x_index + 4]
                y = line[x_index + 4:x_index + 12]
                z = line[x_index + 12:x_index + 20]
                V.append(np.asarray([x, y, z], dtype=float))

    V = np.asarray(V)

    return V


def  kb_frag(j_all, i_insert, j_insert): 
        i_size = i_insert.shape[0]  # shape(0)指的就是N值
        j_size = j_insert.shape[0]

        #print(i_all[0], j_all)
        
        if not i_size == j_size:
            raise Exception("Structures not same size")

        i_coord = copy.deepcopy(i_insert)
        j_coord = copy.deepcopy(j_insert)

        i_cent = centroid(i_coord)
        j_cent = centroid(j_coord)

        i_coord -= i_cent
        j_coord -= j_cent

        U = kabsch(j_coord, i_coord)

        j_all_changed = j_all - j_cent
        j_all_changed = np.dot(j_all_changed, U)
        j_all_changed += i_cent

        j_insertc = np.dot(j_coord, U)
        r = rmsd(i_coord + i_cent, j_insertc+i_cent)
        #print(r)
        
        return j_all_changed, r    

def junction_e(junction_p, bpn, fyp, tn):
    junction_frag = junction_p.split('_')
    junction_ln = len(junction_frag)
    f_splist_list = []
    r_dict = {}
    if junction_ln == 2:
        jf = 'two_way'
    else:
        jf = 'mutli_way'
    
    nn = max([len(x) for x in junction_frag])
    ff_list = []
    for m in junction_frag:
        ff_list.append(glob.glob('./database_frag/%s/%sbp_frag/%s/*.pdb'%(jf, bpn, m,)))
    ffe_list = product(*ff_list)
    fnum = 0
    for m in ffe_list:
        r = assemble_frag(m, junction_frag, junction_p, fnum, bpn, tn)
        r_dict[fnum] = r
        fnum+=1
    if not os.path.exists(fyp):
        os.makedirs(fyp)
    num = 1
    with open(fyp+'seqinfo.txt','w+') as f:
        si = '%s.pdb  X'%(num)
        f.write(si)
    f.close()
    shutil.copy('./frag_e_%s/%s.pdb'%(tn, num),fyp+'%s.pdb'%(num))
    #os.remove('./frag_e/')
    shutil.rmtree('./frag_e_%s/'%(tn))
    #print(r_dict)
    return r_dict[num]

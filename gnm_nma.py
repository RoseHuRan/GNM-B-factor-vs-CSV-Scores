# can run! all is ok!
# half matrix
from pymol import stored
import numpy as np
import copy

def gnm_nma(obj='1crn'):
    pymol.cmd.reinitialize()
    cmd.fetch(obj)
    cmd.remove(obj+' and not (alt ""+A)')  # keep only alternate conformation A                  
    cmd.alter(obj,'alt=""')                # reset PDB information  
    cmd.select('ca', obj+' and n. ca and pol.') 
    model = cmd.get_model('ca')
    size=len(model.atom)

    matrix=np.zeros((size,size))

    calist = []# index of CA atoms

# get index of CA atoms
    for at in model.atom:
        if at.name == 'CA':
            i = at.index
            calist.append(i)
    print('Total residues ',len(calist))

# generate the GNM matrix H
    for i in calist:
        row = calist.index(i)
        cmd.select('cur', obj+' and idx. %d'%i)
        cmd.select('range', obj+' within 12 of cur')
        cmd.select('rangeca', 'range and n. ca')
        cmd.select('nei', 'rangeca and (not cur)')
        stored.NEI=[]
        cmd.iterate('nei', 'stored.NEI.append(index)')
        matrix[row][row] = len(stored.NEI)
        for j in stored.NEI:
            matrix[row][calist.index(j)] = -1

# eigenvector
    eg = np.linalg.eigh(matrix)
    egv = eg[1]

# visualization of normal mode
    obj2=obj+'_2'
    obj3=obj+'_3'
    obj4=obj+'_4'
    cmd.create(obj2, obj)
    cmd.create(obj3, obj)
    cmd.create(obj4, obj)

    i=-1
    for at in model.atom:
        resi = at.resi
        i = i + 1
        if float(egv[i][1]) > 0:
            cmd.color('red', obj+' and resi '+resi)
        else:
            cmd.color('blue', obj+' and resi '+resi)
        
        if float(egv[i][2]) > 0:
            cmd.color('red', obj2+' and resi '+resi)
        else:
            cmd.color('blue', obj2+' and resi '+resi)
        
        if float(egv[i][3]) > 0:
            cmd.color('red', obj3+' and resi '+resi)
        else:
            cmd.color('blue', obj3+' and resi '+resi)
        if float(egv[i][4]) > 0:
            cmd.color('red', obj4+' and resi '+resi)
        else:
            cmd.color('blue', obj4+' and resi '+resi)

    cmd.enable(obj)
    cmd.enable(obj2)
    cmd.enable(obj3)
    cmd.enable(obj4)
    cmd.show_as('cartoon',obj)
    cmd.show_as('cartoon',obj2)
    cmd.show_as('cartoon',obj3)
    cmd.show_as('cartoon',obj4)
    cmd.set('grid_mode', 1)
    cmd.set('grid_slot', 1, obj)
    cmd.set('grid_slot', 2, obj2)
    cmd.set('grid_slot', 3, obj3)
    cmd.set('grid_slot', 4, obj4)
cmd.extend('gnm_nma',gnm_nma)

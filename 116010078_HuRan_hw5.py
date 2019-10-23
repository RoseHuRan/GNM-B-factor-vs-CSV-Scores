import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

def loadPDB(p_filename):
    file = open(p_filename,'r')
    lines = file.readlines()
    ca_lines = []
    for line in lines:
        if line.startswith('ATOM') and ('CA' in line):
            ca_lines.append(line) 
    return ca_lines

def getPosition(p_lines):
    CA_lines_parts = [CA_line.split() for CA_line in p_lines]
    CA_positions = []
    num = 0 
    while num < len(CA_lines_parts):
        if len(CA_lines_parts[num]) < 6:
            num += 1
            continue
        try:
            test = float(CA_lines_parts[num][10])
        except:    
            CA_lines_parts[num][10] = CA_lines_parts[num][9][4:]
        CA_position = []
        if CA_lines_parts[num][10] != '' and (CA_lines_parts[num][4] == 'A'):  # change the chain here
            if num>0 and (CA_lines_parts[num-1][5] != CA_lines_parts[num][5]):
                for index in [5,6,7,8,10]:  
                    CA_position.append(float(CA_lines_parts[num][index]))
                CA_positions.append(CA_position)
        num += 1    
    return CA_positions

def getDistances(p_CA_positions):
    CA_distances = []
    for atom_num in range(len(p_CA_positions)):
        CA_distance = []
        atom1 = p_CA_positions[atom_num]
        for index in range(len(p_CA_positions)):
            if atom_num == index:
                distance = 0
            else:
                atom2 = p_CA_positions[index]    				
                distance = math.sqrt(pow(atom1[1]-atom2[1],2)+pow(atom1[2]-atom2[2],2)+pow(atom1[3]-atom2[3],2))
            CA_distance.append(distance)
        CA_distances.append(CA_distance)
    return CA_distances

def getHmatrix(p_CA_distances, r_c):
    dist = p_CA_distances
    for i in range(len(p_CA_distances)):
        for n in range(len(p_CA_distances[i])):
            if i != n:
                if p_CA_distances[i][n] > r_c:
                    dist[i][n] = 0
                else:
                    dist[i][n] = -1   
        dist[i][i] = -sum(dist[i])
    Hmatrix = np.array(dist)
    return Hmatrix    

def getGNMBfactors(p_Hmatrix):
    GNM_b_factors = np.diag(np.linalg.pinv(p_Hmatrix))
    return GNM_b_factors

def compareBfactors(p_GNM_b_factors, p_consurf_score, p_CA_positions, p_protein, cutoff):
    p_experimental_b_factors = [parts[4] for parts in p_CA_positions]
    N = len(p_experimental_b_factors)
    exp_std = np.std(p_experimental_b_factors)
    exp_average = np.mean(p_experimental_b_factors)
    normalized_experimental_b_factors = (np.array(p_experimental_b_factors) - np.array(N*[exp_average]))/exp_std
    comp_std = np.std(p_GNM_b_factors)
    comp_average = np.mean(p_GNM_b_factors)
    normalized_GNM_b_factors = (np.array(p_GNM_b_factors) - np.array(N*[comp_average]))/comp_std
    gnm_correlation_coef = pearsonr(p_consurf_score, normalized_GNM_b_factors)[0]
    xray_correlation_coef = pearsonr(p_consurf_score, normalized_experimental_b_factors)[0]

# # draw the graph of correlation coefficient, GNM and X-ray B-factor profiles
#     residue_num = [parts[0] for parts in p_CA_positions]
#     plt.figure(figsize=(20,4))
#     plt.title('Conservation scores and GNM B-factor profiles of {} (cut-off={})'.format(p_protein, cutoff))
#     plt.plot(residue_num, normalized_GNM_b_factors, color='lightsalmon', label='GNM B-factors')
#     plt.plot(residue_num, p_consurf_score, color='lightslategray', label='Conservation Scores')
#     plt.legend()
#     plt.xlabel('Residue')
#     plt.show()

#     plt.figure(figsize=(20,4))
#     plt.title('Conservation scores and X-ray B-factor profiles of {}'.format(p_protein))
#     plt.plot(residue_num, normalized_experimental_b_factors, color='cyan', label='X-ray B-factors')
#     plt.plot(residue_num, p_consurf_score, color='lightslategray', label='Conservation Scores')
#     plt.legend()
#     plt.xlabel('Residue')
#     plt.show()
    
    return (gnm_correlation_coef, xray_correlation_coef)

def smooth(p_consurf_score):
    index = 0
    consurf_score_smooth = []
    for consurf_score in p_consurf_score:
        if index == 0:
            consurf_score_smooth.append((p_consurf_score[index] + p_consurf_score[index+1])/2)
        elif index == len(p_consurf_score)-1:
            consurf_score_smooth.append((p_consurf_score[index-1] + p_consurf_score[index])/2)
        elif index == 1 or index == len(p_consurf_score)-2:
            consurf_score_smooth.append((p_consurf_score[index-1] + p_consurf_score[index] + p_consurf_score[index+1])/3)
        else:
            consurf_score_smooth.append((p_consurf_score[index-2] + p_consurf_score[index-1] + p_consurf_score[index] \
             + p_consurf_score[index+1] + p_consurf_score[index+2])/5)
        index += 1
    return consurf_score_smooth 

def main():
    protein_list = ['1b2m', '1ge7', '1eyp', '1cd5', '2hdh', '1ebf', '1dfo', '1ysc', '1gog', '1gpa']
    protein_index = ['3', '53', '104', '179', '219', '338', '401' ,'403', '518', '542']
    head = ''.ljust(12)+'Index'.ljust(9)+'Diameter'.ljust(54)+'GNM'.ljust(51)+'X-ray   '+'GNM_max  '+'cut-off_max'+'\n'
    head += 'Cut-off'.ljust(33)
    for i in range(7, 22):
        head += str(i).ljust(6)
    print(head)  
    for protein in protein_list:
        pdb_file = loadPDB(protein+'.pdb')   # change the pdb file here
        caInfo_all = getPosition(pdb_file)
        caInfo, coef, coef_smooth = [], [], []
        consurf_index = [i[0] for i in getPosition(loadPDB(protein+'_With_Conservation_Scores.pdb'))]
        consurf_score = [i[4] for i in getPosition(loadPDB(protein+'_With_Conservation_Scores.pdb'))]
        consurf_score_smooth = smooth(consurf_score)
        for i in consurf_index:
            for info in caInfo_all:
                if info[0] == i:
                    caInfo.append(info)
                    break
        for radius in range(7,22):
            dist = getDistances(caInfo)
            diameter = max([max(i) for i in dist])
            Hmatrix = getHmatrix(dist, radius)  
            gnm_b_factors = getGNMBfactors(Hmatrix)
            coef.append(compareBfactors(gnm_b_factors, consurf_score, caInfo, protein.upper()+'A', radius))
            coef_smooth.append(compareBfactors(gnm_b_factors, consurf_score_smooth, caInfo, protein.upper()+'A (smoothed)', radius))
        gnm_max = max([i[0] for i in coef])
        gnm_max_smooth = max([i[0] for i in coef_smooth])
        radius_max = [i[0] for i in coef].index(gnm_max)+7
        radius_max_smooth = [i[0] for i in coef_smooth].index(gnm_max_smooth)+7
        gnm_coef = [i[0] for i in coef]
        gnm_coef_smooth = [i[0] for i in coef_smooth]

        # # draw the graph of the correlation coefficient under different cut_off radii
        # plt.title('Correlation Coefficient Under Different Cut_off Radii of {}'.format(protein.upper()+'A'))
        # plt.plot(range(7,22), gnm_coef, color='tomato', label='GNM cor-coef')
        # plt.plot(range(7,22), gnm_coef_smooth, color='mediumseagreen', label='GNM cor-coef (smoothed)')
        # plt.legend()
        # plt.xlabel('Cut-off')
        # plt.ylabel('Correlation Coefficient')
        # plt.show()

        results = '\n'+(protein.upper()+'A').ljust(12)+protein_index[protein_list.index(protein)].ljust(9)+('%.2f' % diameter).ljust(12)
        for i in range(0,14):
            results += ('%.2f' % coef[i][0]).ljust(6)
        results += ('%.2f' % coef[14][0]).ljust(9)+('%.2f' % coef[0][1]).ljust(8)+('%.2f' % gnm_max).ljust(9)+str(radius_max)+'\n'
        results += "(Smooth)".ljust(33)
        for i in range(0,14):
            results += ('%.2f' % coef_smooth[i][0]).ljust(6)
        results += ('%.2f' % coef_smooth[14][0]).ljust(9)+('%.2f' % coef_smooth[0][1]).ljust(8)+('%.2f' % gnm_max_smooth).ljust(9)+str(radius_max_smooth)
        print(results)
        
main()          
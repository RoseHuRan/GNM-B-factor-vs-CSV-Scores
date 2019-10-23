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
        if line.startswith('ATOM') and ('CA' == line.split()[2]) and ('A' == line.split()[4]):
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
        if CA_lines_parts[num][10] != '':  # change the chain here
            if num>=0 and (CA_lines_parts[num-1][5] != CA_lines_parts[num][5]):
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

def getGNMmatrix(p_Hmatrix):
    GNM_matrix = np.linalg.pinv(p_Hmatrix)
    return GNM_matrix

def getGNMcorrelation(GNM_matrix):
    GNM_correlation = 'i    j   score'
    correlation = []
    for i in range(len(GNM_matrix)):
        for j in range(i+1, len(GNM_matrix)):
            GNM_correlation += '\n' + str(i) + '    ' + str(j) + '    ' + str(GNM_matrix[i][j])
            correlation.append(GNM_matrix[i][j])             
    return GNM_correlation, correlation

def loadMSA(p_filename):
    msa_output = open(p_filename+'_msa_output.txt', 'r')
    lines = msa_output.readlines()
    msa_sequence = open(p_filename+'_msa.txt', 'r')
    sequence = msa_sequence.readlines()
    target, msa_output = [], []
    flag = 0
    for i in sequence:
        if p_filename.upper() in i:
            for residue in i.split()[1]:
                if residue != '-':
                    target.append(flag)
                flag += 1            
    for line in lines[1:]:
        i = line.split()[0]
        j = line.split()[1]
        if int(i) in target and int(j) in target:
            msa_output.append(line)
    msa_sequence_output = ''
    for i in msa_output:
        msa_sequence_output += i
    f = open(p_filename + '_sequence_msa.txt', "w")
    f.write(msa_sequence_output)

    return msa_output        

def comparison(p_correlation, p_msa):
    N = len(p_correlation)
    M = len(p_msa)
    correlation_std = np.std(np.array(p_correlation))
    correlation_mean = np.mean(np.array(p_correlation))
    msa_std = np.std(np.array(p_msa))
    msa_mean = np.mean(np.array(p_msa))
    normalized_correlation = (np.array(p_correlation) - np.array(N*[correlation_mean]))/correlation_std
    normalized_msa = (np.array(p_msa) - np.array(M*[msa_mean]))/msa_std
    msa_correlation_coef = pearsonr(normalized_msa, normalized_correlation)[0]

    # draw the graph of the correlation coefficient under different cut_off radii
    plt.title('Correlation Coefficient of relative movement and msa')
    plt.plot(range(0,N), normalized_correlation, color='tomato', label='relative movement')
    plt.plot(range(0,M), normalized_msa, color='mediumseagreen', label='msa')
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('Correlation Coefficient')
    plt.show()
    return msa_correlation_coef

def main():
    pdb_file = loadPDB('1gog.pdb')   # change the pdb file here
    msa_output = loadMSA('1goga')
    caInfo = getPosition(pdb_file)
    dist = getDistances(caInfo)
    Hmatrix = getHmatrix(dist, 22)
    GNM_matrix = getGNMmatrix(Hmatrix)
    GNM_correlation = getGNMcorrelation(GNM_matrix)[0]
    fo = open("1gog_correlation.txt", "w")
    fo.write(GNM_correlation)
    correlation = getGNMcorrelation(GNM_matrix)[1]
    msa = [float(i.split()[2]) for i in msa_output]
    print('The correlation is: ' + str(comparison(correlation, msa)))

main()

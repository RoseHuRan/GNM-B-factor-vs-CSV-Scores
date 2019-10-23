def getMSA(p_filename):
     file = open(p_filename+'_msa.aln','r')
     lines = file.readlines()
     line_list = []
     flag = 1
     for line in lines[2:]:
          if flag:
               line_list.append(line.split())
          if line == '\n':
               flag = 0
          if flag == 0 and line.split() != []:
               for sequence in line_list:
                    if line.split()[0] in sequence:
                         sequence[1] += line.split()[1]
                         break        

     output = ''
     for line in line_list:
          if line != []:
               output += line[0] + ' ' + line[1] + '\n'
     fo = open(p_filename+"_msa.txt", "w")
     fo.write(output)

getMSA('1goga')


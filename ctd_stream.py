#OM Sai Ram
#Sai
#Sai
#Feature calculation app
import streamlit as lit
from peptides import Peptide as pept
import Bio.SeqUtils.ProtParam as proparm
from propy.CTD import CalculateC as calC, CalculateT as calT, CalculateD as calD
import pandas as pd

lit.set_page_config(layout='wide')

lit.write('''
# Welcome to AMPDB Protein Feature calculation toolbox!
Tools to calculate composition features, physicochemical properties and CTD (Composition, Transition, Distribution) descriptors of protein.
''')

my_input = lit.text_input("Please enter your protein sequence/ AMPDB Acc. ID here:")
my_input = my_input.upper()
lit.markdown('<br>', unsafe_allow_html=True)
lit.write('Please select the properties you want to be calculated: ')
col1, col2 = lit.columns(2)
with col1:
    composition = lit.checkbox('Composition')
    physicochemical = lit.checkbox('Physicochemical Properties')
with col2:
    ctddescriptors = lit.checkbox('CTD Descriptors')
    others = lit.checkbox('Other Descriptors')
submit = lit.button('Submit')

if submit:
    lit.info('The results will appear below: ')
#descriptor computation,extraction and printing
if my_input and submit and not composition and not physicochemical and not ctddescriptors and not others:
    lit.error("Please select a descriptor type that you want to be calculated")
elif my_input and submit:
    if 'AMPDB_' not in my_input and my_input.isalpha() is False:
        lit.error("Some non-alphabet is present in the sequence. Please re-check!")
    elif 'AMPDB_' in my_input and my_input.replace('_', '').isalnum() is False:
        lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
    else:
        if 'AMPDB_' in my_input:
            with open('master_dataset.tsv') as file:
                l = ' '
                while(True):
                    i = file.readline()
                    if i=='':
                        lit.error('The AMPDB Acc. ID does not match with our database. Please re-check')
                        my_input = None
                        break
                    j = i.split('\t')
                    if my_input in j[1]:
                        my_input = j[3]
                        break

        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        f = open('descriptor_calculation_output.txt', 'w')
        if my_input and composition:
            lit.subheader('_Composition Features_: '); f.write("-->Composition features: \n")
            pep = pept(my_input)
            propar = proparm.ProteinAnalysis(my_input)
            count = pep.counts()
            counts = [k+':'+str(count[k]) for k in count.keys()]
            scount = [i for i in counts if int(i.split(':')[1]) != 0]
            
            frequency = pep.frequencies()
            frequencies = [k+':'+str(round(frequency[k],3)) for k in frequency.keys()]
            sfrequencie = [i for i in frequencies if float(i.split(':')[1]) != 0.0]

            f.write('Amino acid counts:\n')
            f.write('For all amino acids:\n')
            _ = [f.write(i+'\n') for i in counts]
            f.write("For the amino acids present in input sequence:\n")
            _ = [f.write(i+'\n') for i in scount]
            
            f.write('Amino Acid Frequencies: \n')
            f.write('For all amino acids:\n')
            _ = [f.write(i+'\n') for i in frequencies]
            f.write("For the amino acids present in input sequence:\n")     #"{:.2f}".format(arrs[-1])
            _ = [f.write(i+'\n') for i in sfrequencie]

            all_cnt_freq = [counts[i].split(':')+[frequencies[i].split(':')[1]]+["{:.2f} %".format(float(frequencies[i].split(':')[1])*100)]  for i in range(len(frequencies))]
            nz_cnt_freq = [scount[i].split(':')+[sfrequencie[i].split(':')[1]]+["{:.2f} %".format(float(sfrequencie[i].split(':')[1])*100)]  for i in range(len(sfrequencie))]
            
            lit.write('_Amino acid counts & frequencies_: ') 
            count_col1, count_col2 = lit.columns(2)
            with count_col1:
                lit.write("For all amino acids"); 
                lit.table(pd.DataFrame(all_cnt_freq, columns=['Amino Acids', 'Counts', 'Frequencies', 'Percentages']))#[i.split(':') for i in counts[int(len(counts)/2):]]))
            with count_col2:
                lit.write("For the amino acids present in input sequence")
                lit.table(pd.DataFrame(nz_cnt_freq, columns=['Amino Acids', 'Counts', 'Frequencies', 'Percentages']))
                
            lit.markdown('''<br>''', unsafe_allow_html=True)

            
            
            
##            lit.write('_Amino acid frequencies_: '); 
##            count_col1, count_col2 = lit.columns(2)
##            with count_col1:
##                lit.write("For all amino acids")
##                
##                lit.table(pd.DataFrame([i.split(':') for i in counts], columns=['Amino Acids', 'Counts']))
##            with count_col2:
##                lit.write("For the amino acids present in input sequence")
##                
##                lit.table(pd.DataFrame([i.split(':') for i in scount], columns=['Amino Acids', 'Counts']))
                
            most_common = ', '.join([i+': '+str(count[i]) for i in count.keys() if count[i]==max(count.values())])
            least_common = ', '.join([i+': '+str(count[i]) for i in count.keys() if count[i]==min([i for i in count.values() if i!=0])])
            not_present = ", ".join([i.split(':')[0] for i in counts if int(i.split(':')[1]) == 0])
            phiAA = str(sum([count[i] for i in count.keys() if i in list('RNDCQEHKSTY')]))
            phoAA = str(sum([count[i] for i in count.keys() if i in list('GAMLIVFWP')]))
            basicAA = str(sum([count[i] for i in count.keys() if i in list('HRK')]))
            acidicAA = str(sum([count[i] for i in count.keys() if i in list('DE')]))
            sec_str_frac =  ', '.join(str(round(i,3)) for i in propar.secondary_structure_fraction())
            f.write('Most Common Residue: '+most_common+'\n')
            f.write('Least Common Residue: '+least_common+'\n')
            f.write('Residues not present in the sequence: '+not_present+'\n')
            f.write("No. of Hydrophilic residues: "+phiAA+'\n')
            f.write("No. of Hydrophobic residues: "+phoAA+'\n')
            f.write("No. of Basic residues: "+basicAA+'\n')
            f.write("No. of Acidic residues: "+acidicAA+'\n')
            f.write('Secondary Structure Fraction (Helix, Turn, Sheet): '+sec_str_frac+'\n')

            lit.write('Protein Basic Information Table: ');
            lit.table(
                pd.DataFrame(
                    [
                        ["Most common residue:", most_common],
                        ["Least common residue:", least_common],
                        ["Residues not present in the sequence:", not_present],
                        ["No. of Hydrophilic residues:", phiAA],
                        ["No. of Hydrophobic residues:", phoAA],
                        ["No. of Basic residues:", basicAA],
                        ["No. of Acidic residues:", acidicAA],
                        ["Secondary Structure Fraction (Helix, Turn, Sheet):", sec_str_frac]
                    ],
                    
                    columns=['Data', 'Values']
                )
            )

            lit.markdown('''<br>''', unsafe_allow_html=True)
            


            lit.markdown('''<br>''', unsafe_allow_html=True)
#            lit.write('_Secondary Structure Fraction (Helix, Turn, Sheet)_:      '+', '.join(str(round(i,3)) for i in propar.secondary_structure_fraction()));
            



        if my_input and physicochemical:

            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
        
            lit.subheader('_Physicochemical Features_: '); f.write('\n-->Physicochemical Properties:\n')
            pep = pept(my_input)
            propar = proparm.ProteinAnalysis(my_input)
            #Peptide feature calculation
            f.write('Aliphatic Index: '+str(pep.aliphatic_index())+'\n')
            f.write('Instability Index: '+str(pep.instability_index())+'\n')
            f.write('Hydrophobicity: '+str(pep.hydrophobicity())+'\n')
            f.write('Hydrophobic Moment: '+str(pep.hydrophobic_moment())+'\n')
            f.write('Isoelectric Point: '+str(pep.isoelectric_point())+'\n')
            f.write('Molecular Weight: '+str(pep.molecular_weight())+'\n')
            f.write('Charge (at pH 7): '+str(pep.charge())+'\n')
            f.write('Aromaticity: '+str(propar.aromaticity())+'\n')
            f.write('Molar Extinction Coefficient (Cysteine|Cystine):   '+str(propar.molar_extinction_coefficient()).replace('(','').replace(')','')+'\n')
            try:
                flexibility = ['Flexibility: ', '\n'.join([str(i) for i in propar.flexibility()])]; f.write('Flexibility: '+str(propar.flexibility())+'\n')
            except:
                flexibility = ['Flexibility: ','Cannot be computed for peptide with non-standard amino acid residues']; f.write('Flexibility: Cannot be computed for peptide with non-standard amino acid residues\n')

            lit.table(
                pd.DataFrame(
                    [
                        ['Aliphatic Index:', str(pep.aliphatic_index())],
                        ['Instability Index: ', str(pep.instability_index())],
                        ['Hydrophobicity: ', str(pep.hydrophobicity())],
                        ['Hydrophobic Moment: ', str(pep.hydrophobic_moment())],
                        ['Isoelectric Point: ', str(pep.isoelectric_point())],
                        ['Molecular Weight: ', str(pep.molecular_weight())],
                        ['Charge (at pH 7): ', str(pep.charge())],
                        ['Aromaticity: ', str(propar.aromaticity())],
                        ['Molar Extinction Coefficient (Cysteine|Cystine):   ', str(propar.molar_extinction_coefficient()).replace('(','').replace(')','')],
                        flexibility
                    ],
                    
                    columns=['Data', 'Values']
                )
            )
            
            

        
        if my_input and ctddescriptors:

            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)

            lit.subheader('_CTD Descriptors_: '); f.write('\n-->CTD Descriptors: \n')
            dictC = calC(my_input)
            dictT = calT(my_input)
            dictD = calD(my_input)
            colC, spacerCol, colT = lit.columns([2, 1, 2])
            with colC:
                lit.write("**Composition Descriptors**"); f.write("Composition Descriptors\n")
                for i in dictC.keys():
                    f.write(i.replace('_', '')+': '+str(dictC[i])+'\n')

                lit.table(
                    pd.DataFrame(
                        [[i.replace('_', ''), dictC[i]] for i in dictC.keys()],
                        columns=['Descriptor', 'Value']
                        )
                    )

            with spacerCol:
                _ = [lit.write('\t') for i in range(20)]
            
            with colT:
                lit.write("**Transition Descriptors**"); f.write("Transition Descriptors\n")
                for i in dictT.keys():
                    f.write(i.replace('_', '')+': '+str(dictT[i])+'\n')

                lit.table(
                    pd.DataFrame(
                        [[i.replace('_', ''), dictT[i]] for i in dictT.keys()],
                        columns=['Descriptor', 'Value']
                        )
                    )
                
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.write("**Distribution Descriptors**"); f.write("Distribution Descriptors\n")
            scolD1, scolD2, scolD3 = lit.columns(3)
            listD = [k+': '+str(dictD[k]) for k in dictD.keys()]
            _ = [f.write(i.replace('_', '')+'\n') for i in listD]
            
            with scolD1:
                lit.table(
                    pd.DataFrame(
                        [i.replace('_', '').split(':') for i in listD[:35]],
                        columns=['Descriptors', 'Values']
                        )
                    )
                
            with scolD2:
                lit.table(
                    pd.DataFrame(
                        [i.replace('_', '').split(':') for i in listD[35:70]],
                        columns=['Descriptors', 'Values']
                        )
                    )
                
            with scolD3:
                lit.table(
                    pd.DataFrame(
                        [i.replace('_', '').split(':') for i in listD[70:]],
                        columns=['Descriptors', 'Values']
                        )
                    )
                




        if my_input and others:

            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
        
            lit.subheader('_Other Descriptors_: '); f.write('\n-->Other Descriptors: \n')
            pep = pept(my_input)
            dict1 = pep.descriptors()
            desc = [k+': '+str(dict1[k]) for k in dict1.keys()]; _ = [f.write(k+': '+str(dict1[k])+'\n') for k in dict1.keys()]
            dcol_list = lit.columns(4)
            loop = 0
            for i in dcol_list:
                myList = []
                with i:
                    for j in range(loop, len(desc)):
                        if j%18==0 and j!=0:
                            myList.append(desc[j])
                            loop = j+1
                            break
                        else:
                            myList.append(desc[j])

                    lit.table(pd.DataFrame([i.split(':') for i in myList], columns=['Descriptor', 'Value']))


        f.close()
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.download_button("Download output file", open('descriptor_calculation_output.txt'), file_name='calculation_output.txt')

elif not my_input and submit:
    lit.error("Please enter the input!")

    
print('Data processing/computation complete')

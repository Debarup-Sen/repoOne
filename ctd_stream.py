#OM Sai Ram
#Sai
#Sai
#Feature calculation app
import streamlit as lit
from peptides import Peptide as pept
import Bio.SeqUtils.ProtParam as proparm
from propy.CTD import CalculateC as calC, CalculateT as calT, CalculateD as calD

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
                        my_input = j[6]
                        break

        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        f = open('descriptor_calculation_output.txt', 'w')
        if my_input and composition:
            lit.subheader('_Composition Features_: '); f.write("-->Composition features: \n")
            pep = pept(my_input)
            propar = proparm.ProteinAnalysis(my_input)
            count = pep.counts()
            counts = ['**'+k+'**:  '+str(count[k]) for k in count.keys()]
            scount = [i for i in counts if '0' not in i]
            lit.write('_Amino acid counts_: '); f.write('Amino acid counts:\n')
            count_col1, count_col2 = lit.columns(2)
            with count_col1:
                lit.write("For all amino acids"); f.write('For all amino acids:\n')
                _ = [f.write(i+'\n') for i in counts]
                scol1, scol2 = lit.columns(2)
                with scol1:
                    _ = [lit.write(i) for i in counts[:int(len(counts)/2)]]
                with scol2:
                    _ = [lit.write(i) for i in counts[int(len(counts)/2):]]

            with count_col2:
                lit.write("For the amino acids present in input sequence"); f.write("For the amino acids present in input sequence:\n")
                _ = [f.write(i+'\n') for i in scount]
                scol3, scol4 = lit.columns(2)
                with scol3:
                    _ = [lit.write(i) for i in scount[:int(len(scount)/2)]]
                with scol4:
                    _ = [lit.write(i) for i in scount[int(len(scount)/2):]]
            lit.markdown('''<br>''', unsafe_allow_html=True)
            most_common = [i+': '+str(count[i]) for i in count.keys() if count[i]==max(count.values())][0]
            least_common = [i+': '+str(count[i]) for i in count.keys() if count[i]==min([i for i in count.values() if i!=0])][0]
            not_present = ", ".join([i[2] for i in counts if '0' in i])
            phiAA = str(sum([count[i] for i in count.keys() if i in list('RNDCQEHKSTY')]))
            phoAA = str(sum([count[i] for i in count.keys() if i in list('GAMLIVFWP')]))
            basicAA = str(sum([count[i] for i in count.keys() if i in list('HRK')]))
            acidicAA = str(sum([count[i] for i in count.keys() if i in list('DE')]))
            lit.write("Most common residue: "+most_common); f.write('Most Common Residue: '+most_common+'\n')
            lit.write("Least common residue: "+least_common); f.write('Least Common Residue: '+least_common+'\n')
            lit.write("Residues not present in the sequence: "+not_present); f.write('Residues not present in the sequence: '+not_present+'\n')
            lit.write("No. of Hydrophilic residues: "+phiAA); f.write("No. of Hydrophilic residues: "+phiAA+'\n')
            lit.write("No. of Hydrophobic residues: "+phoAA); f.write("No. of Hydrophobic residues: "+phoAA+'\n')
            lit.write("No. of Basic residues: "+basicAA); f.write("No. of Basic residues: "+basicAA+'\n')
            lit.write("No. of Acidic residues: "+acidicAA); f.write("No. of Acidic residues: "+acidicAA+'\n')

            lit.markdown('''<br>''', unsafe_allow_html=True)
            
            count = pep.frequencies()
            counts = ['**'+k+'**:  '+str(round(count[k],3)) for k in count.keys()]
            scount = [i for i in counts if '0.0' != i[-3:]]
            lit.write('_Amino acid frequencies_: '); f.write('Amino Acid Frequencies: \n')
            count_col1, count_col2 = lit.columns(2)
            with count_col1:
                lit.write("For all amino acids"); f.write("For all amino acids:\n")
                _ = [f.write(i+'\n') for i in counts]
                scol1, scol2 = lit.columns(2)
                with scol1:
                    _ = [lit.write(i) for i in counts[:int(len(counts)/2)]]
                with scol2:
                    _ = [lit.write(i) for i in counts[int(len(counts)/2):]]

            with count_col2:
                lit.write("For the amino acids present in input sequence"); f.write("For the amino acids present in input sequence:\n")
                _ = [f.write(i+'\n') for i in scount]
                scol3, scol4 = lit.columns(2)
                with scol3:
                    _ = [lit.write(i) for i in scount[:int(len(count)/2)]]
                with scol4:
                    _ = [lit.write(i) for i in scount[int(len(count)/2):]]

            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.write('_Secondary Structure Fraction (Helix, Turn, Sheet)_:      '+', '.join(str(round(i,3)) for i in propar.secondary_structure_fraction()));
            f.write('Secondary Structure Fraction (Helix, Turn, Sheet): '+', '.join(str(round(i,3)) for i in propar.secondary_structure_fraction())+'\n')


        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        if my_input and physicochemical:
            lit.subheader('_Physicochemical Features_: '); f.write('\n-->Physicochemical Properties:\n')
            pep = pept(my_input)
            propar = proparm.ProteinAnalysis(my_input)
            #Peptide feature calculation
            lit.write('Aliphatic Index: '+str(pep.aliphatic_index())); f.write('Aliphatic Index: '+str(pep.aliphatic_index())+'\n')
            lit.write('Instability Index: '+str(pep.instability_index())); f.write('Instability Index: '+str(pep.instability_index())+'\n')
            lit.write('Hydrophobicity: '+str(pep.hydrophobicity())); f.write('Hydrophobicity: '+str(pep.hydrophobicity())+'\n')
            lit.write('Hydrophobic Moment: '+str(pep.hydrophobic_moment())); f.write('Hydrophobic Moment: '+str(pep.hydrophobic_moment())+'\n')
            lit.write('Isoelectric Point: '+str(pep.isoelectric_point())); f.write('Isoelectric Point: '+str(pep.isoelectric_point())+'\n')
            lit.write('Molecular Weight: '+str(pep.molecular_weight())); f.write('Molecular Weight: '+str(pep.molecular_weight())+'\n')
            lit.write('Charge (at pH 7): '+str(pep.charge())); f.write('Charge (at pH 7): '+str(pep.charge())+'\n')
            lit.write('Aromaticity: '+str(propar.aromaticity())); f.write('Aromaticity: '+str(propar.aromaticity())+'\n')
            lit.write('Molar Extinction Coefficient (Cysteine|Cystine):   '+str(propar.molar_extinction_coefficient()).replace('(','').replace(')',''));
            f.write('Molar Extinction Coefficient (Cysteine|Cystine):   '+str(propar.molar_extinction_coefficient()).replace('(','').replace(')','')+'\n')
            try:
                lit.write('Flexibility: '+str(propar.flexibility())); f.write('Flexibility: '+str(propar.flexibility())+'\n')
            except:
                lit.write('Flexibility: Cannot be computed for peptide with non-standard amino acid residues'); f.write('Flexibility: Cannot be computed for peptide with non-standard amino acid residues\n')

        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        if my_input and ctddescriptors:
            lit.subheader('_CTD Descriptors_: '); f.write('\n-->CTD Descriptors: \n')
            dictC = calC(my_input)
            dictT = calT(my_input)
            dictD = calD(my_input)
            colC, colT = lit.columns(2)
            with colC:
                lit.write("**Composition Descriptors**"); f.write("Composition Descriptors\n")
                for i in dictC.keys():
                    lit.write(i.replace('_', '')+': '+str(dictC[i])); f.write(i.replace('_', '')+': '+str(dictC[i])+'\n')
            with colT:
                lit.write("**Transition Descriptors**"); f.write("Transition Descriptors\n")
                for i in dictT.keys():
                    lit.write(i.replace('_', '')+': '+str(dictT[i])); f.write(i.replace('_', '')+': '+str(dictT[i])+'\n')
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.write("**Distribution Descriptors**"); f.write("Distribution Descriptors\n")
            scolD1, scolD2, scolD3, scolD4 = lit.columns(4)
            listD = [k+': '+str(dictD[k]) for k in dictD.keys()]; _ = [f.write(i.replace('_', '')+'\n') for i in listD]
            with scolD1:
                _ = [lit.write(i.replace('_', '')) for i in listD[:27]]
            with scolD2:
                _ = [lit.write(i.replace('_', '')) for i in listD[27:53]]
            with scolD3:
                _ = [lit.write(i.replace('_', '')) for i in listD[53:79]]
            with scolD4:
                _ = [lit.write(i.replace('_', '')) for i in listD[79:]]

        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        if my_input and others:
            lit.subheader('_Other Descriptors_: '); f.write('\n-->Other Descriptors: \n')
            pep = pept(my_input)
            dict1 = pep.descriptors()
            desc = [k+': '+str(dict1[k]) for k in dict1.keys()]; _ = [f.write(k+': '+str(dict1[k])+'\n') for k in dict1.keys()]
            dcol_list = lit.columns(4)
            loop = 0
            for i in dcol_list:
                with i:
                    for j in range(loop, len(desc)):
                        if j%19==0 and j!=0:
                            lit.write(desc[j])
                            loop = j+1
                            break
                        else:
                            lit.write(desc[j])
        f.close()
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.markdown('''<br>''', unsafe_allow_html=True)
        lit.download_button("Download output file", open('descriptor_calculation_output.txt'), file_name='calculation_output.txt')

elif not my_input and submit:
    lit.error("Please enter the input!")

    
print('Data processing/computation complete')

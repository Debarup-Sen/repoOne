#Om Sai Ram
#Sai

import streamlit as lit
import subprocess as proc
from skbio.alignment import local_pairwise_align_protein as lalign, global_pairwise_align_protein as galign
from skbio import Protein
import pandas as pd
from io import StringIO


lit.set_page_config(layout='wide')

lit.write("""
# Welcome to the AMPDB Sequence Alignment Toolbox!
*A toolbox for all alignment purposes.*

""")
tool = lit.radio(
    "Select an alignment option:",
    ('BLASTp (Basic Local Alignment Search against AMPDB)', 'MUSCLE (Multiple Sequence Alignment)', 'Needleman-Wunsch (Global Pairwise Alignment)', 'Smith-Waterman (Local Pairwise Alignment)'))

if 'BLASTp' in tool:
    query = lit.text_area('Enter your input protein sequence (in FASTA/multi-FASTA format/plain text sequence format/AMPDB Acc. ID, e.g. AMPDB_111) here').upper()
    file_query = lit.file_uploader("Or, you may upload file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    outfmt = lit.radio(
        "Select an output format:",
        ('Default', 'Pairwise', 'Query-anchored showing identities', 
     'Query-anchored no identities', 'Flat query-anchored showing identities',
     'Flat query-anchored no identities', 'BLAST XML', 'Tabular', 'Tabular with comment lines',
     'Seqalign (Text ASN.1)', 'Seqalign (Binary ASN.1)', 'Comma-separated values',
     'BLAST archive (ASN.1)', 'Seqalign (JSON)')
        )
    outfmt = ('0' if 'Pairwise' in outfmt
              else '1' if 'Query-anchored showing identities' in outfmt
              else '2' if 'Query-anchored no identities' in outfmt
              else '3' if 'Flat query-anchored showing identities' in outfmt
              else '4' if 'Flat query-anchored no identities' in outfmt
              else '5' if 'BLAST XML' in outfmt
              else '6' if outfmt=='7) Tabular'
              else '7' if 'Tabular with comment lines' in outfmt
              else '8' if 'Seqalign (Text ASN.1)' in outfmt
              else '9' if 'Seqalign (Binary ASN.1)' in outfmt
              else '10' if 'Comma-separated values' in outfmt
              else '11' if 'BLAST archive (ASN.1)' in outfmt
              else '12' if 'Seqalign (JSON)' in outfmt
              else 'def' if 'Default' in outfmt
              else None)
    submit = lit.button('Submit')
    if file_query:
        query = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMPDB_' in query:
        with open('master_dataset.tsv') as f:
            l = ' '
            while(True):
                i = f.readline()
                if i=='':
                    lit.error('The AMPDB Acc. ID does not match with our database. Please re-check')
                    query = None
                    break
                j = i.split('\t')
                if query in j[1]:
                    query = '>'+query+'\n'+j[3]
                    break                
    if query and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        open('blast_input.txt', 'w').write(query)
        if outfmt == 'def':
            proc.run(('blastp -query blast_input.txt -db ampdb -out blast_output_def1 -outfmt 0').split())
            proc.run(('blastp -query blast_input.txt -db ampdb -out blast_output_def2 -outfmt 7').split())
        else:
            proc.run(('blastp -query blast_input.txt -db ampdb -out blast_output -outfmt '+outfmt).split())
            
        lit.info("Your output below: [Formats 7-13 show no output when no hits are found]")

        if outfmt == 'def':
            myFile = [i for i in open('blast_output_def2').readlines()]
            lit.text(''.join([i for i in myFile[:5] if '# Fields: ' not in i]))
            try:
                headers = [i for i in myFile if '# Fields: ' in i][0].replace('# Fields: ','').split(',')
                data = [i.strip().split('\t') for i in myFile if '#' not in i]
                myDF = pd.DataFrame(data, columns=headers)
                del myFile, data, headers
                lit.table(myDF)
                lit.markdown('<br><br><u><b>*Alignments:*</b></u><br>', unsafe_allow_html=True)
                myDF.to_csv('blast_output_def2', sep='\t')

                output1 = ''.join(open('blast_output_def1').readlines()[18:])
                crsr = 0
                for i in range(len(output1)):
                    if '>' in output1[i]:
                        crsr = i
                        break
                output1 = ''.join(output1[crsr:])
                lit.text(output1)
                open('blast_output_def1', 'w').write(output1)

                open('blast_output', 'w').write(open('blast_output_def2').read() +'\n\nAlignments:\n'+ open('blast_output_def1').read())
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')
            except:
                lit.text("No hits found")

        elif outfmt == '0':
            myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '1':
            myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '2':
            myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '3':
            myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '4':
            myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '5':
            myFile = ''.join(open('blast_output').readlines()[18:])
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.xml')
            
        elif outfmt == '6' or outfmt == '10':
            lit.text('(Please choose "Tabular with comment lines" to see column headers)')
            myFile = pd.DataFrame([(i.strip().split(',') if ',' in i else i.strip().split('\t')) for i in open('blast_output').readlines()])
            lit.table(myFile)
            if outfmt == '6':
                myFile.to_csv('blast_output', sep='\t')
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.tsv')
            elif outfmt == '10':
                myFile.to_csv('blast_output')
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.csv')

        elif outfmt == '7':
            myFile = [i for i in open('blast_output').readlines()]
            headers = [i for i in myFile if '# Fields: ' in i][0].replace('# Fields: ','').split(',')
            data = [i.strip().split('\t') for i in myFile if '#' not in i]
            myDF = pd.DataFrame(data, columns=headers)
            del myFile, data, headers
            lit.table(myDF)
            myDF.to_csv('blast_output', sep='\t')
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.csv')

        elif outfmt == '8':
            myFile = ''.join(open('blast_output').readlines())
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

        elif outfmt == '9':
            lit.text('Binary output cannot be displayed in browser. Please download file to view output')
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')

        elif outfmt == '11':
            myFile = ''.join(open('blast_output').readlines())
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.asn')

        elif outfmt == '12':
            myFile = ''.join(open('blast_output').readlines())
            lit.text(myFile)
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.json')
##        #########################################################
##        if  outfmt=='9':
##            lit.download_button("Download output file", open('blast_output', 'rb'), file_name='BLAST_out')
##        else:
##            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')
##        #########################################################
    elif submit and not query:
        lit.error("Please enter input sequence!")

if 'MUSCLE' in tool:
    multiseq = lit.text_area('Enter your input protein sequences (in FASTA format/multi-FASTA format/AMPDB Acc. ID, e.g. AMPDB_111 one in each line) here:').upper()
    file_query = lit.file_uploader("Or, you may upload file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    submit = lit.button('Submit')
    if file_query:
        multiseq = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if multiseq and submit:
        if len([i for i in multiseq.split('\n') if i!=''])>=1 and 'AMPDB_' in multiseq:
            multiseq = [i for i in multiseq.replace(' ', '').split('\n') if i!='']
            with open('master_dataset.tsv') as f, open(r'mafft_input.txt', 'w') as g:
                for k in multiseq:
                    f.seek(0,0)
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMPDB Acc. ID {k} does not match with our database. Please re-check')
                            multiseq = None
                            break
                        j = i.split('\t')
                        if k in j[1]:
                            g.write('>'+k+'\n'+j[3]+'\n')
                            break
        else:
            open('mafft_input.txt', 'w').write(multiseq)
    
    if multiseq and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        proc.run('muscle -in mafft_input.txt -out mafft_output'.split())
        lit.info("Your output below (In case you do not see any output, please re-check your input for invalid characters or non-standard residues):")
        lit.text(open(r'mafft_output').read())
        lit.download_button("Download output file", open('mafft_output'), file_name='MAFFT_out')
    elif submit and not multiseq:
        lit.error("Please enter input sequence!")


if 'Needleman-Wunsch' in tool:
    lit.text("FASTA format, plain text sequence format supported.")
    query = myquery = lit.text_area('Enter your query protein sequence here:').upper()
    file_query = lit.file_uploader("Or, you may upload query file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    subject = mysubject = lit.text_area('Enter your subject protein sequence here:').upper()
    file_subject = lit.file_uploader("Or, you may upload subject file")#, label_visibility="collapsed")
    submit = lit.button('Submit')
    if file_query:
        query = myquery = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if file_subject:
        subject = mysubject = StringIO(file_subject.getvalue().decode("utf-8")).read().upper()
    if query and subject and submit:
        if '>' in query:
            query = ''.join(query.split('\n')[1:])
        if '>' in subject:
            subject = ''.join(subject.split('\n')[1:])
        if '\n' in query:
            query = query.replace('\n', '')
        if '\n' in subject:
            subject = subject.replace('\n', '')
        if 'AMPDB_' not in query and query.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            query = None
        elif 'AMPDB_' in query and query.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            query = None
        if 'AMPDB_' not in subject and subject.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            subject = None
        elif 'AMPDB_' in subject and subject.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            subject = None
        if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMPDB_' in query:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMPDB Acc. ID {query} does not match with our database. Please re-check')
                        query = None
                        break
                    j = i.split('\t')
                    if query in j[1]:
                        query = j[3]
                        break
                    
        if subject and len([i for i in subject.split('\n') if i!=''])==1 and 'AMPDB_' in subject:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMPDB Acc. ID {subject} does not match with our database. Please re-check')
                        subject = None
                        break
                    j = i.split('\t')
                    if subject in j[1]:
                        subject = j[3]
                        break
        if query and subject:
            lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
            alignment, score, start_end_positions = galign(Protein(query.strip()), Protein(subject.strip()))
            lit.info("Your output below:")
            lit.text(alignment)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            alignment.write(open('NWFile', 'w'))
            lines = [i for i in open('NWFile').readlines() if i!='']
            lit.text("Full alignment:")
            myquery = '>'+myquery+'\n'+lines[1]
            mysubject = '>'+mysubject+'\n'+lines[3]
            lit.text(myquery)
            lit.text(mysubject)
            lit.text("Score: "+str(score))
            my_write_file = 'AMPDB Needleman-Wunsch Output:\n\nAlignment\n'+myquery+'\n'+mysubject+"\nScore: "+str(score)+'\n'
            open('NWFile', 'w').write(my_write_file)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.download_button("Download output file", open('NWFile'), file_name='NW_out')
    elif submit and (not query or not subject):
        lit.error("Please enter input sequence!")


if 'Smith-Waterman' in tool:
    lit.text("FASTA format, plain text sequence format supported.")
    query = myquery = lit.text_area('Enter your query protein sequence here').upper()
    file_query = lit.file_uploader("Or, you may upload query file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    subject = mysubject = lit.text_area('Enter your subject protein sequence here').upper()
    file_subject = lit.file_uploader("Or, you may upload subject file")#, label_visibility="collapsed")
    submit = lit.button('Submit')
    if file_query:
        query = myquery = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if file_subject:
        subject = mysubject = StringIO(file_subject.getvalue().decode("utf-8")).read().upper()
    if query and subject and submit:
        if '>' in query:
            query = ''.join(query.split('\n')[1:])
        if '>' in subject:
            subject = ''.join(subject.split('\n')[1:])
        if '\n' in query:
            query = query.replace('\n', '')
        if '\n' in subject:
            subject = subject.replace('\n', '')
        if 'AMPDB_' not in query and query.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            query = None
        elif 'AMPDB_' in query and query.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            query = None
        if 'AMPDB_' not in subject and subject.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            subject = None
        elif 'AMPDB_' in subject and subject.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            subject = None
        if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMPDB_' in query:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMPDB Acc. ID {query} does not match with our database. Please re-check')
                        query = None
                        break
                    j = i.split('\t')
                    if query in j[1]:
                        query = j[3]
                        break
        if subject and len([i for i in subject.split('\n') if i!=''])==1 and 'AMPDB_' in subject:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMPDB Acc. ID {subject} does not match with our database. Please re-check')
                        subject = None
                        break
                    j = i.split('\t')
                    if subject in j[1]:
                        subject = j[3]
                        break
        if query and subject:
            lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
            alignment, score, start_end_positions = lalign(Protein(query.strip()), Protein(subject.strip()))
            lit.info("Your output below:")
            lit.text(alignment)
            alignment.write(open('SWFile', 'w'))
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            alignment.write(open('SWFile', 'w'))
            lines = [i for i in open('SWFile').readlines() if i!='']
            lit.text("Full alignment:")
            myquery = '>'+myquery+'\n'+lines[1]
            mysubject = '>'+mysubject+'\n'+lines[3]
            lit.text(myquery)
            lit.text(mysubject)
            lit.text("Score: "+str(score))
            my_write_file = 'AMPDB Smith-Waterman Output:\n\nAlignment\n'+myquery+'\n'+mysubject+"\nScore: "+str(score)+'\n'
            open('SWFile', 'w').write(my_write_file)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.download_button("Download output file", open('SWFile'), file_name='SW_out')
    elif submit and not query:
        lit.error("Please enter input sequence!")

        

lit.write("*Thank you!*")

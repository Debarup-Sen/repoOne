#Om Sai Ram
#Sai

import streamlit as lit
import subprocess as proc
from skbio.alignment import local_pairwise_align_protein as lalign, global_pairwise_align_protein as galign
from skbio import Protein

lit.set_page_config(layout='wide')

lit.write("""
# Welcome to the AMPDB Sequence Alignment Toolbox!
*A toolbox for all your alignment needs.*

""")
tool = lit.radio(
    "What\'s your alignment choice?",
    ('BLAST', 'MUSCLE', 'Needleman-Wunsch', 'Smith-Waterman'))

if tool == 'BLAST':
    query = lit.text_area('Enter your input sequence (in FASTA or plain text sequence format)/AMPDB Acc. ID here').upper()
    lit.markdown('<br>', unsafe_allow_html=True)
    outfmt = lit.radio(
        "What is your output format preference? Available formats:",
        ('1) Pairwise', '2) Query-anchored showing identities', 
     '3) Query-anchored no identities', '4) Flat query-anchored showing identities',
     '5) Flat query-anchored no identities', '6) BLAST XML', '7) Tabular', '8) Tabular with comment lines',
     '9) Seqalign (Text ASN.1)', '10) Seqalign (Binary ASN.1)', '11) Comma-separated values',
     '12) BLAST archive (ASN.1)', '13) Seqalign (JSON)')
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
              else None)
    submit = lit.button('Submit')
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
        proc.run(('blastp -query blast_input.txt -db ampdb -out blast_output -outfmt '+outfmt).split())
        lit.info("Your output below: [Formats 7-13 show no output when no hits are found]")
        if outfmt=='6' or outfmt=='10':
            lit.text('(Please choose "Tabular with comment lines" to see column headers)')
            lit.text(''.join((open('blast_output').readlines())))
        elif outfmt=='7':
            lit.text(''.join((open('blast_output').readlines()[1:])))
        elif outfmt=='9':
            lit.text('Binary output cannot be displayed in browser. Please download file to view output')
        else:
            lit.text(''.join((open('blast_output').readlines()[18:])))
        
        if  outfmt=='9':
            lit.download_button("Download output file", open('blast_output', 'rb'), file_name='BLAST_out')
        else:
            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')
    elif submit and not query:
        lit.error("Please enter input sequence!")

if tool == 'MUSCLE':
    multiseq = lit.text_area('Enter your input sequences (in FASTA format)/AMPDB Acc. IDs (one in each line) here:').upper()
    lit.markdown('<br>', unsafe_allow_html=True)
    submit = lit.button('Submit')
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


if tool == 'Needleman-Wunsch':
    lit.text("FASTA format, plain text sequence format supported.")
    query = myquery = lit.text_area('Enter your query sequence here:').upper()
    subject = mysubject = lit.text_area('Enter your subject sequence here:').upper()
    submit = lit.button('Submit')
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
            mysubject = '>'+mysubject+'\n'+lines[2]
            lit.text(myquery)
            lit.text(mysubject)
            #lit.write(open('NWFile').read())
            lit.text("Score: "+str(score))
            open('NWFile', 'w').write('AMPDB Needleman-Wunsch Output:\n\nAlignment\n')
            alignment.write(open('NWFile', 'a'))
            open('NWFile', 'a').write("Score: "+str(score)+'\n')
            lit.download_button("Download output file", open('NWFile'), file_name='NW_out')
    elif submit and (not query or not subject):
        lit.error("Please enter input sequence!")


if tool == 'Smith-Waterman':
    lit.text("FASTA format, plain text sequence format supported.")
    query = lit.text_area('Enter your query sequence here').upper()
    subject = lit.text_area('Enter your subject sequence here').upper()
    submit = lit.button('Submit')
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
            alignment, score, start_end_positions = lalign(Protein(query), Protein(subject))
            lit.info("Your output below:")
            lit.text(alignment)
            alignment.write(open('SWFile', 'w'))
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.text("Full alignment:")
            lit.write(open('SWFile').read())
            lit.text("Score: "+str(score))
            open('SWFile', 'w').write('AMPDB Smith-Waterman Output:\n\nAlignment\n')
            alignment.write(open('SWFile', 'a'))
            open('SWFile', 'a').write("Score: "+str(score)+'\n')
            lit.download_button("Download output file", open('SWFile'), file_name='SW_out')
    elif submit and not query:
        lit.error("Please enter input sequence!")

        

lit.write("*Thank you!*")

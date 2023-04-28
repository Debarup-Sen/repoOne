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
    ('BLAST', 'MAFFT', 'Needleman-Wunsch', 'Smith-Waterman'))

if tool == 'BLAST':
    query = lit.text_area('Enter your input sequence/AMPDB Acc. ID here')
    lit.markdown('<br>', unsafe_allow_html=True)
    outfmt = lit.radio(
        "What is your output format preference? Available formats:",
        ('Pairwise', 'Query-anchored showing identities', 
     'Query-anchored no identities', 'Flat query-anchored showing identities',
     'Flat query-anchored no identities', 'BLAST XML', 'Tabular', 'Tabular with comment lines',
     'Seqalign (Text ASN.1)', 'Seqalign (Binary ASN.1)', 'Comma-separated values',
     'BLAST archive (ASN.1)', 'Seqalign (JSON)')
        )
    outfmt = ('0' if outfmt=='Pairwise'
              else '1' if outfmt=='Query-anchored showing identities'
              else '2' if outfmt=='Query-anchored no identities'
              else '3' if outfmt=='Flat query-anchored showing identities'
              else '4' if outfmt=='Flat query-anchored no identities'
              else '5' if outfmt=='BLAST XML'
              else '6' if outfmt=='Tabular'
              else '7' if outfmt=='Tabular with comment lines'
              else '8' if outfmt=='Seqalign (Text ASN.1)'
              else '9' if outfmt=='Seqalign (Binary ASN.1)'
              else '10' if outfmt=='Comma-separated values'
              else '11' if outfmt=='BLAST archive (ASN.1)'
              else '12' if outfmt=='Seqalign (JSON)'
              else None)
    submit = lit.button('Submit')
    if query and 'AMPDB_' in query:
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
                    query = '>'+query+'\n'+j[6]
                    break                
    if query and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        open('blast_input.txt', 'w').write(query)
        proc.run(('blastp -query blast_input.txt -db ampdb -out blast_output -outfmt '+outfmt).split())
        lit.info("Your output below:")
        if outfmt=='6' or outfmt=='10':
            lit.text('(Please choose "Tabular with comment lines" to see column headers')
            lit.text(''.join((open('blast_output').readlines())))
        else:
            lit.text(''.join((open('blast_output').readlines()[18:])))
        lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')
    elif submit and not query:
        lit.error("Please enter input sequence!")

if tool == 'MAFFT':
    multiseq = lit.text_area('Enter your input sequences (in FASTA format)/AMPDB Acc. IDs (one in each line) here:')
    lit.markdown('<br>', unsafe_allow_html=True)
    submit = lit.button('Submit')
    if multiseq and submit:
        if 'AMPDB_' in multiseq:
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
                            g.write('>'+k+'\n'+j[6]+'\n')
                            break
        else:
            open('mafft_input.txt', 'w').write(multiseq)
    if multiseq and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        proc.run('mafft --auto --clustalout --reorder --thread -1 mafft_input.txt')
        lit.info("Your output below (In case you do not see any output, please re-check your input for invalid characters or non-standard residues):")
        lit.text(open(r'mafft_output').read())
        lit.download_button("Download output file", open('mafft_output'), file_name='MAFFT_out')
    elif submit and not multiseq:
        lit.error("Please enter input sequence!")


if tool == 'Needleman-Wunsch':
    lit.text("FASTA format, plain text sequence format supported.")
    query = lit.text_area('Enter your query sequence here:')
    subject = lit.text_area('Enter your subject sequence here:')
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
        if query and 'AMPDB_' in query:
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
                        query = j[6]
                        break
                    
        if subject and 'AMPDB_' in subject:
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
                        subject = j[24]
                        break
        if query and subject:
            lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
            alignment, score, start_end_positions = galign(Protein(query), Protein(subject))
            lit.info("Your output below:")
            lit.text(alignment)
            alignment.write(open('NWFile', 'w'))
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.text("Full alignment:")
            lit.write(open('NWFile').read())
            lit.text("Score: "+str(score))
            open('NWFile', 'w').write('AMPDB Needleman-Wunsch Output:\n\nAlignment\n')
            alignment.write(open('NWFile', 'a'))
            open('NWFile', 'a').write("Score: "+str(score)+'\n')
            lit.download_button("Download output file", open('NWFile'), file_name='NW_out')
    elif submit and (not query or not subject):
        lit.error("Please enter input sequence!")


if tool == 'Smith-Waterman':
    query = lit.text_area('Enter your query sequence here')
    subject = lit.text_area('Enter your subject sequence here')
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
        if query and 'AMPDB_' in query:
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
                        query = j[6]
                        break
        if subject and 'AMPDB_' in subject:
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
                        subject = j[24]
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

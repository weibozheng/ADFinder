import sys
import re
import pandas
from Bio import SeqIO
def pointer_check(t,seqr,tag,k_start,k_end):
	if t=='deletion':
		start=int(k_start)
		end=int(k_end)
		len_pointer=end-start+1
		left_check=0
		left_str1=''
		left_str2=''
		right_check=0
		right_str1=''
		right_str2=''
		cost_l=0
		cost_r=0
		for i in range(15):
			if seqr[i:i+1]==dict_chr[tag][end:end+1]:
				left_check+=1
			elif cost_l<1:
				cost_l+=1
				left_check+=1
			else:
				break
			left_str1+=seqr[i:i+1]
		for j in range(15):
			if seqr[len_pointer-1-j:len_pointer-1-j+1]==dict_chr[tag][start-2-j:start-2-j+1]:
				right_check+=1
			elif cost_r<1:
				cost_r+=1
				right_check+=1
			else:
				break
			right_str1=seqr[len_pointer-1-j:len_pointer-1-j+1]+right_str1
		if left_check>=4 or right_check>=4:
			if left_check>=right_check:
				return left_str1,'L'
			elif right_check>left_check:
				return right_str1,'R'
		else:
			return '0','0'
dict_chr=dict()
results=pandas.read_table(sys.argv[2],sep="\t")
for rec in SeqIO.parse(sys.argv[1],"fasta"):
	dict_chr.update({rec.id:str(rec.seq)})
results['pointer'][i]=''
results['pointer_type'][i]=''
for i in results.index:
	this_seq=results['sequence'][i]
	start=results['start'][i]
	end=results['end'][i]
	chromosome=results['chromosome'][i]
	(pointer,type_t)=pointer_check("deletion",this_seq,chromosome,start,end)
	results['pointer'][i]=pointer
	results['pointer_type'][i]=type_t
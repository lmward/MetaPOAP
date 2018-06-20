#Import the math module to add factorials, which are essential for the combinatorial calculations in the False Negative test.
import math
#import the sys module to allow software exit
import sys
#Defining a function for a binomial coefficient, i.e. n choose k, used in the False Negative test.
def ncr(n,r):
	f=math.factorial
	return f(n)/(f(r)*f(n-r))
#Introduction to MetaPOAP, then prompts for user-entered genome and pathway statistics.
print "Welcome to MetaPOAP. You will be prompted for statistics about your MAG and pathway of interest before your MetaPOAP statistics can be calculated.\n"
#Prompt for number of coding sequences. This is estimated post-annotation by services such as RAST, but can also be estimated from bin size assuming average gene length ~1kb.
n=input("How many coding sequences were recovered in your MAG?\n")
if n<0:
	print "I'm sorry, MetaPOAP can't work with negative numbers."
	sys.exit()
if n==0:
	print "I'm sorry, MetaPOAP can't work with empty genomes."
	sys.exit()
#Prompt for genome completeness. We recommend using CheckM to estimate this value.
c=input("How complete is your MAG? Please enter a value between 0.0 and 1.0.\n")
if c<0:
	print "I'm sorry, MetaPOAP can't work with negative numbers."
	sys.exit()
elif c==0:
	print "I'm sorry, MetaPOAP needs completeness greater than 0 and less than 1."
	sys.exit()
elif c>1:
	print "I'm sorry, MetaPOAP needs completeness greater than 0 and less than 1."
	sys.exit()
#Prompt for genome contamination. We use the CheckM estimate of contamination, which reflects the fraction of single-copy genes for which multiple copies were recovered (i.e. this value can be greater than 1)
y=input("What is the estimated contamination in your MAG? This value is likely between 0 and 0.1 for high-quality MAGs, but can be >1.\n")
if y<0:
	print "I'm sorry, MetaPOAP can't work with negative numbers."
	sys.exit()
#Prompt for number of marker genes in the pathway of interest. This value should be conservative, avoiding genes that are difficult to annotate, i.e. have divergent homologs, genes that are multifunctional and could be present for a different purpose, and any other genes that are not good identifiers of the pathway of interest.
g=input("How many marker genes are in your complete pathway of interest? This value should be greater than zero.\n")
#exit if a negative number is entered
if g<0:
	print "I'm sorry, MetaPOAP can't work with negative numbers."
	sys.exit()
#exit if a zero value is entered
if g==0:
	print "I'm sorry, you need to enter a nonzero number for genes in your pathway."
	sys.exit()
#Prompt for marker genes recovered. These should be confidently identified, well-annotated genes that have been recovered in the genome bin.
r=input("How many marker genes from your pathway were recovered in your MAG? This value can be zero.\n")
#exit if a negative number is entered
if r<0:
	print "I'm sorry, MetaPOAP can't work with negative numbers."
	sys.exit()
if r>g:
	print "I'm sorry, you need to enter a number of recovered genes smaller than the total number of marker genes in your pathway."
	sys.exit()
m=input("Your marker genes were recovered across how many contigs? This value can be zero if you recovered no marker genes.\n")
#exit if a negative number is entered
if m==0 and r<0:
	print "I'm sorry, you can't have recovered marker genes that weren't on contigs."
	sys.exit()
if m>r:
	print "I'm sorry, you can't have recovered more contigs with marker genes than marker genes themselves."
	sys.exit()
#Calculation of a few values based on user input. T is the predicted total size of the genome as determined by recovered size and completeness.
T=math.ceil(n/c)
#x is the number of marker genes in the pathway of interest not recovered in the genome bin.
x=g-r
#this is the average number of marker genes on each contig encoding part of the pathway.
z=r/m
#this is the expected number of contigs encoding marker genes that were not recovered
k=x/z
#Q is the estimated fraction of the bin that is made up of contaminated sequences, correcting for the way CheckM reports contamination.
Q=y/(y+1)
#if no genes in the pathway were recovered, a False Positive statistic isn't meaningful and won't be reported.
if r==0:
	print "You didn't recover any genes in the pathway, so a False Positive estimate isn't applicable in this case."
	falsePos="N/A"
else:
	#Calculation of the false positive statistic. Assuming that each contig in the bin has an equal chance of belonging to the contamination fraction, this just multiplies the probability for each of the contigs encoding marker genes recovered to determine how likely it is that they all are contamination.
	falsePos=pow(Q,m)
	#Report of False Positive statistic
	print "False Positive Probability="
	print falsePos

#if genome is complete, or if a full pathway is recovered, false negative estimate is meaningless
if c==1:
	print "Your MAG is complete, so a False Negative estimate is not applicable in this case."
	falseNeg="N/A"
elif g==r:
	print "You recovered the full pathway in your MAG, so a False Negative estimate is not applicable in this case."
	falseNeg="N/A"

#Calculation of the False Negative statistic. Assuming that fragments of the genome are randomly sampled to make up the bin, and genes are approximately all the same length, what are the odds of missing all x of marker genes in the pathway that were not recovered.
#This calculation is performed in multiple steps for clarity and to maintain appropriate type for each variable.
else:
	Pa=float(ncr(n,0))
	a=int(T-n)
	Pb=float(ncr(a,k))
	Pc=float(ncr(T,k))
	falseNeg=float(Pa*Pb/Pc)
	print "False Negative Probability="
	print falseNeg
#Qualitative summary of False Pos and False Neg statistics, presented as confidence whether the pathway of interest is present in the genome.
if falsePos > 0.05:
	print "There's a good chance that the marker genes you recovered are contaminants."
if falseNeg > 0.05:
	print "There's a good chance that the marker genes you did not recover are actually in your source genome."
if falseNeg==0 and falsePos < 0.05:
	print "Your source genome likely contains your pathway of interest."
elif falseNeg > 0.05 and falsePos < 0.05:
	print "Your source genome likely contains your pathway of interest."
elif falsePos > 0.05 and falseNeg < 0.05:
	print "Your source genome probably doesn't contain your pathway of interest."
elif falsePos==0 and falseNeg < 0.05:
	print "Your source genome probably doesn't conntain your pathway of interest."
elif falsePos < 0.05 and falseNeg < 0.5 and r<g:
	print "It looks like your recovered marker genes belong to your source genome, but the missing genes may really be  missing. Your organism might contain a partial pathway."
else:
	print "It's hard to say whether your source genome contains your pathway of interest."

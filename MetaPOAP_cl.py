#Import the math module to add factorials, which are essential for the combinatorial calculations in the False Negative test.
import math
#import the sys module to allow software exit
import sys
#Defining a function for a binomial coefficient, i.e. n choose k, used in the False Negative test.
import argparse
#load command line arguments for MAG and pathway of interest
parser = argparse.ArgumentParser(description='Calculate False Positive and False Negative probabilites for presence of metabolic pathways given genome completeness and assembly. Accepts positional command line arguments (as described below), e.g.: ./MetaPOAP_cl.py 3078 0.87 0.1285 6 6 2')
parser.add_argument('n', metavar='coding sequences', type=int, nargs='+',help='Number of coding sequences recovered in the genome')
parser.add_argument('c', metavar='completeness', type=float, nargs='+',help='Genome completeness between 0 and 1')
parser.add_argument('y', metavar='contamination', type=float, nargs='+',help='Genome contamination, 0 or greater')
parser.add_argument('g', metavar='genes in pathway', type=int, nargs='+',help='Number of marker genes in the pathway of interest')
parser.add_argument('r', metavar='recovered genes', type=int, nargs='+',help='Number of marker genes in the pathway recovered in the genome')
parser.add_argument('m', metavar='marker gene congtis', type=int, nargs='+',help='Number of contigs recovered encoding marker genes')
args=parser.parse_args()
n=args.n[0]
c=args.c[0]
y=args.y[0]
g=args.g[0]
r=args.r[0]
m=args.m[0]
#define combination function
def ncr(n,r):
	f=math.factorial
	return f(n)/(f(r)*f(n-r))
#exit if disallowed values are entered
if n<0:
	sys.exit()
if n==0:
	sys.exit()
if c<0:
	sys.exit()
elif c==0:
	sys.exit()
elif c>1:
	sys.exit()
if y<0:
	sys.exit()
if g<0:
	sys.exit()
if g==0:
	sys.exit()
if r<0:
	sys.exit()
if r>g:
	sys.exit()
#calculate essential intermediate statistics, i.e. prediicted total genome size, number of unrecovered marker genes, and contaminant fraction of genome
T=math.ceil(n/c)
x=g-r
Q=y/(y+1)
#if no marker genes were recovered, a False Positive estimate won't be made.
if r==0:
	falsePos="N/A"
else:
	#Calculation of the false positive statistic. Assuming that each contig in the bin has an equal chance of belonging to the contamination fraction, this just multiplies the probability for each of the contigs encoding marker genes recovered to determine how likely it is that they all are contamination.
	falsePos=pow(Q,m)

#if genome is complete, or if a full pathway is recovered, false negative estimate is meaningless
if c==1:
	falseNeg="N/A"
elif g==r:
	falseNeg="N/A"

#Calculation of the False Negative statistic. Assuming that fragments of the genome are randomly sampled to make up the bin, and genes are approximately all the same length, what are the odds of missing all x of marker genes in the pathway that were not recovered.
#This calculation is performed in multiple steps for clarity and to maintain appropriate type for each variable.
else:
	Pa=float(ncr(n,0))
	a=int(T-n)
	Pb=float(ncr(a,x))
	Pc=float(ncr(T,x))
	falseNeg=float(Pa*Pb/Pc)
#print False Positie and False Negative statistics
print falsePos
print falseNeg

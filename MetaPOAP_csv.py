import csv
import math
import sys
import argparse
#load inut and output files based on command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('i', metavar='inputfile', nargs='+',help='input csv file')
parser.add_argument('o', metavar='outputfile', nargs='+',help='output csv file')
args=parser.parse_args()
i=args.i[0]
o=args.o[0]
#define combination function
def ncr(n,r):
	f=math.factorial
	return f(n)/(f(r)*f(n-r))
#load in and outfiles
infile = open(i, 'rb')
reader = csv.reader(infile, delimiter=',')
outfile = open(o, 'w')
writer=csv.writer(outfile)
output=""
#iterate through each line of the infile
with open(i) as csvfile:
	readCSV=csv.reader(csvfile, delimiter=',')
#treat first line of infile as header, skipping input and printing headers to outfile
	firstline = True
	for row in readCSV:
		if firstline:    
			writer.writerow(["MAG ID","Pathway","Coding sequences","Completeness","Contamination","Marker genes in pathway","Marker genes recovered","Total contigs","Contigs with marker genes","False Positive estimate","False Negative estimate"])
			firstline = False
			continue
#read each column of the infile and keep values with appropriate type
		ID=row[0]
		pathway=row[1]
		n=row[2]
		n=int(n)
		c=row[3]
		c=float(c)
		y=row[4]
		y=float(y)
		g=row[5]
		g=int(g)
		r=row[6]
		r=int(r)
		s=row[7]
		s=int(s)
		m=row[8]
		m=int(m)
		e="Error"
#check that values are allowed ranges, if not output error and move on to next line of infile
		if n<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if n==0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if c<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		elif c==0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		elif c>1:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if y<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if g<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if g==0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if r<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if r>g:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if s==0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if s<0:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if m>s:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if m>r:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
		if m>g:
			writer.writerow([ID,pathway,n,c,y,g,r,s,m,e,e])
			continue
#calculate intermediatte values
		T=math.ceil(n/c)
		x=g-r
		Q=y/(y+1)
#if estimates aren't meaningful, pass N/A, otherwise calculated False Negative and False positive values
		if r==0:
			falsePos="N/A"
		else:
			falsePos=pow(Q,m)
		if c==1:
			falseNeg="N/A"
		elif g==r:
			falseNeg="N/A"
		else:
			Pa=float(ncr(n,0))
			a=int(T-n)
			Pb=float(ncr(a,x))
			Pc=float(ncr(T,x))
			falseNeg=float(Pa*Pb/Pc)
#write output to outfile
		writer.writerow([ID,pathway,n,c,y,g,r,s,m,falsePos,falseNeg])
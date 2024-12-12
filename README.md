Small codes:

A) For VolcanoPlot.R:

Usage: 

	Rscript VolcanoPlot.R [options]

Options:

	-m CHARACTER, --myMatrix=CHARACTER
		matrix file: exported from perseus. 
  		The first, second and third column of this file should be: Significant, -Log(P)value and Difference respectively. 
    		Remove - sign in -Log(P)value before use

	-c CHARACTER, --myCurve=CHARACTER
		curve file: exported from perseus

	-p CHARACTER, --proteins=CHARACTER
		list of proteins you want to highlight. 
  		Note that name of this protein should match with name present in your matrix file. 
    		Seperate proteins of each category by semi-column e.g. Jun,FOS as AP1 category and NFYA is in another category, 
      		then use option as follows: --proteins JUN,FOS;NFYA. This program can handle only 5 complex. 
		Let me know if you need to annotate more that 5 complex.

	-b CHARACTER, --bait=CHARACTER
		Name of bait protein. 
  		Note that name of this protein should match with name present in your matrix file.

	-t CHARACTER, --title=CHARACTER
		Title of image

	-h, --help
		Show this help message and exit

Example:
	
 	Rscript VolcanoPlot.R \
  	-m MS-1matrix.txt \
   	-c MS-1curve.txt \
    	--proteins "WDR5,DPY30,ASH2L,KMT2D,PAXIP1,NCOA6,KMT2C,RBBP5;U2AF1,HSPA8" \
     	--bait UTY \
      	--title "UTY nuclear"

B) For genometrack.R

Usage: 
	Rscript genomttrack.R [option]

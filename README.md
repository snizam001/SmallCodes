Small codes:

A) For VolcanoPlot.R:

	This script will generate volcano plot if you have perseus matrix and curve data.
 
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

B) For genomeTrack.R

	This script will generate genome track files. 

Usage: 
	Rscript genomeTrack.R [option]

Options:

	-i CHARACTER, --inputFiles=CHARACTER
		Input bigwig file. Give full path and provide as comma separated values

	-o CHARACTER, --outputImagefile=CHARACTER
		Name of the output image files

	-c CHARACTER, --chr=CHARACTER
		chromosome example chr1 (always use "chr")

	-s CHARACTER, --start=CHARACTER
		starting location of genomic region (hg38)

	-e CHARACTER, --end=CHARACTER
		ending location of genomic region (hg38)

	-g CHARACTER, --group=CHARACTER
		provide group of bw file. Suppose you have three input file and you want to group first two samples in one and third as another, then give it as follows: --group 1,1;2

	-h, --help
		Show this help message and exit

Example:
	
 	infiles="MV411-JUN-NT1_expr.ScaledSpikeIn.bw,MV411-JUN-G4_expr.ScaledSpikeIn.bw,MV411-MEN-NT1_expr.ScaledSpikeIn.bw,MV411-MEN-G4_expr.ScaledSpikeIn.bw,MV411-MLL1-NT1_expr.ScaledSpikeIn.bw,MV411-MLL1-G4_expr.ScaledSpikeIn.bw"
 
	Rscript genomeTrack.R  \
 	-i $infiles \
	-o genomicTrack1.jpeg \
	--chr chr6 \
	--start 170554302 \
	--end 170554303 \
	--group "1,1;2,2;3,3;4,4"

Rscript genomeTrack.R 

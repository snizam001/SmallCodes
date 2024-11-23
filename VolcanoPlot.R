#!/usr/bin/env Rscript
# example Rscript MSplot.R -m MS-1matrix.txt -c MS-1curve.txt -p "WDR5,DPY30,ASH2L,KMT2D,PAXIP1,NCOA6,KMT2C,RBBP5;U2AF1,HSPA8" -b UTY -t "UTY nuclear"



if(!require(ggplot2)){
    print ("ggplot2 is not avaiable in your system")
        BiocManager::install("ggplot2")
        }


if(!require(ggrepel)){
    print ("ggrepel is not avaiable in your system")
        BiocManager::install("ggrepel")
        }

if(!require(optparse)){
        install.packages("optparse",repos = "http://cran.us.r-project.org")}
library(optparse)

# note: delete - sign of -log(pvalue) from matrix file
# note: you can use these colors according to your taste: https://nanx.me/ggsci/reference/pal_npg.html

option_list = list(
make_option(c("-m", "--myMatrix"),
    type="character",
    default=NA,
    help="matrix file: exported from perseus. This first, second and third column of this file should be: Significant, -Log(P)value and Difference respectively. Remove - sign in -Log(P)value before use",
    metavar="character"
           ),

make_option(c("-c", "--myCurve"),
    type="character",
    default=NA,
    help="curve file: exported from perseus",
    metavar="character"),

make_option(c("-p", "--proteins"),
    type="character",
    default=NA,
    help="list of proteins you want to highlight. Note that name of this protein should match with name present in your matrix file. Seperate proteins of each category by semi-column e.g. Jun,FOS as AP1 category and NFYA is in another category, then use option as follows: --proteins JUN,FOS;NFYA. This program can handle only 5 complex. Let me know if you need to annotate more that 5 complex.",
    metavar="character"),

make_option(c("-b", "--bait"),
    type="character",
    default=NA,
    help="Name of bait protein. Note that name of this protein should match with name present in your matrix file.",
    metavar="character"),

make_option(c("-t", "--title"),
    type="character",
    default="",
    help="Title of image",
    metavar="character")

);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Rscript MSplot.R -m MS-1matrix.txt -c MS-1curve.txt -p "WDR5,DPY30,ASH2L,KMT2D,PAXIP1,NCOA6,KMT2C,RBBP5;U2AF1,HSPA8" -b UTY -t "UTY nuclear"


inputFile=opt$myMatrix
myCurve=opt$myCurve
proteins=opt$proteins
myBait=opt$bait
title=opt$title 

#---
pp = unlist(strsplit(proteins,";"))


myColors = c("#3C5488","#E64B35","#8491B4","darkgreen","#DC0000")
# change this, this are target proteins you want to highlight in figure, name should match with the name of protein in matrix file

#-----------------------------

mybait = myBait

data=read.table(inputFile, header = T, sep = "\t")

colnames(data)[1]="Significant"
colnames(data)[2]="X.Log.P.value."
colnames(data)[3]="Difference"

myPvalueMax = round(max(data[,"X.Log.P.value."], rm.na = T ) + 2, 1 )
myDifferenceMax = round(max(data[,"Difference"], rm.na = T ) + 2, 1)
myDifferenceMin = round(min(data[,"Difference"], rm.na = T ) + 2, 1)
mySig = data[data$Significant=="+",]

myC = read.table(myCurve, header = T)

#--
p <- ggplot(data, aes(Difference, X.Log.P.value.)) +
  geom_point(color = "grey90") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Differences (experiment - control)", y = "-log10(adjusted p value)", title = title)
p = p + geom_line(data=myC,aes(x=x,y=y))
p = p + xlim(myDifferenceMin, myDifferenceMax) + ylim (0, myPvalueMax)
#--
p <- p + geom_point(data = mySig, aes(x = Difference, y = X.Log.P.value.), shape=20, colour = "grey40", fill = "grey40", size = 2) 

#--------------
TotalProteins=c()

for(i in c(1:length(pp))){
	TotalProteins = c(TotalProteins,
										unlist(strsplit(pp[i],","))
										)
}

filtered = data[data$Gene.names %in% TotalProteins, ]

filtered = rbind(filtered,data[data$Gene.names %in% mybait, ])

filtered_d = filtered
#--------------

mycolors=rep("", length(TotalProteins)+1)

mycolors = ifelse(
	filtered$Gene.names %in% mybait , "black",
	ifelse (
		filtered$Gene.names %in% unlist(strsplit(pp[1],",")), myColors[1],
		ifelse(
			filtered$Gene.names %in% unlist(strsplit(pp[2],",")), myColors[2],
			ifelse(
				filtered$Gene.names %in% unlist(strsplit(pp[3],",")), myColors[3],
				ifelse(
					filtered$Gene.names %in% unlist(strsplit(pp[4],",")), myColors[4],
					ifelse(
						filtered$Gene.names %in% unlist(strsplit(pp[5],",")), myColors[5],
						"cyan"
					)
				)
			)
		)
	)
)

#--------------
p <- p + geom_point(data = filtered, aes(x = Difference, y = X.Log.P.value.), shape=20, colour = mycolors, fill = "darkblue", size = 2)  
p = p + geom_text_repel(data = filtered_d, aes(x = Difference, y = X.Log.P.value., label = Gene.names), max.overlaps = Inf, colour =mycolors)


jpeg(paste(title, ".jpeg", sep = ""), unit= "in", res=300, height=4.5, width = 4.5)
p

dev.off()


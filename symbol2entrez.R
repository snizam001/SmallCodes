#!/usr/bin/env Rscript
#ChIPpeakAnno
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require(ChIPpeakAnno)){
        BiocManager::install("ChIPpeakAnno")
    library(ChIPpeakAnno)
}


if(!require(optparse)){
        install.packages("optparse")
    library(optparse)
}

if(!require(org.Hs.eg.db)){
        install.packages("org.Hs.eg.db")
    library(org.Hs.eg.db)
}

if(!require(limma)){
        BiocManager::install("limma")
    library(limma)
}





option_list = list(
    
    make_option(c("-i", "--input"), 
                type="character", 
                default=NA,
                help="Input list of the gene symbols", 
                metavar="character")
)


parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)

file=opt$input
d= read.table(file)

geneList = limma::alias2Symbol (as.character(d[,1]),species="Hs")

entrezIDs = convert2EntrezID(geneList,orgAnn="org.Hs.eg.db",ID_type = "gene_symbol" )
print(entrezIDs)

write.table(entrezIDs,'entrez.txt', quote=F,sep="\t", row.names =F , col.names =F )

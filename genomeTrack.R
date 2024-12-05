#!/usr/bin/env Rscript
if(!require(optparse)){
        install.packages("optparse",repos = "http://cran.us.r-project.org")}
library(optparse)

#outputImagefile="check.jpeg"
#inputFiles="/media/sheikh/Swift0/CutAndRun/Publication_hg38/bamcoverage/enh-30.bw,/media/sheikh/Swift0/CutAndRun/Publication_hg38/bamcoverage/NFYA-10.bw,/media/sheikh/Swift0/CutAndRun/Publication_hg38/bamcoverage/NFYA-90.bw"
#chromosomewa="chr1"
#chromosome_forgene=gsub(pattern="chr",
#                            replacement="",
#                            x=chromosomewa)
#starts=155650443
#ends=155650642
#Groups="1,1;2;"




option_list = list(
make_option(c("-i", "--inputFiles"),
    type="character",
    default=NA,
    help="Input bigwig file. Give full path and provide as comma separated values",
    metavar="character"),

make_option(c("-o", "--outputImagefile"),
    type="character",
    default=NA,
    help="Name of the output image files",
    metavar="character"),

make_option(c("-c", "--chr"),
    type="character",
    default="NA",
    help="chromosome example chr1 (always use \"chr\")",
    metavar="character"),

make_option(c("-s", "--start"),
    type="character",
    default="NA",
    help="starting location of genomic region (hg38)",
    metavar="character"
           ),

make_option(c("-e", "--end"),
    type="character",
    default="NA",
    help="ending location of genomic region (hg38)",
    metavar="character"
           ),

make_option(c("-g", "--group"),
    type="character",
    default="NA",
    help="provide group of bw file. Suppose you have three input file and you want to group first two samples in one and third as another, then give it as follows: --group 1,1;2",
    metavar="character"
           )
);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
outputImagefile=opt$outputImagefile
inputFiles=opt$inputFiles
chromosomewa=opt$chr
chromosome_forgene=gsub(pattern="chr",
                            replacement="",
                            x=chromosomewa)
starts=as.numeric(opt$start)
ends=as.numeric(opt$end)
Groups=opt$group

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

if(!require(rtracklayer)){
    print ("rtracklayer is not avaiable in your system")
        BiocManager::install("rtracklayer")
        }
library(rtracklayer)

if(!require(AnnotationHub)){
    print ("AnnotationHub is not avaiable in your system")
        BiocManager::install("AnnotationHub")
        }
library(AnnotationHub)

if(!require(Rsamtools)){
    print ("Rsamtools is not avaiable in your system")
        BiocManager::install("Rsamtools")
        }
library(Rsamtools)

if(!require(Gviz)){
    print ("Gviz is not avaiable in your system")
        BiocManager::install("Gviz")
        }
library(Gviz)

if(!require(areaplot)){
        install.packages("areaplot",repos = "http://cran.us.r-project.org")}
library(areaplot)

if(!require(EnsDb.Hsapiens.v86)){
    print ("EnsDb.Hsapiens.v86 is not avaiable in your system")
        BiocManager::install("EnsDb.Hsapiens.v86")
        }
library(EnsDb.Hsapiens.v86)


# reading files

files=unlist(
        strsplit(
            inputFiles,
            ","))


myGroup = if(!is.na(Groups)) {
    unlist(
        strsplit(
            Groups,
            ";"))
}


color=c('#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4",'#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4",'#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4",'#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4",'#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4",'#c60000','grey50','orange','grey72','hotpink4','darkblue','red4','black',"khaki4")
#color=wes_palette(n=length(files), name="AsteroidCity1")

genes(EnsDb.Hsapiens.v86) -> gene

data(geneModels)


# checking given range

if ((ends-starts) <= 100000) {
    midPoint = starts + round(
                        ends-starts, 0
                            )
    ends = midPoint + 50000
    starts = midPoint - 50000

    if(starts < 0) {starts = 0}

    cat("----- Given genomic location was < 100,000 base-pairs. Range is expanded to 50,000 on both side from mid point of given range.\n")
}

subset=GRanges(seqnames=chromosomewa,ranges=IRanges(start = starts, end = ends))


#-- 
# checking if file path exists and is not empty

error = 0;
for(i in c(1:length(files)))
    {
        if (file.exists(files[i]) && file.info(files[i])$size > 0) {} 
        else {
          cat(paste(files[i], ": does not exist or is empty.\n", sep = ""))
          error = error + 1
        }
    }

if(error > 0) {quit(save = "no")} 

# reading files-bigwig

bwData = list()
bwDataSubset = list()

for(i in c(1:length(files)))
    {
        print(paste("-- reading ", files[i],sep=""))
        bwData[[i]] = import(files[i], format = 'bw')
        bwDataSubset[[i]]=subsetByOverlaps(bwData[[i]],subset)
    }

# generating image

jpeg(outputImagefile,
    unit="cm",
    res=300,
    height=0.3*((length(files)) + 3 ),
    width=8.5)

par(mar=c(0,4.1,0,2.1))

opar<-par(mfrow=c(
    length(files) + 2, 1
    )
)

# first one is empty
plot(c(0,1),c(0,1), axes=F, xlab = "", ylab ="", main = "", type= "n")
#-----
j = 0
for (i in c(1:length(myGroup))) {

    data4plot = list()
    aa = unlist(
        strsplit(
            myGroup[[i]],
            ","))

    maxValue = 0
    for(k in c(1:length(aa))) {
        j = j + 1
        data4plot[[k]] = bwDataSubset[[j]]

        if (maxValue < max(score(data4plot[[k]]) )) {maxValue = max(score(data4plot[[k]]))}
    }

    maxValue = round( maxValue + (0.1 * maxValue), 0 )
    col=color[i]

    for(k in c(1:length(aa))) {
        areaplot(
            (start(data4plot[[k]])+end(data4plot[[k]]))/2,
            score(data4plot[[k]]),
            col=col, border = col,
            axes=F,xlab="",
            ylab="",
            xlim=c(starts,ends),
            ylim=c(0,maxValue))

        axis(2,at=maxValue,labels=as.character(maxValue),las=2)
        abline(h=0)
    }

}

#----
subset_gene=GRanges(seqnames=chromosome_forgene,ranges=IRanges(start = starts, end = ends))
file = subsetByOverlaps(gene,subset_gene)
diff=10/(2*nrow(data.frame(file)))
plot((start(file)+end(file))/2,rep(10,nrow(data.frame(file))),type='n',axes=F,xlab="",ylab="",xlim=c(starts,ends),ylim=c(0,10))
y0=diff

for(i in c(1:nrow(data.frame(file)))){
    code=ifelse(strand(file)[i] == "+",2,1)
    arrows(start(file)[i],y0,end(file)[i],y0,code=code,col='darkgreen',length=0.025,lwd=1)
    y0=y0+diff

}

dev.off()


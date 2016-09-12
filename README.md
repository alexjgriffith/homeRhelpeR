# homeRhelpeR
## Homer Helper for R
### Alexander Griffith
### Sept 11 2016


---
### Installation

``` sh
## Get the source from github
cd /tmp  # optional
git clone https://github.com/alexjgriffith/homeRhelpeR.git
cd homeRhelpeR

## Install 
R CMD INSTALL .

## Once it installs correctly remove the repository
cd ../
rm -rf homeRhelpeR
```

---
### Using Homer Wrapper
Please note that you need homer2 installed on your system visit [homer install](http://homer.salk.edu/homer/introduction/install.html) to aquire the source code. For this example my executable is found in "~/binaries/homer/bin/homer2". 

```R
library(homeRhelpeR)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

## location of executable on your system
homer="~/binaries/homer/bin/homer2"
## Bed files for the foreground and background
fg<-read.table("foreground.bed",header=T,comment=NULL)
bg<-read.table("background.bed",header=T,comment=NULL)
bed<-rbind(fg,bg)
## find the neucleotide info for the bed files
fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed$chr,
	start=bed$start,end=bed$end)
## Which indicies of fasta belong to the foreground and background	
len<-dim(bed)[1]
foreground<-background<-rep(FALSE,len)
foreground[1:dim(fg)[1]]<-TRUE
background[(dim(fg)[1]+1):len]<-TRUE

## Apply the homer2 function looking for 25 motifs of length 6
motifs<-homerWrapper(fasta,
	foreground,
	background,
	homer,
	opts="-S 25 -len 6")

## save the motifs
write.homer(motifs,"test_6.motifs")
	
```

---
### Grouping Motifs
```R
library(homeRhelpeR)
library(Biostrings)

## find the neucleotide info for the bed files used to generate the motifs
fgbed<-read.table("foreground.bed",header=T,comment=NULL)
bgbed<-read.table("background.bed",header=T,comment=NULL)
fg<-getSeq(BSgenome.Hsapiens.UCSC.hg19,fgbed$chr,
	start=fgbed$start,end=fgbed$end)
bg<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bgbed$chr,
	start=bgbed$start,end=bgbed$end)

filenames<-c("test_6.motifs","test_7.motifs","test_8.motifs")
motifsList<-sapply(filenames,read.homer)

motifs<-unique(sapply(motifsList,getNames))

align<-findPairwise(motifsList)

names(align)<-motifs

ds<-combn2dist(align)

mots<-combineClusterdMotifs(ds,sapply(motifs,genPwm))

homerList<-motifs2Homer(fg,bg,mots)

write.homer(homerList,"combn.motifs")

```


---
### Calling Stamp on One set of motifs

```R
library(homeRhelpeR)
library(httr)
library(xlsx)

homerList<-read.homer("combn.motifs")

res<-callStamp(homerList,"JASPAR_Fams")

write.homer.db(res)

saveasSingleXLSX(res)

```

---
### Calling Stamp for may sets of motifs

```R
library(homeRhelpeR)
library(httr)
library(xlsx)

filenames<-c("test_6.motifs","test_7.motifs","test_8.motifs")
homerList<-lapply(filenames,read.homer)

res<-callMultiStamp(homerList,filenames,"JASPAR_Fams")

write.homer.db(res)

saveasMultiXLSX(res)

```


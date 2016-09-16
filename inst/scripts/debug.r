library(Biostrings)
library('BSgenome.Hsapiens.UCSC.hg19')
library(stringr)

motifs<-c("CANNTG","CCACATG","GATTAA","GATAA","HGATAA")


load("~/Dropbox/Data/ccca_20.RData")

fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,ccca$bed$chr,start=(ccca$bed$start+ccca$bed$end)/2-150,width=300)

fg<-fasta[ccca$reg[,"Leukemia"]]
bg<-fasta[ccca$reg[,"NONE"]]


## save data for package

with(ccca,{
    write.table(file="inst/extdata/leukemia.bed",x=bed[reg[,"Leukemia"],],quote=FALSE,row.names=FALSE)
    write.table(file="inst/extdata/erythroid.bed",x=bed[reg[,"Erythroid"],],quote=FALSE,row.names=FALSE)
    write.table(file="inst/extdata/hsc.bed",x=bed[reg[,"HSC"],],quote=FALSE,row.names=FALSE)
    write.table(file="inst/extdata/ecfc.bed",x=bed[reg[,"ECFC"],],quote=FALSE,row.names=FALSE)
    write.table(file="inst/extdata/background.bed",x=bed[reg[,"NONE"],],quote=FALSE,row.names=FALSE)
})


## Examples

homerLoc="~/binaries/homer/bin/homer2"

foregroundFile<-system.file("/extdata/leukemia.bed",package="homeRhelpeR")
backgroundFile<-system.file("/extdata/background.bed",package="homeRhelpeR")


fg<-read.table(foregroundFile,header=T,comment="")

bg<-read.table(backgroundFile,header=T,comment="")
bed<-rbind(fg,bg)

fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed$chr,
	start=bed$start,end=bed$end)

len<-dim(bed)[1]
foreground<-background<-rep(FALSE,len)
foreground[1:dim(fg)[1]]<-TRUE
background[(dim(fg)[1]+1):len]<-TRUE

## Apply the homer2 function looking for 25 motifs of length 6
motifs<-lapply(c(6,7,8), function(i) homerWrapper(fasta,
                                      foreground,
                                      background,
                                      homerLoc,
                                      opts=paste("-S 2 -len ",i)))


smots<-sapply(do.call(mergeMotifLists,motifs),getCons.homer)


combineClusterdMotifs(clusterMotifs(smots),smots,3)

motifs2Homer(fasta[foreground],fasta[background],combineClusterdMotifs(clusterMotifs(smots),smots,3))


str(motifs2Homer(fg,bg,combineClusterdMotifs(clusterMotifs(motifs),motifs,3)))


homerFiles<-strsplit("ECFC_sqrt_PCA Leukemia_PCA Erythroid_PCA HSC_PCA ECFC_BED Leukemia_BED Erythroid_BED HSC_BED"," ")[[1]]

fileLocs<-paste0("~/Dropbox/Data/homer-paper/homer_",sapply(homerFiles,rep,3),"_",c(6,7,8),".txt")
homerData<-lapply(fileLocs,read.homer)

names(homerData)<-fileLocs


library(stringr)
library(httr)

source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")


mergeHomer(homerData[[1]])

par(mfrow=c(2,2))


names(MOTS[[1]])[MOTS[[1]]==1]


PCA<-lapply(contexts,function(cont){
    reg<-grepl(cont,fileLocs) & grepl("PCA",fileLocs)
    mots<-unique(Filter(function(x) nchar(x)>4, trimMotif(sapply(mergeMotifLists(homerData[reg]),getCons.homer))))
    print(mots)
    combineClusterdMotifs(clusterMotifs(mots),mots,0.7)
})

names(PCA)<-contexts

motifs<-mapply(function(cont,x) motifs2Homer(ccca$fasta[ccca$reg[,cont]],ccca$fasta[ccca$reg[,"NONE"]],x),contexts,PCA,SIMPLIFY=FALSE)

stamp<-local({    
    stampdb<-mapply(function(mots,name){
        rmna<-mots[sapply(mots,function(x) !all(is.na(x)))]
        homerStamp.motifList(rmna,c("JASPAR_Fams","TRANSFAC_Fams"),name)}
                   ,motifs,contexts,SIMPLIFY=FALSE)
    names(stampdb)<-contexts
    save(stampdb,file="~/Dropbox/Data/PCAStamp.RData")
    stampdb
})

ts<-homerStamp.motifList(motifs[[1]],c("JASPAR_Fams"),"test")


ts

ts

library(xlsx)

saveMultiXLS(stamp,motifs,"~/Desktop/test.xlsx")

#myms<-msa(mots,type="dna",cluster="upgma")


library(httr)
#library("homeRhelpeR")

## Fix the issue with the homerStamp method (stamp returning motifs
## out of order
a<-stamp[[1]][,1:6]

b<-stamp[[1]][,7:10]

pb<-switchAttr(b[!duplicated(b[,1]),],"motif","tmotif")

merge(a,pb)[,c("motif","name2")]


a<-ts[,1:6]

b<-ts[,7:10]

pb<-switchAttr(b[!duplicated(b[,1]),],"motif","tmotif")

merge(a,b)


[,c("motif","name2")]

switchAttr<-function(obj,name,oldid,attr="names"){
    oldattr<-attributes(obj)
    oldnames<-oldattr[[attr]]
    if(is.null(oldnames)){
        warning(attr," is not a valid attribute of obj")
        return(obj)
    }
    if(! any(oldnames==oldid)){
        warning(oldid," not in obj ",attr)
        return(obj)
    }
    oldnames[oldnames==oldid]<-name
    oldattr[[attr]]<-oldnames
    `attributes<-`(obj,oldattr)        
}




n<-lapply(mots,function(m)    
unlist(lapply(getKmers(m,4),IUPACtoBase,TRUE)))

kmers<-function(n){
    IUPACtoBase(do.call(paste0,as.list(rep("N",n))),TRUE)
}

library(stringr)

reg<-grepl("Erythroid",fileLocs) & grepl("PCA",fileLocs)
mots<-unique(Filter(function(x) nchar(x)>4, trimMotif(sapply(mergeMotifLists(homerData[reg]),getCons.homer))))

n<-str_match_all(sapply(mots,IUPACtoBase,TRUE),paste0(kmers(4),collapse="|"))

Xs<-outer(n,n,Vectorize(function(a,b) 1-length(intersect(a,b))/length(union(a,b))))

colnames(Xs)<-mots
rownames(Xs)<-mots

dist<-as.dist(Xs)

hc<-hclust(dist,"complete")
a<-cutree(hc,h=0.8)





dist<-clusterMotifs(mots)

combineClusterdMotifs(dist,mots,0.8)

hc<-hclust(1/(clust+1),"average")

plot(hc,hang=-1)

lines(rep(0.2,56))

cutree(hc,h=0.2)

library("xlsx")


library(combinat)

file<-("~/Desktop/test.xlsx")



stampWrapper<-function(motifFile,comparator,name=motifFile,...){
    postFormB<-function(motif,comparator){
        motif_file<-path.expand(motif)
        list(input="",
             MAX_FILE_SIZE="10000000",
             motiffile=upload_file(motif_file),
             mtrim="on",
             ICbox="0.4",
             num_hits="5",
             match_db=comparator,
             MAX_FILE_SIZE="1000000",
             metric="PCC",
             align="SWU",
             mulalign="IR",
             tree="UPGMA",
             submit="Submit",
             metric_des="Pearson%27s+Correlation+Coefficient%3A%0D%0A-------------------------------%0D%0A%28Recommended%29%0D%0AMin%3A+-1%2C+Max%3A+%2B1+%28per+column%29.%0D%0A ",
             align_des="Ungapped+Smith-Waterman%3A%0D%0A------------------------------%0D%0A%28Recommended%29%0D%0ALocal+alignment+for+comparisons+between+short+motifs.",
             mulalign_des="Iterative+Refinement%3A%0D%0A------------------------------%0D%0AUpdates+multiple+alignment+by+leave-one-out+iteration.",
             tree_des="UPGMA%3A%0D%0A------------------------------%0D%0A%28Recommended%29%0D%0ADistance-based+tree+construction.+")
    }
    ua<- "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:33.0) Gecko/20100101 Firefox/33.0"
    body<-postFormB(motifFile,comparator)
    a<-POST("http://www.benoslab.pitt.edu/stamp/run.php" ,body=body,user_agent(ua))
    datalink<-strsplit(strsplit(grep("parent.location",strsplit(content(a,"text"),"\n")[[1]],value=TRUE)," ")[[1]][5],"'")[[1]][2]
    b<-GET(datalink)
    motifData<-data.frame()
    for(i in (seq(100)-1)){
        c<-grep(paste("\"motif",i,"\"",sep=""),strsplit(content(b,"text"),"\n")[[1]],value=TRUE)
        if(length(c)==0)
            break
        eScores<-as.numeric(unlist(lapply(grep("width=\"75\"",strsplit(c,"</TD>")[[1]],value=TRUE),function(x) strsplit(x,">")[[1]][2]))[2:6])
        ann<-unlist(lapply(grep("width=\"175\"",strsplit(c,"</TD>")[[1]],value=TRUE),function(x) strsplit(x,">")[[1]][2]))[2:6]
    motif<-unlist(strsplit(grep("<strong>",strsplit(c,"</strong></")[[1]],value=TRUE),">"))
        t<-data.frame(motif=motif[length(motif)],rank=seq(5),ann=ann,escore=eScores,dataset=name,comparator=comparator)
        motifData<-rbind(motifData,t)
    }
    motifData
}

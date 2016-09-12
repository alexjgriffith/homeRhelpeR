logPvalueFH<-function(filename){
    homer<-loadPWM(as.character(filename)) # rename loadHomer and
                                           # remove jaspar capabilites
    pwms<-strsplit(do.call(rbind,sapply(homer[,2],strsplit,","))[,3],":")
    as.numeric(sapply(sapply(strsplit(sapply(pwms,"[[",2),"T"),"[",1),function(x) substr(x,4,nchar(x))))
}

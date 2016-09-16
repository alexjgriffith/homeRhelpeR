#' @export
buildFM<-function(stamp,n=3,pvalue=1e-3){
    comparators<-unique(as.vector(stamp[,"comparator"]))    
    reg<-as.vector(stamp[,"comparator"])==comparators[1] & stamp[,"rank"]==1
    motifInfo<-stamp[reg,c("motif","order","pvalues")]
    complement<-sapply(as.vector(motifInfo[,"motif"]),function(x)consensusIUPAC(complement(IUPACtoBase(x))))
annotations<-do.call("cbind",lapply(comparators,function(comp)do.call(cbind,lapply(seq(n),function(i){
    reg<-as.vector(stamp[,"comparator"])==comp & stamp[,"rank"]==i
    stamp[reg,c("escore","ann")]
}))))
    ret<-cbind(motifInfo,complement=complement,annotations)
    #ret<-ret[ret$pvalues<pvalue,]
    ret[order(ret$pvalues),]
}


#' @export
saveXLS<-function(stamp,file,name="stamp"){
    require("xlsx")
    wb<-createWorkbook()
    ts<-createSheet(wb,name)        
    srow=1
    bf<-buildFM(stamp)
    addDataFrame(bf,sheet=ts,row.names = FALSE,startRow = srow)        
    saveWorkbook(wb,path.expand(file))
}


#' @export
saveMultiXLS<-function(stamp,lists,file){
    require("xlsx")
    wb<-createWorkbook()

    lapply(names(stamp),function(i){
        ts<-createSheet(wb,paste(i,"ann",sep="_"))
        bf<-buildFM(stamp[[i]])
        addDataFrame(bf,sheet=ts,row.names = FALSE)
        th<-createSheet(wb,paste(i,"header",sep="_"))
        header<-as.data.frame(do.call(rbind,lapply(lists[[i]],function(x) do.call(c,x$header))))
        header$pvalue<-as.numeric(as.vector(header$pvalue))
        header$score<-as.numeric(as.vector(header$score))
        header<-header[exp(header$pvalue)<1e-3,]
        addDataFrame(header[order(header$pvalue),],sheet=th,row.names = FALSE)

        tb<-createSheet(wb,paste(i,"homer",sep="_"))

        body<-as.data.frame(do.call(rbind,Filter(function(x) !is.null(x),lapply(motifs[[1]], function(mots) {if(exp(mots$pvalues)<1e-3) {rbind(mots$header[1:4],mots$pwm)} else {NULL}}))[order(header$pvalue)]))
        addDataFrame(body,sheet=tb,row.names = FALSE,col.names=FALSE)

        
    })
    saveWorkbook(wb,path.expand(file))
}

checkHomerContents<-function(conts){    
    required<-c("consensus","name","score")
    optional<-c("pvalue","gapped","pstring","stats")
    if(is.null(conts))
       return (required)
    ## ensure that the required contents are in conts
    if(! any(is.na(match(required,conts)))){
        rem<-Filter(function(x){
            if(x %in% optional){
                TRUE}
            else{
                warning(x," is not in ",paste(optional,collapse=" "))
                FALSE
            }
        }
       ,setdiff(conts,required))
        return(c(required,rem))
    }
    else{
        stop("consensus,name, and threshold must be included in contents")
    }
}



#' Read Homer File
#' 
#' http://homer.salk.edu/homer/motif/creatingCustomMotifs.html
#' @export
read.homer<-function(fileLocation)
{
    
    selectHomerHeader<-function(table){
        required<-c("consensus","name","score")
        optional<-c("pvalue","gapped","pstring","stats")
        contents<-c(required,optional);
        whichH<-substr(table[,1],1,1)==">"
        header<-table[whichH,]
        max<-length(contents)
        diff=length(contents)-length(header)
        if(diff>4)
            stop("To few header elements")
        else if(diff>0)
            warning(paste(contents[(max-diff+1):max],collapse=", ")," are missing.")
        else if (diff< 0)
            warning("More header elements than accounted for")
        
        ret<-apply(header,1,function(x){
            x<-as.list(x)
            names(x)<-contents[1:(max-diff)]
            x}
            )

        ret
        
    }
    selectHomerBody<-function(table){
        whichH<-which(substr(table[,1],1,1)==">")
        pwmRanges<-cbind(whichH+1,c(whichH[2:length(whichH)]-1,dim(table)[1]))
        pwms<-apply(pwmRanges,1,function(x) {y<-table[x[1]:x[2],c(1,2,3,4)];colnames(y)<-c("A","C","G","T");y})
        lapply(pwms,function(pwm) as.matrix(apply(pwm,2,as.numeric)))
    }

    table<-read.table(fileLocation,fill=TRUE)
    header<-selectHomerHeader(table)
    names<-sapply(header,function(x) x$name)
    names(header)<-names
    #contents<-names(header)
    pwms<-selectHomerBody(table)
    names(pwms)<-names
    ranks<-findRank(header)
    homers<-mapply(homer,header,pwms,ranks,SIMPLIFY=FALSE)
                                        #names(homers)<-sapply(homers,function(x) x$header$name
    #print(lapply(homers,function(x) x$rank))
    #print(addMotifHeader.homer(homers))
    class(homers)<-"motifList"
    homers
}

findRank<-function(list){
    if (all(sapply(list,function(x) !is.null(list$pvalue)))){
        ranks<-rank(sapply(list,function(x) x$pvalue),ties.method="first")
    }
    else {
        ranks<-seq(length(list))
    }
    ranks
}

#' @export
homer<-function(header,pwm,rank=0,contents=names(header)){
    ret<-list(header=header,pwm=pwm,contents=contents,
              pvalues=header$pvalue,name=header$name,score=header$score,
              motif=pwm2consensus(t(pwm)),
              pstring=header$pstring,
              rank=rank
              )
    class(ret)<-"homer"
    return(ret)
}


#' @export
getPwm.homer<-function(x,...){
    x$pwm
}

#' @export
getCons.homer<-function(x,...){
    x$motif
}

#' @export
getName.homer<-function(x,...){
    x$name
}

#' @export
homer2matrix<-function(x){
    do.call(rbind,lapply(x,as.list.homer))    
}

as.list.homer<-function(x){
    list(getName.homer(x),
    paste(x$header,collapse=","),
    getPwm.homer(x))
}
#' @export
mergeMotifLists<-function(...){
    x<-as.list(...)
    ret<-as.list(do.call(c,x))
    names(ret)<-sapply(ret,function(x) x$names)
    class(ret)<-"motifList"
    ret
}

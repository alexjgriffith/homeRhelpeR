checkHomerContents<-function(conts){    
    required<-c("consensus","name","threshold")
    optional<-c("lpvalue","gapped","occurence","stats")
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
read.homer<-function(fileLocation,contents=NULL)
{
    
    selectHomerHeader<-function(table,contents=NULL){
        whichH<-substr(table[,1],1,1)==">"
        header<-NULL
        if (is.null(contents)){
            header<-table[whichH,c(1,2,3)]
            colnames(header)<-c("consensus","name","threshold")
            header[,1]<-sapply(header[,1],function(x){
                substr(as.character(x),2,nchar(as.character(x)))})
        }
        else if (length(contents)==dim(table)[2]){
            contents<-checkHomerContents(contents)
            header<-table[whichH,]
            colnames(header)<-contents
        }
        else {
            warning("Inconsistancy between header width and contents passed into read.homer")
            contents<-checkHomerContents(contents)
            header<-table[whichH,1:length(contents)]
            colnames(header)<-contents
        }
        apply(header,1,as.list)
        
    }
    selectHomerBody<-function(table){
        whichH<-which(substr(table[,1],1,1)==">")
        pwmRanges<-cbind(whichH+1,c(whichH[2:length(whichH)]-1,dim(table)[1]))
        pwms<-apply(pwmRanges,1,function(x) table[x[1]:x[2],c(1,2,3,4)])
        lapply(pwms,function(pwm) as.matrix(apply(pwm,2,as.numeric)))
    }

    table<-read.table("motifs.motifs",fill=TRUE)
    header<-selectHomerHeader(table,contents)
    names<-sapply(header,function(x) x$name)
    names(header)<-names
    contents<-checkHomerContents(contents)
    pwms<-selectHomerBody(table)
    names(pwms)<-names
    homers<-mapply(homer,header,pwms,SIMPLIFY=FALSE,MoreArgs=list(contents))
                                        #names(homers)<-sapply(homers,function(x) x$header$name
    homers
}

homer<-function(header,pwm,contents=NULL){
    ret<-list(header=header,pwm=pwm,contents=contents)
    class(ret)<-"homer"
    return(ret)
}


getPwm.homer<-function(x,...){
    x$pwm
}

getCons.homer<-function(x,...){
    x$header$concensus
}

getName.homer<-function(x,...){
    x$header$name
}

homer2matrix<-function(x){
    do.call(rbind,lapply(x,as.list.homer))    
}

as.list.homer<-function(x){
    list(getName.homer(x),
    paste(x$header,collapse=","),
    getPwm.homer(x))
}


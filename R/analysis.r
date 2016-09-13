## combineClusterdMotifs
## mergeBase
## mergeConsensus
## consensus2pwm
## addMotifHeader
## motifList2Homer
## scorePWMHomer
## count
## countMotifs
## countMotifsOne
## consensus2pwm
## counts2PWM
## count2freq
## homerProbs
## homerFind
## addMotifHeader.homer
## motifs2Homer
## write.homer
## print.homer
## trimMotif

#' @export
trimMotif<-function(motifs){
    gsub("^N*|N*$","",motifs)
}

#' @export
clusterMotifs<-function(consMots){
    require(Biostrings)
    comb<-t(combn(consMots,2))
    as<-apply(comb,1,function(x) pairwiseAlignment(x[1],x[2],scoreOnly=TRUE))
    Xs<-matrix(0,ncol=length(consMots),nrow=length(consMots))
    Xs[lower.tri(Xs)]<-as
    rownames(Xs)<-consMots
    colnames(Xs)<-consMots
    as.dist(Xs)
}

#' @export
combineClusterdMotifs<-function(distance,names,k=length(names)/2){
    require(msa)
    ns<-unique(names)
    if(length(ns)!=length(names))
        stop("names must share dimensions with distance")
    clust<-hclust(distance,method="ward.D2")
    regs<-cutree(clust,k=k)
    sapply(seq(k),
           function(i){
               if(length(ns[regs==i])>1){                   
                   log<-capture.output({samp<-consensusMatrix(msa(ns[regs==i],type="dna"))})
                   mergeConsensus(samp)
                   
               }
               else
                   ns[regs==i]
               })
}


mergeBase<-function(string){
    temp<-paste0(sort(unique(Filter(function(x) !x %in% c("[","]",".","+","-"),strsplit(string,"")[[1]]))),collapse="")
    if(nchar(temp)==1)
        temp
    else
        paste0("[",temp,"]")
}

mergeConsensus<-function(samp){
    require(msa)
    sampNames<-rownames(samp)
    consensusIUPAC(paste(sapply(seq(dim(samp)[2]),function(i) mergeBase(IUPACtoBase(paste0(Filter(function(x) !x %in% c("[","]",".","+","-"),sampNames[samp[,i]>0]),collapse="")))),collapse=""))
}

#' @export
consensus2pwm<-function(string){
    chars<-strsplit(string,"")[[1]]
    n<-length(chars)
    do.call(rbind,lapply(chars,function(x){
        a<-matrix(rep(0.001,4),ncol=4)
        colnames(a)<-c("A","C","T","G")
        i<-Filter(function(i) !i %in% c("]","["), strsplit(IUPACtoBase(x),"")[[1]])
        a[,i]<-1/length(i)-(4-length(i))*0.001
        a
    }))                      
}


addMotifHeader<-function(motif,pwm){
    paste0(">",motif,"\t","0_",motif,"\t",100,"\t",-20,"\n",paste(apply(pwm,1,paste,collapse=" "),collapse="\n"),"\n")
}

#' @export
motifList2Homer<-function(motif){
    paste(sapply(motif,function(m)
        addMotifHeader(m,consensus2pwm(m))),collapse="")
}


#' @export
scoreMotifHomer<-function(motif){
    spl<-strsplit(motif,"")[[1]]
    sum(sapply(spl,function(x) log(0.997/(0.25*length(IUPACtoBase(x,rl=TRUE))))))
    
}

count<-function(list,ops=sort(unique(list))){    
    ret<-sapply(ops,function(x) sum(list==x))
    names(ret)<-ops
    ret
}

countMotifs<-function(sequence,motifs){
    ops<-IUPACtoBase(motifs)
    oprs<-IUPACtoBase(motifs,TRUE)
    oprsrev<-sapply(oprs,complement)
    ret<-rbind(count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,ops)),unique)),oprs),
               count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,complement(ops))),unique)),oprsrev))
    rownames(ret)<-c("motif","complement")
    ret
}

countMotifsOne<-function(sequence,motifs){
    oprsfor<-IUPACtoBase(motifs,TRUE)
    oprsrev<-sapply(oprsfor,complement)
    ops<-paste(oprsrev,oprsfor,collapse="|",sep="|")
    sum(grepl(ops,sequence))
    
}

consensus2pwm<-function(string){
    chars<-strsplit(string,"")[[1]]
    n<-length(chars)
    do.call(rbind,lapply(chars,function(x){
        a<-matrix(rep(0.001,4),ncol=4)
        colnames(a)<-c("A","C","T","G")
        i<-Filter(function(i) !i %in% c("]","["), strsplit(IUPACtoBase(x),"")[[1]])
        a[,i]<-1/length(i)-(4-length(i))*0.001
        a
    }))                      
}

counts2PWM<-function(counts){
    weights<-colSums(counts)
    mat<-matrix(0,ncol=4,nrow=nchar(names(weights)[1]))
    colnames(mat)<-c("A","C","G","T")
    mapply(function(name,weight){
        cs<-strsplit(name,"")[[1]]

        for(i in seq(length(cs))){
            mat[i,cs[i]]<<-mat[i,cs[i]]+weight
                                 }
    },
           names(weights),weights)
    mat
}

count2freq<-function(mat){
    round(t(apply(mat,1,function(x){
        freq<-x/sum(x)
        nil<-freq==0
        freq[nil]<-0.001
        freq[which.max(freq)]<-freq[which.max(freq)]-0.001*sum(nil)
        freq
    })),3)
}

homerProbs<-function(a,b,na,nb,a1,b1){
    fgmcr<-do.call("/",as.list(rowSums(a)))
    bgmcr<-do.call("/",as.list(rowSums(b)))
    fgc<-sum(a1)
    bgc<-sum(b1)
    fgr<-fgc/na
    bgr<-bgc/nb
    
    list(fgc=fgc,fgr=fgr,fgt=na,bgc=bgc,bgr=bgr,bgt=nb,fgmcr=fgmcr,bgmcr=bgmcr)
        
}

homerFind<-function(fg,bg,motif){
    a<-countMotifs(fg,motif)
    b<-countMotifs(bg,motif)
    a1<-countMotifsOne(fg,motif)
    b1<-countMotifsOne(bg,motif)
    na<-length(fg) # mot and compl
    nb<-length(bg)
    probs<-homerProbs(a,b,na,nb,a1,b1)
    pwm<-as.matrix(count2freq(counts2PWM(a)))
    colnames(pwm)<-c("A","C","G","T")
    square<-rbind(c(probs$fgc,probs$fgt),c(probs$bgc,probs$bgt))
    score<-scoreMotifHomer(motif)
    test<-chisq.test(square)
    p.value<-test$p.value
    pstring=paste0("T:",probs$fgc,".0(",round(probs$fgr*100,2),"%)",
                   "B:",probs$bgc,".0(",round(probs$bgr*100,2),"%)",
                   "P:1e",floor(log(p.value,10)))
    header<-list(consensus=paste0(">",motif),name=motif,score=score,pvalue=log(p.value),gapped=0,pstring=pstring)
    list(header=header,pwm=pwm)
}

addMotifHeader.homer<-function(x,...){
    with(x,{
        
        paste0(paste(c(header,"\n"),collapse="\t"),paste(apply(pwm,1,paste,collapse="\t"),collapse="\n"),"\n")
    })
}

#' @export
motifs2Homer<-function(fg,bg,mots){
    summary<-lapply(mots,function(mot)homerFind(fg,bg,mot))    
    ranks<-rank(sapply(summary,function(x) x$header$pvalue),ties.method="first")
    pwm<-lapply(summary,function(x) x[["pwm"]])
    header<-lapply(summary,function(x) x[["header"]])
    contents<-c("consensus","name","score","pvalue","gapped","pstring")
    ret<-mapply(homer,header,pwm,ranks,SIMPLIFY=FALSE,MoreArgs=list(contents=contents))
    class(ret)<-"motifList"
    ret
}

#' @export
write.motifList<-function(x,file){
    cat(sapply(x,addMotifHeader.homer),file=file,sep="")

}

#' @export
print.motifList<-function(x,...){
                                        #lapply(x,print
    cat(sapply(x,addMotifHeader.homer),sep="")
}

#' @export
summary.motifList<-function(object,...){
    do.call(rbind,lapply(object,function(x)c(x[c("motif","rank")],exp(as.numeric(x["pvalues"])))))
}

#' @export
print.homer<-function(x,...){

        cat(paste0(addMotifHeader.homer(x)),sep="")
}


#' @export
homerStamp.default<-function(homerFile,dbs,...){
    homerMotifs<-read.homer(homerFile)
    ann<-lapply(dbs,function(db) stampWrapper(homerFile,db,...))
    pcasum<-do.call(rbind,lapply(homerMotifs,function(x)c(x[c("name","motif","rank")],exp(as.numeric(x["pvalues"])))))
    
    colnames(pcasum)<-c("name","motif","rank","pvalues")

    cbind(do.call(rbind,ann),name2=unlist(lapply(pcasum[,"name"],rep,5)),
          tmotif=unlist(lapply(pcasum[,"motif"],rep,5)),
          order=unlist(lapply(pcasum[,"rank"],rep,5)),
          pvalues=unlist(lapply(pcasum[,"pvalues"],rep,5))
          )
}

#' @export
#'
homerStamp.motifList<-function(homerFile,dbs,...){
    homerLoc<-tempfile()
    write.motifList(homerFile,homerLoc)
    homerStamp.default(homerLoc,dbs,...)
}

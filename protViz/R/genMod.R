#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/genMod.R $
# $Id: genMod.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $


genMod<-function(sequences, modificationPattern, nModification=2){

    R<-list()
    for (peptideSequence in sequences){
        n<-nchar(peptideSequence)
        m<-length(modificationPattern)
        idx<-1:n
        pidx<-1:m
        acids<-(strsplit(peptideSequence, ''))[[1]]

        hits<-idx[acids %in% modificationPattern] 

        # no modification is always part of the result
        result<-paste(rep(0,n), collapse="")

        if (length(hits) < nModification)
            nModification <- length(hits)

        if (length(hits) > 1){
            for (ii in 1:nModification){
                c<-combn(hits, ii)

                result<-c(result, (apply(c, 2, 
                    FUN=function(x){
                        cand=rep(0,n); 
                        cand[x] <- pidx[ modificationPattern %in% acids[x] ]  - 1
                        return(paste(cand, collapse=""))
                    }) ))
                }
        }else if (length(hits) == 1){
            cand=rep(0,n);
            cand[hits]<-1
            result<-c(result, paste(cand, collapse=""))
        }

        R[[length(R)+1]] <- unique(result)
    }
    return(R)
}

#peptideSequence<-'HTLMQIDSVK'
#ptm.0<-list(pattern="xxx",
#    mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)
#
#ptm.STY<-list(pattern="[STY]", 
#    mono=79.966331, avg=79.9799, desc="Phosphorylation", 
#        unimodAccID=NA)
#
#ptm.M<-list(pattern="[M]", 
#    mono=15.994915, avg=NA, desc="Oxidation or Hydroxylation", 
#        unimodAccID=35)
#
#my.modification=as.data.frame(rbind(ptm.0, ptm.STY, ptm.M))$pattern
#
#genMod(peptideSequence, my.modification[2:3])
# tt<-genMod(c<-combn(c(1,4,6),2), n=10)

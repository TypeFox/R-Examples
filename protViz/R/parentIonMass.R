#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/parentIonMass.R $
# $Id: parentIonMass.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $

#fetuinPeptides<-c('LCPGR', 'GSVIQK', 'QYGFCK', 'AHYDLR', 'EVVDPTK', 'CNLLAEK', 'ALGGEDVR', 'HTLNQIDSVK', 'QDGQFSVLFTK', 'CDSSPDSAEDVR', 'TPIVGQPSIPGGPVR', 'LCPDCPLLAPLNDSR', 'QQTQHAVEGDCDIHVLK', 'HTFSGVASVESSSGEAFHVGK', 'EPACDDPDTEQAALAAVDYINK', 'AQFVPLPVSVSVEFAVAATDCIAK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK', 'VVHAVEVALATFNAESNGSYLQLVEISR', 'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR', 'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR')

parentIonMass <- function(sequence, fixmod=NULL, NTerm=1.007825) {
    if (!is.character(sequence)) {
        stop ("argument x must be a character")
    }else if (length(fixmod)>0 & length(fixmod)!=26){
        stop ("length of fixmod must be 26. One value for each letter from 1 to 26.")
    }else if (length(fixmod)==26){
           out <- .C("computeParentIonMass2",
                   n=as.integer(length(sequence)),
                   pepSeq=as.character(sequence),
                   pim=as.double(rep(0,length(sequence))),
                   M_=as.double(fixmod),
                   M_term_=as.double(NTerm)
                   )
    }else{
           out <- .C("computeParentIonMass",
                   n=as.integer(length(sequence)),
                   pepSeq=as.character(sequence),
                   pim=as.double(rep(0,length(sequence)))
                   )
    }
    return(out$pim)
}

# Christian Trachsel
# Christian Panse <cp@fgcz.ethz.ch>
# 2013-03-22 1mass.min
# takes as input a vector of peptide sequences and computes 
# 2+,3+,4+ mass histograms

.gasPhaseFragQuantiles<-function(peptides, mass.min=300, mass.max=1250, fixmod=NULL, NTerm=1.007825){

    pim <- parentIonMass(peptides, fixmod, NTerm)
    pim2 <- (pim + 1.00782)/2
    pim3 <- (pim + 2*1.00782)/3
    pim4 <- (pim + 3*1.00782)/4


    # breaks for Histogram
    bb<-seq(0,max(pim), by=10);

    hist(pim2, bb, xlim=c(mass.min,mass.max), prob=FALSE, xlab='mass bin size 10 Dalton')
    hist(pim3, bb, xlim=c(mass.min,mass.max), prob=FALSE, xlab='mass bin size 10 Dalton')
    hist(pim4, bb, xlim=c(mass.min,mass.max), prob=FALSE, xlab='mass bin size 10 Dalton')

    hist(c(pim2, pim3, pim4), bb, xlim=c(mass.min,mass.max), prob=FALSE)

    ss<-(sort(c(pim2, pim3, pim4)))

    sss<-ss[mass.min < ss & ss < mass.max]

    q <- quantile(sss)

    hist(sss, breaks=q, prob=FALSE, xlab="mass")

    axis(3, q, round(q))

    # not covert % by mass range
    100 * (length(ss[ ! ((mass.min < pim2 & pim2 < mass.max) 
        | (mass.min < pim3 & pim3 < mass.max) 
        | (mass.min < pim4 & pim4 < mass.max)) ]) / length(peptides))

}

# s<-read.table("p1239_sequences.txt")
# pdf("p1239-pep-charge-mass-hist.pdf", 19,12)
# par(mfrow=c(2,3),mar=c(6,4,6,1)); 
# gasPhaseFragQuantiles(as.character(s$V1))
# dev.off()


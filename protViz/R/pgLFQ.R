#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/pgLFQ.R $
# $Id: pgLFQ.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $

pgLFQtNpq <- function(QuantitativeValue, 
    peptide, protein, 
    N=3, plot=TRUE, FUN=asinh){

# to determine the top N; NOTE: wont be used for quantification
    s<-as.data.frame(cbind(Q=rowSums(QuantitativeValue),
        peptide=peptide,
        protein=protein, id=1:length(peptide)))

    r<-numeric()
    r.protein<-character()

# todo replace by lapply
    for (prot in unique(s$protein) ){
        qq<-(s[s$prot == prot,])
        qqq<-tapply(as.double(as.character(qq$Q)), paste(qq$peptide), sum)

        if (length(qqq) >= N){
            idx<-rev(order(qqq))[1:N]
                id<-as.integer(as.character(s[s$peptide %in% names(qqq[idx]), ]$id))
                r <- rbind(r, as.vector(as.double(colSums(QuantitativeValue[id,]))))
                r.protein<-c(r.protein, prot)
        }
    }

    if (plot)
        image(cor(FUN(r)),col=c(gray.colors(100), 'white'), 
            main=paste("top",N, "/ #Qprot", length(r.protein)), 
            axes='F') 

    row.names(r)=r.protein
    return (r)
}


pgLFQaov<-function(data, groups, names, idx=1:nrow(data), 
    plot=FALSE, 
    FUN=function(x){return(x)}){

    a<-rep(NA, nrow(data))

    for (i in idx){
        x.aov<-aov(FUN(as.numeric(data[i,])) ~ groups)
        a[i]<-round(summary(x.aov)[[1]][["Pr(>F)"]][1],4)

        boxplot.color='#ffcccc'
        if (a[i] < 0.05)
            boxplot.color='lightgreen'

        if (plot == TRUE){
            boxplot(FUN(as.numeric(data[i,]))~groups, 
                col=boxplot.color,
                sub=paste("ANOVA Pr(>F):",a[i]),
                main=names[i])
         }
    }
    return(a)
}


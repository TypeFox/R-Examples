#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/iTRAQ2GroupAnalysis.R $
# $Id: iTRAQ2GroupAnalysis.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $


iTRAQ2GroupAnalysis<-function(data, 
    group1, 
    group2, 
    INDEX,
    FUN=function(x){return(x)},
    plot=FALSE){

    g1<-numeric()
    g2<-numeric()

    for (i in group1){
        #g1 <- cbind(g1, asinh(tapply(data[,i], INDEX, sum, na.rm=TRUE)))
        g1 <- cbind(g1, FUN(tapply(data[,i], INDEX, sum, na.rm=TRUE)))

    }

    for (i in group2){
        #g2 <- cbind(g2, asinh(tapply(data[,i], INDEX, sum, na.rm=TRUE)))
        g2 <- cbind(g2, FUN(tapply(data[,i], INDEX, sum, na.rm=TRUE)))
    }


    R<-numeric()
    for (i in 1:nrow(g1)){
        boxplot.color='#ffcccc'
        tt.p_value<-t.test(as.numeric(g1[i,]) , as.numeric(g2[i,]))$p.value
        if (tt.p_value < 0.05)
            boxplot.color='lightgreen'

        if (row.names(g1)[i] != row.names(g2)[i])
            boxplot.color='cyan'

        if (plot == TRUE){
            b<-boxplot(as.numeric(g1[i,]), as.numeric(g2[i,]),
                main=paste(row.names(g1)[i]),
                sub=paste("t-Test: p-value =", round(tt.p_value,2)),
                col=boxplot.color)
            box()
            axis(1, 1:2, c('1','2'))

            points(rep(1,b$n[1]), as.numeric(g1[i,]), col='blue')
            points(rep(2,b$n[2]), as.numeric(g2[i,]), col='blue')
        }
        R<-rbind(R, cbind(name=row.names(g1)[i], p_value=as.numeric(round(tt.p_value,3)), Group1=data[i,group1], Group2=data[i,group2]))
    }
    return (as.data.frame(R))
}
# qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, group1=c(3,4,5,6), group2=7:10, INDEX=paste(iTRAQ$prot,iTRAQ$peptide), plot=FALSE)


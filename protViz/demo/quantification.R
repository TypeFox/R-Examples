#R
library(lattice)
data(fetuinLFQ)

cv<-1-1:7/10
t<-trellis.par.get("strip.background")
t$col<-(rgb(cv,cv,cv))
trellis.par.set("strip.background",t)

my.xlab="Fetuin concentration spiked into experiment [fmol]"
my.ylab<-"Abundance"

xyplot(abundance~conc|prot*method, data=fetuinLFQ$apex, groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)",
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,data=fetuinLFQ$empai,groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)" ,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,data=fetuinLFQ$t3pq,groups=prot,
    aspect=1,
    main="package:protViz - data(fetuinLFQ)" ,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
        if (groups[subscripts][1] == "Fetuin")  {
            panel.text(min(fetuin.t3pq$conc),
                max(fetuin.t3pq$abundance),
                paste("R-squared:", 
                round(summary(lm(x~y))$r.squared,2)),
                cex=0.75,
                pos=4)
        }
            },
    xlab=my.xlab,
    ylab=my.ylab
)

# iTRAQ
data(iTRAQ)
barchart(pl<-unique(iTRAQ[,1:2]),
    main="package:protViz - data(iTRAQ) / Overview" ,
    col='grey',
    xlab='number of unique peptides',
)

par(mfrow=c(4,4),mar=c(5,5,5,5));
for (i in 3:10){
    hist(asinh(iTRAQ[,i]),
    main=names(iTRAQ)[i])
            qqnorm(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
            qqline(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
}

par(mfrow=c(8,8),mar=c(1,1,1,1))
for (i in 3:10){
    for (j in 3:10){
        if (i==j){
            hist(asinh(iTRAQ[,i]), main=names(iTRAQ)[i])
        }
        else if (i<j){
            qqplot(asinh(iTRAQ[,i]) , asinh(iTRAQ[,j]),axes=F);box()
        } else{plot(0,0,type='n',axes=FALSE, xlab="",ylab="")}
        
    }
}

par(mfrow=c(1,1),mar=c(5,5,5,5));
b<-boxplot(asinh(iTRAQ[,c(3:10)]), axes=FALSE,
    main="package:protViz - data(iTRAQ) / QC" )
box()
axis(1,3:10-2,names(iTRAQ)[3:10])
axis(2)


par(mfrow=c(1,5),mar=c(6,5,5,1))
qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, group1=c(3,4,5,6), group2=7:10, INDEX=iTRAQ$prot, plot=TRUE, FUN=asinh)

op<-par(mfrow=c(3,5), mar=c(5,4,5,1))
qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, group1=c(3,4,5,6), group2=7:10, INDEX=paste(iTRAQ$prot,iTRAQ$peptide), plot=TRUE,FUN=asinh)
par(op)


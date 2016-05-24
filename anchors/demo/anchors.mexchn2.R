#######################################################################
## Replication of Figure 2 in King et al 2004
## 
## Author  : Jonathan Wand 
## Created : 2002-08-01 
##
## Modified
## - 2008-05-09 : JW
##   anchors 3.0 syntax
##
#######################################################################
cat("Replication of Figure 2 in King et al 2004\n")

data(mexchn)
dim(mexchn)

fo <- list(self = xsayself ~ 1,
           vign = cbind( xsay5,xsay4,xsay3,xsay2,xsay1) ~ 1)

ra  <- anchors( fo, data=mexchn, method="C")
raC <- anchors( fo, data=mexchn, subset = china == 1, method="C")
raM <- anchors( fo, data=mexchn, subset = china == 0, method="C")
summary(raC)
summary(raM)

mexchn.cut <- trim.data( mexchn, ra)
dim(mexchn.cut)

tb <- xtabs( ~  china + xsayself, data=mexchn.cut)
tbp <- tb / apply(tb,1,sum)

op <- par(no.readonly=TRUE)
par(mfrow=c(1,2))

barplot( tbp, beside=TRUE, col=c("grey","white"),xlab="y",ylim=c(0,.5),cex.names=.8,
        main="Raw responses: xsayself")

barplot(raM, raC, ties="uniform", col=c("grey","white"),xlab="C",ylim=c(0,.5),cex.names=.8,
     main="C (Ties Uniformly Allocated)")

par(op)

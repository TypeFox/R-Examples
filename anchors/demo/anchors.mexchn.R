#######################################################################
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2006-11-03
##
## Modified
## - 2008-05-09 : JW
##   anchors 3.0 syntax
##
#######################################################################
cat("Repl from King and Wand (2007) Fig1\n")

data(mexchn)
par(mfrow=c(2,2))

## these equations will be used throughout:
fo <- list(self = xsayself ~ 1,
           vign = cbind( xsay5,xsay4,xsay3,xsay2,xsay1) ~ 1,
           cpolr= ~ age + male + educyrs)

raC <- anchors( fo, data=mexchn, subset = china == 1, method="C")
raM <- anchors( fo, data=mexchn, subset = china == 0, method="C")
summary(raC)
summary(raM)

################### anchors
## Fig 1a
barplot(raC, raM, ties="uniform", col=c("grey","white"),xlab="C",ylim=c(0,.5),cex.names=.8,
     main="Uniformly Allocated")



## Fig 1b
umC <- fitted( raC, ties="uniform", average=TRUE)
umM <- fitted( raM, ties="uniform", average=TRUE)
um  <- rbind(umC, umM)

umc <- cbind(um[,1],um[, 2],
             um[,3]+um[, 4], 
             um[,5]+um[, 6],
             um[,7]+um[, 8],
             um[,9]+um[,10], 
             um[,11])

ump <- cbind(um[,1],
             um[,2],
             um[,3],
             um[,5],
             um[,7],
             um[,9],
             um[,11])/umc

bb <- barplot( umc,beside=TRUE,col=c("grey","white"),xlab="C",axes=TRUE,las=1,ylim=c(0,.5),cex.names=.8,
        names.arg=c("1","2","3,4","5,6","7,8","9,10","11"), main="Uniformly Allocated, Collapsed")
for (i in 3:(ncol(bb)-1)) {
  offset1 <- ump[1,i]-.5
  lines( c( bb[1,i]+offset1, bb[1,i]+offset1), c(0,umc[1,i]))
  offset2 <- ump[2,i]-.5
  lines( c( bb[2,i]+offset2, bb[2,i]+offset2), c(0,umc[2,i]))
}

################### cpolr
## note this differs from King and Wand by having separate cpolr models
## for each country, but uses the same covariates (King and Wand only allow for a mean
## shift by country)

## fig 1c
barplot(raC, raM, ties="cpolr", col=c("grey","white"),xlab="C",ylim=c(0,.5),cex.names=.8,
     main="Censored Ordered Probit")

umC <- fitted( raC, ties="cpolr", average=TRUE)
umM <- fitted( raM, ties="cpolr", average=TRUE)
um  <- rbind(umC, umM)

## Fig 1d
umc <- cbind(um[,1],um[, 2],
             um[,3]+um[, 4], 
             um[,5]+um[, 6],
             um[,7]+um[, 8],
             um[,9]+um[,10], 
             um[,11])

ump <- cbind(um[,1],
             um[,2],
             um[,3],
             um[,5],
             um[,7],
             um[,9],
             um[,11])/umc

bb <- barplot( umc,beside=TRUE,col=c("grey","white"),xlab="C",axes=TRUE,las=1,ylim=c(0,.5),cex.names=.8,
        names.arg=c("1","2","3,4","5,6","7,8","9,10","11"),main="Censored Ordered Probit, Collapsed")
for (i in 3:(ncol(bb)-1)) {
  offset1 <- ump[1,i]-.5
  lines( c( bb[1,i]+offset1, bb[1,i]+offset1), c(0,umc[1,i]))
  offset2 <- ump[2,i]-.5
  lines( c( bb[2,i]+offset2, bb[2,i]+offset2), c(0,umc[2,i]))
}

# ## 
# dev.off()
# system(paste("epstopdf",fname))

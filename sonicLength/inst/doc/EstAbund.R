### R code from vignette source 'EstAbund.Rnw'

###################################################
### code chunk number 1: makeUpData
###################################################
require(sonicLength)
set.seed(123)
chr <- sample(c(1:23,"X","Y"), 200, repl=TRUE )
pos <- sample(100000L,200)
strand <- sample(c("+","-"),200,repl=TRUE)
lens <- lapply(1:200, rfrag )
nlens <- sapply( lens, length )
loc.dframe <- data.frame( Chromosome=chr, Position=pos, Ort=strand )
len.dframe <- unique(cbind( loc.dframe[ rep( 1:200, nlens ) , ], length=unlist(lens) ))
rbind( head( len.dframe ), tail( len.dframe ) )



###################################################
### code chunk number 2: plotUnique
###################################################

id <- with( len.dframe , paste(Chromosome, Position, Ort ) )
id.counts <- table(factor(id,unique(id))) 
plot( as.vector(nlens), as.vector(id.counts),
     xlab='number of sonicants',ylab='distinct lengths',
     xlim=c(1,100), ylim=c(1,100))
abline(a=0,b=1,col='gray')


###################################################
### code chunk number 3: strOptions
###################################################

options( str = strOptions(strict.width="wrap") )



###################################################
### code chunk number 4: loadSL
###################################################

fit <- estAbund(id, len.dframe$length)

str( fit )

plot(nlens, fit$theta[unique(id)],
     xlab="actual sonicant count",
     ylab="estimated sonicant count"
     )
abline(a=0,b=1,col='gray')



###################################################
### code chunk number 5: repSim
###################################################

len.dframe$repl <- 1
lens2 <- lapply(1:200,rfrag, rate=0.03 )
nlens2 <- sapply( lens2, length )
len.dframe2 <- unique(cbind( loc.dframe[ rep( 1:200, nlens2 ) , ], length=unlist(lens2), repl=2 ))
lens3 <- lapply(1:200,rfrag, rate=0.04 )
nlens3 <- sapply( lens, length )
len.dframe3 <- unique(cbind( loc.dframe[ rep( 1:200, nlens3 ) , ], length=unlist(lens3), repl=3 ))

len.dframe <- rbind(len.dframe,len.dframe2,len.dframe3)

fit2 <- with(len.dframe, estAbund( paste(Chromosome, Position, Ort ),
                                   length, repl ))




###################################################
### code chunk number 6: strfit2
###################################################

str( fit2 )



###################################################
### code chunk number 7: phiPlot
###################################################

with( fit2, 
     plot( lframe$x[ lframe$orig ], phi, pch=16,col=lframe$strata[ lframe$orig ],
          xlab="length",
          ylab='phi',
          xlim=c(1,200)
          ))
legend("topright",col=1:3,pch=rep(16,3),legend=paste("replicate",1:3))




###################################################
### code chunk number 8: sonicPlot
###################################################

with( fit2,     
     plot(3*nlens , theta[unique(id)],
     xlab="actual sonicant count",
     ylab="estimated sonicant count"
     )
     )

abline(a=0,b=1,col='gray')




###################################################
### code chunk number 9: jackSim
###################################################

fit2 <- with(len.dframe, 
             estAbund( paste(Chromosome, Position, Ort ),
                      length, repl, jackknife=TRUE,
                      theta.var=TRUE))

head.names <- function(x) head( names(x) )

sapply( fit2, head.names )



###################################################
### code chunk number 10: jackLen
###################################################

length(fit2$jackknife)



###################################################
### code chunk number 11: jackStr
###################################################

str(fit2$jackknife[[1]])



###################################################
### code chunk number 12: maxprop
###################################################
nreps <- length( fit2$jackknife )
argmax.theta <- which.max( fit2$theta )
## function to extract the estimate
maxpr <- function(x) prop.table( x$theta )[ argmax.theta ]
## apply to full sample
maxpr.total <- maxpr( fit2 )
## make pseudo observations
pseudo.maxpr <- nreps * maxpr.total - 
  (nreps-1) * sapply( fit2$jackknife, maxpr )  

## report results
likeStdErr <- sqrt( fit2$var.theta$prop[ argmax.theta ] )
jackStdErr = sd ( pseudo.maxpr ) / sqrt(nreps) 
all.res <- c( uncorrected = unname(maxpr.total), 
             corrected = mean(pseudo.maxpr),
             likeStdErr = likeStdErr,
             jackStdErr = jackStdErr)
cbind(all.res)
            



###################################################
### code chunk number 13: demoSimSonic
###################################################

more.data <- simSonic(fit2$theta,fit2$phi)

fit3 <- do.call(estAbund, more.data )

plot( fit2$theta[names(fit3$theta)], fit3$theta )
abline(a=0,b=1)

str(fit3)




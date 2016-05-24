### R code from vignette source 'both_papers.Rnw'

###################################################
### code chunk number 1: both_papers.Rnw:409-409
###################################################



###################################################
### code chunk number 2: both_papers.Rnw:410-413
###################################################
set.seed(0)
require(BACCO)
data(toy)


###################################################
### code chunk number 3: both_papers.Rnw:416-418
###################################################
toy <- latin.hypercube(20,6)
head(toy)


###################################################
### code chunk number 4: both_papers.Rnw:427-429
###################################################
f <- function(x) sum( (0:6)*x)
expectation <- apply(regressor.multi(toy), 1, f)


###################################################
### code chunk number 5: both_papers.Rnw:437-440
###################################################
toy.scales <- rep(1,6)
toy.sigma <- 0.4
A <- toy.sigma*corr.matrix(xold=toy, scales=toy.scales)


###################################################
### code chunk number 6: both_papers.Rnw:448-449
###################################################
d <-  as.vector(rmvnorm(n=1 , mean=expectation , sigma=A))


###################################################
### code chunk number 7: both_papers.Rnw:457-460
###################################################
x.unknown <- rep(0.5 , 6)
jj <- interpolant(x.unknown, d, toy, scales=toy.scales, g=TRUE)
print(drop(jj$mstar.star))


###################################################
### code chunk number 8: both_papers.Rnw:468-469
###################################################
print(jj$betahat)


###################################################
### code chunk number 9: both_papers.Rnw:473-474
###################################################
print(jj$beta.marginal.sd)


###################################################
### code chunk number 10: both_papers.Rnw:500-502
###################################################
scales.optim <- optimal.scales(val=toy, scales.start=rep(1,6), d=d, give=FALSE)
print(scales.optim)


###################################################
### code chunk number 11: both_papers.Rnw:515-517
###################################################
interpolant(x.unknown, d, toy, scales=toy.scales  , g=FALSE)
interpolant(x.unknown, d, toy, scales=scales.optim, g=FALSE)


###################################################
### code chunk number 12: both_papers.Rnw:564-639
###################################################

data(results.table) 
data(expert.estimates)
            # Decide which column we are interested in:
     output.col <- 26
     wanted.cols <- c(2:9,12:19)
            # Decide how many to keep;
            # 30-40 is about the most we can handle:
     wanted.row <- 1:27

            # Values to use are the ones that appear in goin.test2.comments:
     val <- results.table[wanted.row , wanted.cols]

            # Now normalize val so that 0<results.table[,i]<1 for all i:

     normalize <- function(x){(x-mins)/(maxes-mins)}
     unnormalize <- function(x){mins + (maxes-mins)*x}

     mins  <- expert.estimates$low 
     maxes <- expert.estimates$high
     jj <- t(apply(val,1,normalize))

     val <- as.matrix(jj)


            ## Answer is the 19th (or 20th or ... or 26th)
     d  <- results.table[wanted.row ,  output.col]

scales.optim <- 
  c(0.054095, 0.007055, 0.034944, 10.772536, 0.085691, 0.144568, 0.033540, 
    0.641465, 0.235039, 0.046189, 0.949328, 0.055576, 0.058894, 0.098077,
    0.045411, 0.167629)

     A <- corr.matrix(val,scales=rep(1,ncol(val)))
     Ainv <-  solve(A)

            ## Now try to find the best correlation lengths:
            ## the idea is to minimize the maximum absolute deviation
            ## from the known points.

d.observed <- results.table[, output.col]
A <- corr.matrix(val,scales=scales.optim)
Ainv <- solve(A)

design.normalized <- as.matrix(t(apply(results.table[,wanted.cols],1,normalize)))
jj.preds <- interpolant.quick(design.normalized, d, val, Ainv, give.Z=TRUE,
                                 scales=scales.optim)
d.predicted <- jj.preds$mstar.star
d.errors <- jj.preds$Z

jj <- range(c(d.observed,d.predicted))
par(pty="s")
plot(d.observed, d.predicted, pch=16, asp=1,
     xlim=jj,ylim=jj,
     xlab=expression(paste(temperature," (",{}^o,C,"), model"   )),
     ylab=expression(paste(temperature," (",{}^o,C,"), emulator"))
     )

errorbar <- function(x,y,delta,serifwidth=7,...){
  if(abs(delta)<1e-5){return()}
  lines(x=c(x,x),y=c(y-delta,y+delta), ...)
  lines(x=c(x-serifwidth/2,x+serifwidth/2),
        y=c(y+delta,y+delta), ...
        )
  lines(x=c(x-serifwidth/2,x+serifwidth/2),
        y=c(y-delta,y-delta), ...
        )
}
  
for(i in 1:length(d.observed)){
  errorbar(x=d.observed[i],y=d.predicted[i],delta=qt(0.975,df=11)*d.errors[i],serifwidth=0.1,col="red")
}

abline(0,1)



###################################################
### code chunk number 13: both_papers.Rnw:898-899
###################################################
args(ht.fun)


###################################################
### code chunk number 14: both_papers.Rnw:1249-1264
###################################################
load.the.files <- TRUE
#library(goldstein)
if(load.the.files){
  load("e10000")
  load("temps.jj")
  o <- order(temps.jj)
   
  load("probs.from.prior")
  j0 <- j0-min(j0)
}
x.pp <- cumsum(exp(j0[o] ))
plot(temps.jj[o],x.pp/max(x.pp),ylab="cumulative probability",xlab="temperature (C)",
     type="l",lty=1,main="Prior and posterior CDF for temperature in Northern Europe")
points(sort(temps.jj),(1:9000)/9000,type="l",lty=2)
legend(0,0.95,legend=c("posterior","prior"),lty=1:2)


###################################################
### code chunk number 15: SetTheBib
###################################################
bib <- system.file( "doc", "uncertainty.bib", package = "emulator" )
bib <- sub('.bib$','',bib)


###################################################
### code chunk number 16: usethebib
###################################################
cat( "\\bibliography{",bib,"}\n",sep='')



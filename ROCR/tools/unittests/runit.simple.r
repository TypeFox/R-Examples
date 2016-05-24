
library(RUnit)

library(ROCR)
#source("prediction.R")
#source("performance.R")
#source("performance_measures.R")
#source("zzz.R")
 
some.predictions <- c(0.02495517, 0.92535646,
                      0.86251887, 0.80946685,
                      0.70922858, 0.69762824,
                      0.50604485, 0.25446810,
                      0.10837728, 0.07250349)
some.labels <- c(0,1,1,0,1,1,0,1,0,0)

tp.reference <- c(0, 1, 2, 2, 3, 4, 4, 5, 5, 5, 5)
fp.reference <- c(0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5)

pp.reference <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
np.reference <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
p.reference <- rep(5, 11)
n.reference <- rep(5, 11)

tn.reference <- n.reference-fp.reference
fn.reference <- p.reference-tp.reference

# manually calculated reference measures
rpp.reference <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
rnp.reference <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)

tpr.reference <- c(0.0, 0.2, 0.4, 0.4, 0.6, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0)
fpr.reference <- c(0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.4, 0.4, 0.6, 0.8, 1.0)
acc.reference <- c(0.5, 0.6, 0.7, 0.6, 0.7, 0.8, 0.7, 0.8, 0.7, 0.6, 0.5)
err.reference <- c(0.5, 0.4, 0.3, 0.4, 0.3, 0.2, 0.3, 0.2, 0.3, 0.4, 0.5)

rec.reference <- tpr.reference
sens.reference<- tpr.reference
fnr.reference <- c(1.0, 0.8, 0.6, 0.6, 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0) 
tnr.reference <- c(1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.6, 0.6, 0.4, 0.2, 0.0)
spec.reference<- tnr.reference

ppv.reference <- c(0/0, 1/1, 2/2, 2/3, 3/4, 4/5, 4/6, 5/7, 5/8, 5/9, 5/10)
npv.reference <- c(5/10, 5/9, 5/8, 4/7, 4/6, 4/5, 3/4, 3/3, 2/2, 1/1, 0/0)
prec.reference<- ppv.reference

fall.reference <- fpr.reference
miss.reference <- fnr.reference
pcfall.reference <- c(0/0, 0/1, 0/2, 1/3, 1/4, 1/5, 2/6, 2/7, 3/8, 4/9, 5/10)
pcmiss.reference <- c(5/10, 4/9, 3/8, 3/7, 2/6, 1/5, 1/4, 0/3, 0/2, 0/1, 0/0)

auc.reference <- 0.84

cal.reference <- c()
ind <- rev(order(some.predictions))
sorted.predictions <- some.predictions[ind]
sorted.labels <- some.labels[ind]
for (i in 1:8) {
    mean.pred <- mean( sorted.predictions[i:(i+2)] )
    frac.pos <- sum( sorted.labels[i:(i+2)] ) / 3
    cal.reference <- c(cal.reference, abs( mean.pred - frac.pos ))
}

prbe.reference<- 0.8
prbe.reference.x <- 0.69762824
rch.reference.x <- fpr.reference[c(1,3,6,8,11)]
rch.reference.y <- tpr.reference[c(1,3,6,8,11)]
mxe.reference <- -(1/length(some.predictions)) * sum(some.labels*log(some.predictions) +
                                                     (1-some.labels)*log(1-some.predictions))
rmse.reference <- sqrt((1/length(some.predictions)) * sum((some.predictions-some.labels)^2))

phi.reference <- (tp.reference*tn.reference-fp.reference*fn.reference) /
                 sqrt(p.reference*n.reference*pp.reference*np.reference)
mat.reference <- phi.reference

my.log2 <- function( x ) {
    ans <- log2(x)
    ans[ ans==-Inf ] <- 0
    ans
}

mi.reference <- (tn.reference * my.log2( tn.reference / (n.reference*np.reference)) +
                 fn.reference*my.log2(fn.reference/(np.reference*p.reference)) +
                 fp.reference*my.log2(fp.reference/(n.reference*pp.reference)) +
                 tp.reference*my.log2(tp.reference/(p.reference*pp.reference))) / length(some.labels) + log2(length(some.labels))

chisq.reference <-
(((pp.reference*p.reference/length(some.predictions)) -
tp.reference)^2 / (pp.reference*p.reference/length(some.predictions))
+ ((pp.reference*n.reference/length(some.predictions)) -
fp.reference)^2 / (pp.reference*n.reference/length(some.predictions))
+ ((np.reference*p.reference/length(some.predictions)) -
fn.reference)^2 / (np.reference*p.reference/length(some.predictions))
+ ((np.reference*n.reference/length(some.predictions)) -
tn.reference)^2 / (np.reference*n.reference/length(some.predictions)))

odds.reference <- (tp.reference*tn.reference) / (fn.reference*fp.reference)
                                             
lift.reference <- (tp.reference/p.reference) / (pp.reference/(p.reference+n.reference))

f.reference <- 1 / (0.5 * ((1/prec.reference) + (1/rec.reference)))

sar.reference <- 1/3 * (acc.reference + auc.reference + (1-rmse.reference))
cost.reference <- (fpr.reference * n.reference/length(some.labels) * 1 +
                   fnr.reference * p.reference/length(some.labels) * 1)
  
.get.performance.measures <- function(pred) {
    
    tpr <- performance(pred, "tpr")@y.values[[1]]
    fpr <- performance(pred, "fpr")@y.values[[1]]
    acc <- performance(pred, "acc")@y.values[[1]]
    err <- performance(pred, "err")@y.values[[1]]

    rec <- performance(pred, "rec")@y.values[[1]]
    sens<- performance(pred, "sens")@y.values[[1]]
    fnr <- performance(pred, "fnr")@y.values[[1]]
    tnr <- performance(pred, "tnr")@y.values[[1]]
    spec<- performance(pred, "spec")@y.values[[1]]
    ppv <- performance(pred, "ppv")@y.values[[1]]
    prec<- performance(pred, "prec")@y.values[[1]]
    npv <- performance(pred, "npv")@y.values[[1]]
   
    fall<- performance(pred, "fall")@y.values[[1]]
    miss<- performance(pred, "miss")@y.values[[1]]
    pcfall <- performance(pred, "pcfall")@y.values[[1]]
    pcmiss <- performance(pred, "pcmiss")@y.values[[1]]
    rpp <- performance(pred, "rpp")@y.values[[1]]
    rnp <- performance(pred, "rnp")@y.values[[1]]
    
    auc <- performance(pred, "auc")@y.values[[1]]
    prbe<- performance(pred, "prbe")@y.values[[1]]
    rch <- performance(pred, "rch")@y.values[[1]]

    mxe <- performance(pred, "mxe")@y.values[[1]]
    rmse<- performance(pred, "rmse")@y.values[[1]]

    phi <- performance(pred, "phi")@y.values[[1]]
    mat <- performance(pred, "mat")@y.values[[1]]
    mi  <- performance(pred, "mi")@y.values[[1]]
    chisq<- performance(pred, "chisq")@y.values[[1]]
    odds<- performance(pred, "odds")@y.values[[1]]
    lift<- performance(pred, "lift")@y.values[[1]]
    f   <- performance(pred, "f")@y.values[[1]]
    sar <- performance(pred,"sar")@y.values[[1]]
    ecost  <- performance(pred, "ecost")@y.values[[1]]
    cost  <- performance(pred, "cost")@y.values[[1]]
    return(list(tpr=tpr, fpr=fpr, acc=acc, err=err,
                rec=rec, sens=sens, fnr=fnr, tnr=tnr,
                spec=spec, ppv=ppv, prec=prec, npv=npv, 
                fall=fall, miss=miss, pcfall=pcfall, pcmiss=pcmiss, rpp=rpp, rnp=rnp,
                auc=auc, prbe=prbe, rch=rch, mxe=mxe, 
                rmse=rmse, phi=phi, mat=mat, mi=mi, chisq=chisq, odds=odds,
                lift=lift, f=f, sar=sar, ecost=ecost, cost=cost))

}


testEcost <- function() {
    ecost.x.reference <- c(0,1/3,0.5,1)
    ecost.y.reference <- c(0,0.2,0.2,0)
    
    pred <- prediction(some.predictions, some.labels)
    perf <- performance(pred, "ecost")
    ecost.x <- perf@x.values[[1]]
    ecost.y <- perf@y.values[[1]]
    
    checkEquals( ecost.x, ecost.x.reference )
    
    checkEquals( ecost.y, ecost.y.reference )
}

testCal <- function() {
    pred <- prediction(some.predictions, some.labels)
    cal <- performance(pred, "cal", window.size=floor(length(pred@predictions[[1]])/3))@y.values[[1]]
    cal.x <- performance(pred, "cal", window.size=floor(length(pred@predictions[[1]])/3))@x.values[[1]]
    cal.x.reference <- rev(sort( some.predictions ))[2:(length(some.predictions)-1)]
    
    checkEquals( cal, cal.reference)
    checkEquals( cal.x, cal.x.reference)
}

testCost <- function() {
    pred <- prediction(some.predictions, some.labels)

    for (cost.fp in rnorm(50)) {
        cost.fn <- rnorm(1)
        
        perf <- performance(pred, "cost", cost.fp=cost.fp, cost.fn=cost.fn)
        cost <- perf@y.values[[1]]
        my.cost.reference <- (fpr.reference * n.reference/length(some.labels) * cost.fp +
                              fnr.reference * p.reference/length(some.labels) * cost.fn)

        checkEquals( cost, my.cost.reference)
    }
}

testRch <- function() {
    pred <- prediction(some.predictions, some.labels)
    perf <- performance( pred, "rch")
    rch.x <- perf@x.values[[1]]
    rch.y <- perf@y.values[[1]]

    checkEquals( rch.x, rch.reference.x )
    checkEquals( rch.y, rch.reference.y )
    
}

testPerformanceMeasuresReference <- function() {
    pred <- prediction(some.predictions, some.labels)
    measures <- .get.performance.measures(pred)
    attach(measures)

    checkEquals(tpr, tpr.reference)
    checkEquals(fpr, fpr.reference)
    checkEquals(acc, acc.reference)    
    checkEquals(err, err.reference)
    checkEquals(rec, rec.reference)
    checkEquals(sens, sens.reference)
    checkEquals(fnr, fnr.reference)
    checkEquals(tnr, tnr.reference)
    checkEquals(spec, spec.reference)

    checkEquals(ppv, ppv.reference)
    checkEquals(prec,prec.reference)
    checkEquals(npv, npv.reference)

    checkEquals(fall, fall.reference)
    checkEquals(miss,miss.reference)
    checkEquals(pcfall, pcfall.reference)
    checkEquals(pcmiss,pcmiss.reference)
    checkEquals(rpp, rpp.reference)
    checkEquals(rnp,rnp.reference)
    
    checkEquals(auc, auc.reference)
    checkEquals(prbe, prbe.reference)

    checkEquals(mxe, mxe.reference)

    checkEquals(rmse, rmse.reference)
    checkEquals(phi, phi.reference)
    checkEquals(mat, mat.reference)

    checkEquals(mi, mi.reference)

    checkEquals(chisq, chisq.reference)
    checkEquals(odds, odds.reference)
    checkEquals(lift, lift.reference)

    checkEquals(f, f.reference)
    checkEquals(sar,sar.reference)

    checkEquals(cost, cost.reference)
}

testRMSE <- function() {
    pred <- prediction(c(0, 0, 1, 1),
                       ordered(c(0, 0, 1, 1)))
    rmse <- performance(pred, "rmse")@y.values[[1]]
    checkEquals(rmse, 0)

    pred <- prediction(c(0.0, 0.0, 1.0, 1.0),
                       ordered(c(1, 1, 0, 0), levels=c(1,0)))
    rmse <- performance(pred, "rmse")@y.values[[1]]

    checkEquals(rmse, 1)

    pred <- prediction(c(0.0, 0.0, 1.0, 1.0),
                       ordered(c(2, 2, 3, 3)))
    rmse <- performance(pred, "rmse")@y.values[[1]]

    checkEquals( rmse, 2)

    pred <- prediction(c(-0.5, 0.2, 2.5, 0.3),
                       ordered(c(-1, -1, 1, 1)))
    rmse <- performance(pred, "rmse")@y.values[[1]]

    checkEquals( rmse, sqrt(1/4*(0.5^2 + 1.2^2 + 1.5^2 + 0.7^2)))
    
}

testPRBE <- function() {
    pred <- prediction(some.predictions, some.labels)
    prbe.y <- performance(pred, "prbe")@y.values[[1]]
    prbe.x <- performance(pred, "prbe")@x.values[[1]]
    checkEquals(prbe.y, prbe.reference)
    checkEquals(prbe.x, prbe.reference.x)
}
    
testPredictionInterface <- function() {
    pred <- prediction(seq(0, 1, length=10),
                       c(rep(0,5), rep(1,5)))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                1)    
    pred <- prediction(seq(1, 0, length=10),
                       c(rep(0,5), rep(1,5)))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                0)    
    pred <- prediction(seq(0, 1, length=10),
                       factor(c(rep(0,5), rep(1,5))))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                1)    
    pred <- prediction(seq(0, 1, length=10),
                       ordered(c(rep(0,5), rep(1,5))))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                1)
    pred <- prediction(seq(0, 1, length=10),
                       ordered(c(rep(0,5), rep(1,5)), levels=c(1,0)))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                0)    
    pred <- prediction(seq(0, 1, length=10),
                       ordered(c(rep("A",5), rep("B",5))))
    checkEquals(performance(pred, "auc")@y.values[[1]],
                1)

    checkException(pred <- prediction(seq(0, 1, length=10),
                                      c(rep(0,5), rep(1,5)), label.ordering=c(1,2)))
    checkException(pred <- prediction(list(c(0.1,0.3,0.7,1),
                                           c(0,0.2,0.8,1)),
                                      list(factor(c(0,0,1,1)),
                                           factor(c(1,1,2,2))))) 
    checkException(pred <- prediction(list(c(0.2,0.3,0.7,1),
                                           c(0,0.2,0.8,1)),
                                      list(factor(c(0,0,1,1)),
                                           ordered(c(0,0,1,1)))))
    pred <- prediction(list(c(0,0.3,0.7,1),
                            c(0,0.2,0.8,1)),
                       list(factor(c(0,0,1,1)),
                            factor(c(0,0,1,1))))
    checkEquals(performance(pred, "auc")@y.values,
                list(1, 1))

    pred1 <- prediction(data.frame(c(0,0.3,0.7,1),
                            c(0,0.2,0.8,1)),
                       data.frame(factor(c(0,0,1,1)),
                            factor(c(0,0,1,1))))
    checkEquals( pred, pred1)
    
    pred2 <- prediction(cbind(c(0,0.3,0.7,1),
                              c(0,0.2,0.8,1)),
                        cbind(c(0,0,1,1),
                              c(0,0,1,1)))
    checkEquals(pred, pred2)
                
    
}

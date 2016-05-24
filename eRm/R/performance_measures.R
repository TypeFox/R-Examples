## ------------------------------------------------------------------------
## classical machine learning contingency table measures
## ------------------------------------------------------------------------

.performance.accuracy <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      list( cutoffs, (tn+tp) / length(predictions) )
  }

.performance.error.rate <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      list( cutoffs, (fn+fp) / length(predictions) )
  }

.performance.false.positive.rate <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      list( cutoffs, fp / n.neg )
  }

.performance.true.positive.rate <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      list( cutoffs, tp / n.pos )
  }

.performance.false.negative.rate <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      list( cutoffs, fn / n.pos )
  }

.performance.true.negative.rate <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      list( cutoffs, tn / n.neg )
  }

.performance.positive.predictive.value <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      ppv <- tp / (fp + tp)
      list( cutoffs, ppv )
  }

.performance.negative.predictive.value <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
    
      npv <- tn / (tn + fn)
      list( cutoffs, npv )
  }

.performance.prediction.conditioned.fallout <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      ppv <- .performance.positive.predictive.value(predictions, labels,
                                                    cutoffs, fp, tp, fn, tn,
                                                    n.pos, n.neg, n.pos.pred,
                                                    n.neg.pred)[[2]]
      list( cutoffs, 1 - ppv )
  }

.performance.prediction.conditioned.miss <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      npv <- .performance.negative.predictive.value(predictions, labels,
                                                    cutoffs, fp, tp, fn, tn,
                                                    n.pos, n.neg, n.pos.pred,
                                                    n.neg.pred)[[2]]
      list( cutoffs, 1 - npv )
  }

## ------------------------------------------------------------------------
## ...not actually performance measures, but very useful as a second axis
## against which to plot a "real" performance measure
## (popular example: lift charts)
## ------------------------------------------------------------------------

.performance.rate.of.positive.predictions <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      list( cutoffs, n.pos.pred / (n.pos + n.neg) )
  }

.performance.rate.of.negative.predictions <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      list( cutoffs, n.neg.pred / (n.pos + n.neg) )
  }


## ------------------------------------------------------------------------
## Classical statistical contingency table measures
## ------------------------------------------------------------------------

.performance.phi <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      list(cutoffs,
           (tn*tp - fn*fp) / sqrt(n.pos * n.neg * n.pos.pred * n.neg.pred) )
  }

.performance.mutual.information <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      n.samples <- n.pos + n.neg
      mi <- c()
      for (k in 1:length(cutoffs)) {
          kij <- rbind( c(tn[k],fn[k]), c(fp[k],tp[k]) )

          ki.j. <- rbind(c(n.neg * n.neg.pred[k], n.neg.pred[k] * n.pos),
                         c(n.neg * n.pos.pred[k], n.pos * n.pos.pred[k]))

          log.matrix <- log2( kij / ki.j.)
          log.matrix[kij/ki.j.==0] <- 0
          
          mi <- c(mi,  log2(n.samples) + sum( kij * log.matrix) / n.samples  )
      }

      list( cutoffs, mi )
  }


.performance.chisq <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      chisq <- c()
      for (i in 1:length(cutoffs)) {
          A <- rbind( c( tn[i], fn[i]), c(fp[i], tp[i]) )
          chisq <- c(chisq, chisq.test(A,correct=F)$statistic )
      }
      list( cutoffs, chisq )
  }

.performance.odds.ratio <- 
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      

    list( cutoffs, tp * tn / (fn * fp) )
    
}

## ------------------------------------------------------------------------
## Other measures based on contingency tables
## ------------------------------------------------------------------------

.performance.lift <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      
      n.samples <- n.pos + n.neg
      list( cutoffs, (tp / n.pos) / (n.pos.pred / n.samples) )
  }

.performance.f <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred, alpha) {

      prec <- .performance.positive.predictive.value(predictions, labels,
                                                     cutoffs, fp, tp, fn, tn,
                                                     n.pos, n.neg, n.pos.pred,
                                                     n.neg.pred)[[2]]
      list( cutoffs,  1/ ( alpha*(1/prec) + (1-alpha)*(1/(tp/n.pos))  ) )
  }

.performance.rocconvexhull <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      
      x <- fp / n.neg
      y <- tp / n.pos

      finite.bool <- is.finite(x) & is.finite(y)
      x <- x[ finite.bool ]
      y <- y[ finite.bool ]
      if (length(x) < 2) {
          stop("Not enough distinct predictions to compute ROC convex hull.")
      }

      ## keep only points on the convex hull
      ind <- chull(x, y)
      x.ch <- x[ind]
      y.ch <- y[ind]

      ## keep only convex hull points above the diagonal, except (0,0)
      ## and (1,1)
      ind.upper.triangle <- x.ch < y.ch
      x.ch <- c(0, x.ch[ind.upper.triangle], 1)
      y.ch <- c(0, y.ch[ind.upper.triangle], 1)

      ## sort remaining points by ascending x value
      ind <- order(x.ch)
      x.ch <- x.ch[ind]
      y.ch <- y.ch[ind]

      list( x.ch, y.ch )
  }

## ----------------------------------------------------------------------------
## Cutoff-independent measures
## ----------------------------------------------------------------------------

.performance.auc <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred, fpr.stop) {
      
      x <- fp / n.neg
      y <- tp / n.pos

      finite.bool <- is.finite(x) & is.finite(y)
      x <- x[ finite.bool ]
      y <- y[ finite.bool ]
      if (length(x) < 2) {
          stop(paste("Not enough distinct predictions to compute area",
                     "under the ROC curve."))
      }

      if (fpr.stop < 1) {
        ind <- max(which( x <= fpr.stop ))
        tpr.stop <- approxfun( x[ind:(ind+1)], y[ind:(ind+1)] )(fpr.stop)
        x <- c(x[1:ind], fpr.stop)
        y <- c(y[1:ind], tpr.stop)
      }
      
      ans <- list()
      auc <- 0
      for (i in 2:length(x)) {
          auc <- auc + 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
      }
      ans <- list( c(), auc)
      names(ans) <- c("x.values","y.values")
      return(ans)
  }

.performance.precision.recall.break.even.point <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      pred <- prediction( predictions, labels)
      perf <- performance( pred, measure="prec", x.measure="rec")
      x <- rev(perf@x.values[[1]])
      y <- rev(perf@y.values[[1]])
      alpha <- rev(perf@alpha.values[[1]])

      finite.bool <- is.finite(alpha) & is.finite(x) & is.finite(y)
      x <- x[ finite.bool ]
      y <- y[ finite.bool ]
      alpha <- alpha[ finite.bool ]

      if (length(x) < 2) {
          stop(paste("Not enough distinct predictions to compute",
                     "precision/recall intersections."))
      }
      intersection.cutoff <- c()
      intersection.pr <- c()
      
      ## find all intersection points by looking at all intervals (i,i+1):
      ## if the difference function between x and y has different signs at the
      ## interval boundaries, then an intersection point is in the interval;
      ## compute as the root of the difference function
      if ( (x[1]-y[1]) == 0) {
          intersection.cutoff <- c( alpha[1] )
          intersection.pr <- c( x[1] )
      }

      for (i in (1:(length(alpha)-1))) {
          if ((x[i+1]-y[i+1]) == 0) {
              intersection.cutoff <- c( intersection.cutoff, alpha[i+1] )
              intersection.pr <- c( intersection.pr, x[i+1] )
          } else if ((x[i]-y[i])*(x[i+1]-y[i+1]) < 0 ) {
              ans <- uniroot(approxfun(c(alpha[i], alpha[i+1] ),
                                       c(x[i]-y[i], x[i+1]-y[i+1])),
                             c(alpha[i],alpha[i+1]))
              intersection.cutoff <- c(intersection.cutoff, ans$root)
              intersection.pr <- c(intersection.pr, ans$f.root)
          }
      }

      list( rev(intersection.cutoff), rev(intersection.pr) )
  }


.performance.calibration.error <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred, window.size) {

      if (window.size > length(predictions)) {
          stop("Window size exceeds number of predictions.")
      }
      if (min(predictions)<0 || max(predictions)>1) {
          stop("Calibration error needs predictions between 0 and 1")
      }
      
      pos.label <- levels(labels)[2]
      neg.label <- levels(labels)[1]

      ordering <- rev(order( predictions ))
      predictions <- predictions[ ordering ]
      labels <- labels[ ordering ]

      median.cutoffs <- c()
      calibration.errors <- c()

      for (left.index in 1 : (length(predictions) - window.size+1) ) {
          right.index <- left.index + window.size - 1
          pos.fraction <-
            sum(labels[left.index : right.index] == pos.label) / window.size
          mean.prediction <- mean( predictions[ left.index : right.index ] )

          calibration.errors <- c(calibration.errors,
                                  abs(pos.fraction - mean.prediction))
          median.cutoffs <- c(median.cutoffs,
                              median(predictions[left.index:right.index]))
      }
      list( median.cutoffs, calibration.errors )
  }


.performance.mean.cross.entropy <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

      if (! all(levels(labels)==c(0,1)) ||
          any(predictions<0) || any(predictions>1) ) {
          stop(paste("Class labels need to be 0 and 1 and predictions between",
                     "0 and 1 for mean cross entropy."))
      }
      
      pos.label <- levels(labels)[2]
      neg.label <- levels(labels)[1]
    
      list( c(), - 1/length(predictions) *
           (sum( log( predictions[which(labels==pos.label)] ))  +
            sum( log( 1 - predictions[which(labels==neg.label)] ))) )
  }


.performance.root.mean.squared.error <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      ## convert labels from factor to numeric values
      labels <- as.numeric(levels(labels))[labels]
      if (any(is.na(labels))) {
          stop("For rmse predictions have to be numeric.")
      }
      list( c(),  sqrt( 1/length(predictions) *
                      sum( (predictions - labels)^2 ))  )
  }

## ----------------------------------------------------------------------------
## Derived measures:
## ----------------------------------------------------------------------------

.performance.sar <- function( predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {

    pred <- prediction( predictions, labels)
    perf.acc <- performance( pred, measure="acc")
    perf.rmse <- performance( pred, measure="rmse")
    perf.auc <- performance( pred, measure="auc")

    list(cutoffs,
         1/3 * (perf.acc@y.values[[1]] +
                (1 - perf.rmse@y.values[[1]]) +
                perf.auc@y.values[[1]]))
}

## ----------------------------------------------------------------------------
## Measures taking into account actual cost considerations
## ----------------------------------------------------------------------------

.performance.expected.cost <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred) {
      
      ## kick out suboptimal values (i.e. fpr/tpr pair for which another one
      ## with same fpr and higher tpr exists, 
      ## or one for which one with same tpr but lower fpr exists

      if (n.neg==0 || n.pos==0) {
          stop(paste("At least one positive and one negative sample are",
                     "needed to compute a cost curve."))
      }
      fpr <- fp / n.neg
      tpr <- tp / n.pos

      ## sort by fpr (ascending), in case of ties by descending tpr
      ind <- order(fpr,-tpr)
      
      fpr <- fpr[ind]
      tpr <- tpr[ind]
      ## for tied fprs, only the one with the highest tpr is kept
      ind <- !duplicated(fpr)
      fpr <- fpr[ind]
      tpr <- tpr[ind]

      ## for tied tprs, only keep the one with the lowest fpr
      ind <- order(-tpr,fpr)
      fpr <- fpr[ind]
      tpr <- tpr[ind]
      ind <- !duplicated(tpr)
      fpr <- fpr[ind]
      tpr <- tpr[ind]

      if (!any(0==fpr & 0==tpr)) {
          fpr <- c(0,fpr)
          tpr <- c(0,tpr)
      }
      if (!any(1==fpr & 1==tpr)) {
          fpr <- c(fpr,1)
          tpr <- c(tpr,1)
      }
      
      ## compute all functions
      f <- list()
      for (i in 1:length(fpr)) {
          f <- c(f, .construct.linefunct( 0, fpr[i], 1, 1-tpr[i] ))
      }
      
      ## compute all intersection points
      x.values <- c()
      y.values <- c()
      for (i in 1:(length(fpr)-1)) {
          for (j in (i+1):length(fpr)) {
              ans <- .intersection.point( f[[i]], f[[j]] )
              if (all(is.finite(ans))) {
                  y.values.at.current.x <- c()
                  for (k in 1:length(f)) {
                      y.values.at.current.x <- c(y.values.at.current.x,
                                                 f[[k]](ans[1]))
                  }
                  if (abs(ans[2] - min(y.values.at.current.x )) <
                      sqrt(.Machine$double.eps)) {
                      
                      x.values <- c(x.values, ans[1])
                      y.values <- c(y.values, ans[2])
                  }
              }
          }
      }

      if (!any(0==x.values & 0==y.values)) {
          x.values <- c(0,x.values)
          y.values <- c(0,y.values)
      }
      if (!any(1==x.values & 0==y.values)) {
          x.values <- c(x.values,1)
          y.values <- c(y.values,0)
      }

      ind <- order( x.values)
      list( x.values[ind], y.values[ind] )
  }


.performance.cost <-
  function(predictions, labels, cutoffs, fp, tp, fn, tn,
           n.pos, n.neg, n.pos.pred, n.neg.pred, cost.fp, cost.fn) {
      
    n.samples <- n.pos + n.neg
    cost <- ((n.pos / n.samples) * (fn / n.pos) * cost.fn +
             (n.neg / n.samples) * (fp / n.neg) * cost.fp)
    list( cutoffs, cost )
}



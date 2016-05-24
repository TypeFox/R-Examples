pamr.options <- list(debug=TRUE, #whether to turn on debugging or not
                     err.file=ifelse(.Platform$OS.type=="windows", "C:/pamrtrace.txt", "pamrtrace.txt"),
                     image.file=ifelse(.Platform$OS.type=="windows", "C:/pamrimage.Rdata", "pamrimage.Rdata"),
                     reserved.class.label="Unspecified")

##
## Our error handler
##
pamr.xl.error.trace <- function() {
  err.message <- geterrmessage()
  if (!is.null(pamr.options$image.file)) {
    save.image(pamr.options$image.file)
  }
  if (!is.null(pamr.options$err.file)) {
    sink(pamr.options$err.file)
    print(err.message)
    traceback()
    sink()
  }
  winDialog(type="ok", message=err.message)
}

##
## Upon loading, if we are in a windows environment, we use the windows
## dialog mechanism to display errors. Useful for debugging COM apps
##
.onLoad <- function(lib, pkg) {

# Rob changed this next line on  apr 10, 2005, requested by Uwe Ligges

#  if ( .Platform$OS.type == "windows"  ) {
  if ( .Platform$OS.type == "windows" && interactive() ) {
    options(error=pamr.xl.error.trace)
  }

}

##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath){
  if ( .Platform$OS.type == "windows") {
    options(error=NULL)
  }
}


pamr.xl.get.threshold.range  <- function(fit) {
  return(range(fit$threshold))
}

pamr.xl.get.soft.class.labels  <- function(fit, survival.times, censoring.status) {
  proby <-  pamr.surv.to.class2(survival.times, censoring.status,
                                n.class=fit$ngroup.survival)$prob
  group <-apply(proby,1,which.is.max)
  return(group)
}


pamr.xl.compute.offset <- function(data, offset.percent=50, prior=prior){
  x <- data$x
  y <- data$y
  n.class <- table(y)
  if(min(n.class)==1){stop("Error: each class must have >1 sample")}
  norm.cent <-NULL
  n <- sum(n.class)
  xtest <- x
  ntest <- ncol(xtest)
  K <- length(prior)
  p <- nrow(x)
  Y <- model.matrix( ~ factor(y) - 1, data = list(y = y))
  dimnames(Y) <- list(NULL, names(n.class))
  centroids <- scale(x %*% Y, FALSE, n.class)
  sd <- rep(1, p)
  xdif <- x - centroids %*% t(Y)
  sd <- (xdif^2) %*% rep(1/(n - K), n)
  sd <- drop(sqrt(sd))
  offset  <- quantile(sd, offset.percent/100)
  return(offset)
}

pamr.xl.get.offset  <- function() {
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    pamr.xl.data <- get("pamr.xl.data", envir=.GlobalEnv)
    x.train <- get("x.train", envir=.GlobalEnv)
    if (exists("x.train")) {
        temp=x.train$offset
        if(is.null(temp)){temp=0}
        return (temp)
    } else {
        return (pamr.xl.compute.offset(pamr.xl.data,
                                       offset.percent=pamr.xl.training.parameters$offset.percent,
                                       prior=pamr.xl.training.parameters$prior))
  }
}

pamr.xl.derive.adjusted.prior  <- function(prior, data) {
  ## Check this next code in if statement. For survival setting, it is always uniform
  ## and so the check may not be needed. Anyway, needs cleaning....
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    if (pamr.xl.survival.setting) {
    s  <-  pamr.xl.get.uniform.prior(data, nclasses=pamr.xl.training.parameters$ngroup.survival)
    return (list(prior=s, prior.name="Uniform Prior"))
  } else {
    s  <- pamr.xl.get.sample.prior(data)
    temp <- prior - s
    if (sum(temp*temp) < pamr.xl.training.parameters$epsilon) {
      return (list (prior=s, prior.name="Sample Prior"))
    } else {
      s  <-  pamr.xl.get.uniform.prior(data)
      temp  <- prior - s
      if (sum(temp*temp) < pamr.xl.training.parameters$epsilon) {
        return (list (prior=s, prior.name="Uniform Prior"))
      } else {
        return (list (prior=prior, prior.name="Custom Prior"))
      }
    }
  }
}

#pamr.xl.get.default.training.parameters <- function(data) {
#  if (pamr.xl.survival.setting) {
#    return (list(offset.percent=50,
#                 prior=pamr.xl.get.uniform.prior(data, nclasses=2),
#                 prior.name="Uniform Prior",
#                 sign.contrast="both",
#                 epsilon=1e-7,
#                 ngroup.survival=2,
#                 survival.method="Kaplan Meier"))
#  } else {
#    return (list(offset.percent=50,
#                 prior=pamr.xl.get.sample.prior(data),
#                 prior.name="Sample Prior",
#                 sign.contrast="both",
#                 epsilon=1e-7,
#                 ngroup.survival=2,
#                 survival.method="Kaplan Meier"))
#  }
#}



pamr.xl.get.default.training.parameters <- function(data) {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.regression.setting <- get("pamr.xl.regression.setting", envir=.GlobalEnv)

    if (pamr.xl.survival.setting) {
    return (list(offset.proportion=0.5,
                 offset.percent=50,
                 prior=NULL,
                 prior.name=NULL,
                 sign.contrast="both",
                 epsilon=1e-7,
                 ngroup.survival=2,
                 decorrelate=FALSE,
                 n.components=1))
                }
 if (pamr.xl.regression.setting) {

    return (list(offset.proportion=0.5,
                offset.percent=50,
                 prior=NULL,
                 prior.name=NULL,
                 sign.contrast="both",
                 epsilon=1e-7,
                 ngroup.survival=NULL,
                 decorrelate=FALSE,
                 n.components=1))
 }
if(!pamr.xl.survival.setting & !pamr.xl.regression.setting){
    return (list(offset.percent=50,
                 prior=pamr.xl.get.sample.prior(data),
                 prior.name="Sample Prior",
                 sign.contrast="both",
                 epsilon=1e-7,
                 ngroup.survival=NULL,
                 n.components=NULL))

  }
}

## Return the uniform prior on class labels
pamr.xl.get.uniform.prior  <- function(data, nclasses=NULL) {
  if (is.null(nclasses)) {
    w <- table(data$y)
    n  <- length(w)
  } else {
    n = nclasses
  }
  return(rep(1.0/n, n))
}

## Return the sample proportion prior on class labels
pamr.xl.get.sample.prior  <- function(data) {
  w <- table(data$y)
  return(w/sum(w))
}

pamr.xl.get.class.names  <- function() {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    pamr.xl.data <- get("pamr.xl.data", envir=.GlobalEnv)

    if (pamr.xl.survival.setting) {
    return(as.character(1:pamr.xl.training.parameters$ngroup.survival))
  } else {
    return(names(table(pamr.xl.data$y)))
  }
}



#pamr.xl.get.class.labels  <- function() {
#  if (pamr.xl.survival.setting) {
#    return(rep(" ", length(pamr.xl.survival.times)))
#  } else {
#    return(pamr.xl.data$y)
#  }
#}

pamr.xl.get.class.labels  <- function() {
    pamr.xl.data <- get("pamr.xl.data", envir=.GlobalEnv)
    return(pamr.xl.data$y)
}


pamr.xl.get.number.of.classes  <- function() {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    pamr.xl.data <- get("pamr.xl.data", envir=.GlobalEnv)

    if (pamr.xl.survival.setting) {
        return(pamr.xl.training.parameters$ngroup.survival)
    } else {
        return(length(names(table(pamr.xl.data$y))))
    }
}

#pamr.xl.process.data <- function(use.old.version=FALSE) {
#
#  res <- list(x=pamr.xl.raw.data, y=pamr.xl.class.labels, genenames=pamr.xl.gene.names,
#              geneid=pamr.xl.gene.ids, samplelabels=pamr.xl.sample.labels,
#              batchlabels=pamr.xl.batch.labels, survival.time=pamr.xl.survival.times,
#              censoring.status=pamr.xl.censoring.status)
#
#  if (pamr.xl.data.has.missing.values) {
#    if (use.old.version) {
#      res <- pamr.knnimpute.old(res, k = pamr.xl.knn.neighbors)
#    } else {
#      res <- pamr.knnimpute(res, k = pamr.xl.knn.neighbors)
#    }
#  }
#  return(res)
#}

pamr.xl.process.data <- function(use.old.version=FALSE) {

# in this new version, the outcome is always stored in y
# the survival times component  is no longer used. Superpc now handles
# both the surival and regression problems

    pamr.xl.class.labels <- get("pamr.xl.class.labels", envir=.GlobalEnv)
    pamr.xl.survival.times <- get("pamr.xl.survival.times" , envir=.GlobalEnv)
    pamr.xl.raw.data <- get("pamr.xl.raw.data", envir=.GlobalEnv)
    pamr.xl.gene.names <- get("pamr.xl.gene.names", envir=.GlobalEnv)
    pamr.xl.gene.ids <- get("pamr.xl.gene.ids", envir=.GlobalEnv)
    pamr.xl.sample.labels <- get("pamr.xl.sample.labels", envir=.GlobalEnv)
    pamr.xl.batch.labels <- get("pamr.xl.batch.labels", envir=.GlobalEnv)
    pamr.xl.censoring.status <- get("pamr.xl.censoring.status", envir=.GlobalEnv)
    pamr.xl.data.has.missing.values <- get("pamr.xl.data.has.missing.values", envir=.GlobalEnv)
    pamr.xl.knn.neighbors <- get("pamr.xl.knn.neighbors", envir=.GlobalEnv)

 if(!is.null(pamr.xl.class.labels)){
    y=pamr.xl.class.labels
   }
  if(is.null(pamr.xl.class.labels)){
    y=pamr.xl.survival.times
   }




 res <- list(x=pamr.xl.raw.data, y=y, genenames=pamr.xl.gene.names,
              geneid=pamr.xl.gene.ids, samplelabels=pamr.xl.sample.labels,
              batchlabels=pamr.xl.batch.labels,
              censoring.status=pamr.xl.censoring.status)

  if (pamr.xl.data.has.missing.values) {
    if (use.old.version) {
      res <- pamr.knnimpute.old(res, k = pamr.xl.knn.neighbors)
    } else {
      res <- pamr.knnimpute(res, k = pamr.xl.knn.neighbors)
    }
  }
  return(res)
}


pamr.xl.compute.cv.confusion  <- function (fit, cv.results, threshold) {
  threshold.rank  <- which(rank(abs(cv.results$threshold - threshold))==1)
  t.threshold  <- cv.results$threshold[threshold.rank]
  true  <- cv.results$y
  predicted  <- cv.results$yhat[, threshold.rank]
  tt <- table(true, predicted)
  tt1 <- tt
   diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  overall.err  <- round(sum(tt1)/sum(tt), 3)
  return(list(confusion.matrix=tt, overall.error=overall.err, threshold=round(t.threshold, 5)))
 }
pamr.xl.compute.confusion  <- function (fit, threshold) {
  ii <- (1:length(fit$threshold))[fit$threshold >= threshold]
  ii <- ii[1]
  predicted <- fit$yhat[, ii]
  if(!is.null(fit$y)){
    true <- fit$y[fit$sample.subset]
    tt <- table(true, predicted)
  } else {
    true <- fit$proby[fit$sample.subset,]
    ytemp<- apply(true,1,which.is.max)
    temp <- c(predicted,names(table(ytemp)))
    nams <- names(table(temp))
    Yhat <- model.matrix( ~ factor(temp) - 1,
                         data = list(y = temp))
    Yhat <- Yhat[1:length(predicted),]
    tt <- matrix(NA,nrow=length(fit$prior),ncol=length(fit$prior))
    for(i in 1:length(fit$prior)){
      for(j in 1:length(fit$prior)){
        tt[i,j] <- sum(true[,i]*Yhat[,j])
      }
    }
    dimnames(tt) <- list(names(table(ytemp)),nams)
  }
  tt1 <- tt
  diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  overall.err  <- round(sum(tt1)/sum(tt), 3)
  return(list(confusion.matrix=tt, overall.error=overall.err))
}

pamr.xl.is.a.subset  <- function(a, y) {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)

    if (pamr.xl.survival.setting) {
    x  <- as.character(1:pamr.xl.training.parameters$ngroup.survival)
  } else {
    x  <- a$y
  }
  if (nlevels(factor(x)) == nlevels(factor(c(x, y[!is.na(y)])))) {
    return (1)  # True
  } else {
    return (0)  # False
  }
}

pamr.xl.listgenes.compute  <- function (fit, data, threshold, fitcv=NULL,  genenames = FALSE) {
  x <- data$x[fit$gene.subset, fit$sample.subset]
  if (genenames) {
    gnames <- data$genenames[fit$gene.subset]
  }
  if (!genenames) {
    gnames <- NULL
  }
  geneid <- data$geneid[fit$gene.subset]
  if(!is.null(fit$y)){
       nc <- length(fit$y)
      }
 if(is.null(fit$y)){
       nc <- ncol(fit$proby)
      }
 clabs <- colnames(fit$centroids)

  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
  d <- (cen - fit$centroid.overall)[aa, ]/fit$sd[aa]

  gene.order <- order(-apply(abs(d), 1, max))
  d <- round(d, 4)
  g <- gnames[aa]
  g1 <- geneid[aa]
  if (is.null(gnames)) {
    gnhdr <- NULL
  }
  if (!is.null(gnames)) {
    gnhdr <- "name"
  }

if(!is.null(fitcv)){
nfold=length(fitcv$cv.objects)

ind=matrix(F,nrow=nrow(x),ncol=nfold)
ranks=NULL
for( ii in 1:nfold){
        cen=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=0, type="centroid")
         dtemp <- (cen - fitcv$cv.objects[[ii]]$centroid.overall)[,, drop=FALSE]/fitcv$cv.objects[[ii]]$sd
          r <- apply(abs(dtemp), 1, max)
        ranks=cbind(ranks,rank(-abs(r)))

        junk=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=threshold, type="nonzero")
        ind[junk,ii]=T
}

av.rank=apply(ranks,1,mean)
av.rank=round(av.rank[aa],2)
prop=apply(ind[aa,,drop=F],1,sum)/nfold
}

  options(width = 500)
  schdr <- paste(clabs, "score", sep = " ")


if(is.null(fitcv)){
res <- cbind(as.character(g1), g, d)[gene.order, ]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr))

}
if(!is.null(fitcv)){
  res <- cbind(as.character(g1), g, d, av.rank, prop)[gene.order, ]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr, "av-rank-in-CV", "prop-selected-in-CV"))
}


  return(list(gene.headings = dimnames(res)[[2]],
              gene.ids = res[ , 2],   # This was switched with gene.names.
              gene.names = res[ , 1],
              gene.scores = res[ , -(1:2)]))
  ##print(res, quote = FALSE)
}
pamr.xl.plot.test.probs.compute  <- function(fit, new.x, newx.classes, missing.class.label,
	threshold, sample.labels=NULL) {

    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    pamr.xl.test.survival.times <- get("pamr.xl.test.survival.times", envir=.GlobalEnv)
    pamr.xl.test.censoring.status <- get("pamr.xl.test.censoring.status", envir=.GlobalEnv)
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)

    predicted.probs  <- pamr.xl.predict.test.probs(fit, new.x, threshold=threshold)
  py  <- pamr.xl.predict.test.class.only(fit, new.x, threshold=threshold)

  if (pamr.xl.survival.setting & !is.null(pamr.xl.test.survival.times) ) {
    proby <-  pamr.surv.to.class2(pamr.xl.test.survival.times, pamr.xl.test.censoring.status,
                                  n.class=fit$ngroup.survival)$prob
    group <-apply(proby,1,which.is.max)
    order.classes  <- order(group)
    actual.classes  <- group[order.classes]
  } else {
    order.classes  <- order(newx.classes)
    actual.classes <- newx.classes[order.classes]
    actual.classes[is.na(actual.classes)] <- missing.class.label
  }
  pp  <- predicted.probs[, order.classes]
  ny  <- py$predicted[order.classes]
  n  <- length(ny)
  ss  <- sample.labels
  if (!is.null(ss)) {
    ss  <- ss[order.classes]
  }
  if (pamr.xl.survival.setting) {
    training.classes <- levels(factor(as.character(1:pamr.xl.training.parameters$ngroup.survival)))
  } else {
    training.classes  <- levels(factor(fit$y))
  }

  return (list(x = 1:n,
               y = t(pp),
               x.label = "Sample",
               y.label = "Predicted Test Probabilities",
               y.names = training.classes,
               y.lines = cumsum(table(actual.classes)) + 0.5,
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               panel.names = levels(factor(actual.classes)),
               x.names = ss))
}



pamr.xl.plot.training.error.compute  <- function(trained.object) {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)

    if (pamr.xl.survival.setting) {
    n  <- length(trained.object$survival.time)
  } else {
    n  <- length(trained.object$y)
  }
  return (list(x = trained.object$threshold,
               y = trained.object$errors/n,
               y.ytop = trained.object$nonzero,
               x.label = "Threshold",
               y.label = "Training Error"))
}
pamr.xl.plotcen.compute  <- function(fit, data, threshold) {
  genenames <- data$genenames[fit$gene.subset]
  x <- data$x[fit$gene.subset, fit$sample.subset]
  clabs <- colnames(fit$centroids)
  scen <- pamr.predict(fit, data$x, threshold = threshold, type = "cent")
  dif <- (scen - fit$centroid.overall)/fit$sd
  if(!is.null(fit$y)){
       nc <- length(unique(fit$y))
  }
   if(is.null(fit$y)){
      nc <- ncol(fit$proby)
}
  o <- drop(abs(dif) %*% rep(1, nc)) > 0
  d <- dif[o,  ]
  nd <- sum(o)
  genenames <- genenames[o]
  xx <- x[o,  ]
  oo <- order(apply(abs(d), 1, max))
  d <- d[oo,  ]
  genenames <- genenames[oo]
  ##win.metafile()
  par(mar = c(1, 5, 1, 1), col = 1)
  plot(rep(2, nd) + d[, 1], 1:nd, xlim = c(0, 2*nc+1), ylim = c(1, nd + 3),
       type = "n", xlab = "", ylab = "", axes = FALSE)
  box()
  abline(h = seq(nd), lty = 3, col = 7)
  jj <- rep(0, nd)
  for(j in 1:nc) {
    segments(jj + 2 * j, seq(nd), jj + 2 * j + d[, j], seq(nd), col
             = j + 1, lwd = 4)
    lines(c(2 * j, 2 * j), c(1, nd), col = j + 1)
    text(2 * j, nd + 2, label = clabs[j], col = j + 1)
  }
  g <- substring(genenames, 1, 20)
  text(rep(0, nd), seq(nd), label = g, cex = 0.4, adj = 0, col = 1)
  if (.Platform$OS.type == "windows") {
      savePlot("", type="wmf")
  }
  dev.off()
#  pamr.plot.y <<- matrix(d, nrow=dim(d)[1])
#  pamr.plot.x <<- seq(nd)
#  pamr.plot.seriesnames <<- dimnames(d)[[2]]
#  pamr.plot.genenames <<- genenames

  return(TRUE)
}
pamr.xl.plotcv.compute  <- function(aa) {
    pamr.xl.survival.setting <- get("pamr.xl.survival.setting", envir=.GlobalEnv)
    n <- nrow(aa$yhat)
  y <- aa$y
  if(!is.null(aa$newy)) {
    y <- aa$newy[aa$sample.subset]
  }
  nc <- length(table(y))
  nfolds <- length(aa$folds)
  err <- matrix(NA, ncol = ncol(aa$yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(aa$yhat), nrow = n)
  ni <- rep(NA, nfolds)
  for(i in 1:nfolds) {
    ii <- aa$folds[[i]]
    ni[i] <- length(aa$folds[[i]])
    err[i,  ] <- apply(temp[ii,  ] != aa$yhat[ii,  ], 2, sum)/ni[i]
  }
  se <- sqrt(apply(err, 2, var)/nfolds)

  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(aa$threshold)-1)
  for(i in 1:(length(aa$threshold) - 1)) {
    s <- pamr.confusion(aa, aa$threshold[i], extra = FALSE)
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }
  if (pamr.xl.survival.setting) {
    p.values <- aa$pvalue.survival
  } else {
    p.values  <- NULL
  }

  return (list(x = aa$threshold,
               y = aa$error,
               x.label = "Threshold",
               y.label = "Misclassification Error",
               y.se = se,
               p.values = p.values,
               y.ytop = aa$size,
               cv.err = t(err2),
               cv.legend = dimnames(table(y))[[1]]))

}
pamr.xl.plotcvprob.compute  <- function(fit, data, threshold) {
  ii <- (1:length(fit$threshold))[fit$threshold > threshold]
  ii <- ii[1]
  ss <- data$samplelabels
  pp <- fit$prob[,  , ii]
  if(is.null(fit$newy)) {
    y <- fit$y[fit$sample.subset]
  }
  if(!is.null(fit$newy)) {
    y <- fit$newy[fit$sample.subset]
  }
  o <- order(y)
  y <- y[o]
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,  ]
  n <- nrow(ppp)
  nc <- length(unique(y))


#  axis(2, labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", ""))
#  if (!is.null(ss)) {
#    pamr.plot.x.names <<- ss
#  }

  return (list(x = 1:n,
               y = ppp,
               x.label = "Sample",
               y.label = "CV Probabilities",
               y.names = levels(y),
               y.lines = cumsum(table(fit$y)),
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               x.names = ss))

#   for(j in 1:nc) {
#     points(1:n, ppp[, j], col = j + 1)
#   }
#   for(j in 1:(nc - 1)) {
#     abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
#   }
#   h <- c(0, table(y))
#   for(j in 2:(nc + 1)) {
#     text(sum(h[1:(j - 1)]) + 0.5 * h[j], 1.02, label = levels(y)[j -
#                                                  1], col = j)
#   }
#   abline(h = 1)
#   if(!is.null(ss)) {
#     text(1:length(ss), 1.1, labels = ss, srt = 90, cex = 0.7)
#   }
  ##if(!is.null(ss)){axis(3,labels=ss,at=1:length(ss),srt=90)}
}
pamr.xl.predict.test.class<- function(fit, newx, threshold, test.class.labels) {
  predicted  <- pamr.predict(fit, newx, threshold, type="class")
  return(list(confusion.matrix=table(test.class.labels, predicted), predicted=as.vector(predicted)))
}

pamr.xl.predict.test.surv.class <- function(fit, newx, threshold, survival.times, censoring.status) {
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    predicted  <- pamr.predict(fit, newx, threshold, type="class")
  soft.probs  <- pamr.surv.to.class2(survival.times, censoring.status,
                                     n.class=pamr.xl.training.parameters$ngroup.survival)$prob
  w  <- pamr.test.errors.surv.compute(soft.probs, predicted)
  return(list(confusion.matrix=w$confusion, predicted=as.vector(predicted)))
}

pamr.xl.predict.test.class.only  <- function(fit, newx, threshold) {
  return(list(predicted=as.vector(pamr.predict(fit, newx, threshold, type="class"))))
}

pamr.xl.predict.test.probs  <- function(fit, newx, threshold) {
  predicted  <- pamr.predict(fit, newx, threshold, type="posterior")
  return(t(predicted))
}

pamr.xl.test.data.impute  <- function(x, k, use.old.version=FALSE) {
    pamr.xl.knn.neighbors <- get("pamr.xl.knn.neighbors", envir=.GlobalEnv)
    if (use.old.version) {
    res <- pamr.knnimpute.old(list(x=x), k = pamr.xl.knn.neighbors)
  } else {
    res  <- pamr.knnimpute(list(x=x), k = pamr.xl.knn.neighbors)
  }
  return(res$x)
}

pamr.xl.test.errors.surv.compute <- function(fit, newx, threshold=fit$threshold, survival.times, censoring.status) {
    pamr.xl.training.parameters <- get("pamr.xl.training.parameters", envir=.GlobalEnv)
    prediction.errs  <- vector(mode="numeric", length=length(threshold))
  soft.probs  <- pamr.surv.to.class2(survival.times, censoring.status,
                                     n.class=pamr.xl.training.parameters$ngroup.survival)$prob
  for (i in 1:length(threshold)) {
    predicted  <- pamr.predict(fit, newx, threshold=threshold[i], type="class")
    w  <- pamr.test.errors.surv.compute(soft.probs, predicted)
    prediction.errs[i]  <- w$error
  }
  return(list(x=threshold, y=prediction.errs, x.label="Threshold", y.label="Test Error"))
}



pamr.xl.test.errors.compute  <- function(fit, newx, newx.classes, threshold=fit$threshold,
                                         prior = fit$prior,  threshold.scale = fit$threshold.scale,
                                         ...) {
  n  <- length(which(!is.na(newx.classes)))
## Note: n is assumed to be nonzero! Check before calling!
  actual.classes  <- newx.classes
  prediction.errs  <- vector(mode="numeric", length=length(threshold))

  for(i in 1:length(threshold)){
    t <- pamr.predict(fit,newx,threshold=threshold[i],type="class",...)
    prediction.errs[i]  <- length(which(t != actual.classes)) / n
  }

  return(list(x=threshold, y=prediction.errs, x.label="Threshold", y.label="Test Error"))

}

pamr.xl.transform.class.labels  <- function(x) {
  y  <- x
  y[is.na(y)]  <- " "
  return(y)
}

pamr.xl.transform.data <- function(data) {

    pamr.xl.take.cube.root <- get("pamr.xl.take.cube.root", envir=.GlobalEnv)
    pamr.xl.batch.labels.present <- get("pamr.xl.batch.labels.present", envir=.GlobalEnv)
    pamr.xl.center.columns <- get("pamr.xl.center.columns", envir=.GlobalEnv)
    pamr.xl.scale.columns <- get("pamr.xl.scale.columns",  envir=.GlobalEnv)

    if (pamr.xl.take.cube.root) {
    data$x = pamr.cube.root(data$x)
  }

  if (pamr.xl.batch.labels.present) {
    data <- pamr.batchadjust(data)
  }

  if (pamr.xl.center.columns && pamr.xl.scale.columns) {
    data$x = scale(data$x, center=TRUE, scale=TRUE)
  } else if (pamr.xl.center.columns) {
    data$x = scale(data$x, center=TRUE, scale=FALSE)
  } else if (pamr.xl.scale.columns) {
    data$x = scale(data$x, center=FALSE, scale=TRUE)
  }

  return (data)
}

pamr.xl.transform.test.data <- function(test.x) {
  res <- test.x

    pamr.xl.take.cube.root <- get("pamr.xl.take.cube.root", envir=.GlobalEnv)
    pamr.xl.center.columns <- get("pamr.xl.center.columns", envir=.GlobalEnv)
    pamr.xl.scale.columns <- get("pamr.xl.scale.columns",  envir=.GlobalEnv)

  if (pamr.xl.take.cube.root) {
    res = pamr.cube.root(res)
  }

  if (pamr.xl.center.columns && pamr.xl.scale.columns) {
    res = scale(res, center=TRUE, scale=TRUE)
  } else if (pamr.xl.center.columns) {
    res = scale(res, center=TRUE, scale=FALSE)
  } else if (pamr.xl.scale.columns) {
    res = scale(res, center=FALSE, scale=TRUE)
  }

  return (res)
}

pamr.xl.plotsurvival<- function(fit, data, threshold) {
  group  <- pamr.predict(fit, data$x, threshold=threshold)
  ## plots Kaplan-Meier curves stratified by "group"
#  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(fit$survival.time, fit$censoring.status)~as.factor(group))
  ## win.metafile()
  plot(junk, col=2:(2+n.class-1) ,xlab= "Time", ylab="Probability of survival", main="Survival Plot")
  legend(.8*max(fit$survival.time),.9, col=2:(2+n.class-1), lty=rep(1,n.class),
         legend=as.character(1:n.class))
  if (.Platform$OS.type == "windows") {
      savePlot("", type="wmf")
  }
  dev.off()
  return(TRUE)
}

pamr.xl.plotsurvival.test <- function(fit, newx, survival.time, censoring.status, threshold) {
  group  <- pamr.predict(fit, newx, threshold=threshold)
  ## plots Kaplan-Meier curves stratified by "group"
  #require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status)~as.factor(group))
  ## win.metafile()
  plot(junk, col=2:(2+n.class-1) ,xlab= "Time", ylab="Probability of survival", main="Test Survival Plot")
  legend(.8*max(survival.time),.9, col=2:(2+n.class-1), lty=rep(1,n.class),
         legend=as.character(1:n.class))
  if (.Platform$OS.type == "windows") {
      savePlot("", type="wmf")
  }
  dev.off()
  return(TRUE)
}

pamr.xl.plotsurvival.strata <- function(fit, data) {
  group <-apply(fit$proby,1,which.is.max)
  #require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(data$survival.time, data$censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(data$survival.time, data$censoring.status) ~ as.factor(group))

  pv <- 1-pchisq(2*(junk2$loglik[2]-junk2$loglik[1]),df=n.class-1)

  if(!is.null(fit$cutoffs.survival)){
    labels <- rep(NULL,n.class)
    labels[1] <- paste("(1)   ","<= ", round(fit$cutoffs.survival[1],2),sep="")
    if(n.class>2){
      for(i in 2:(n.class-1)){
        labels[i] <- paste("(",as.character(i),")  ", " > ",
                           round(fit$cutoffs.survival[i-1],2), "  & <= ",
                           round(fit$cutoffs.survival[i],2), sep="")
      }}
    labels[n.class] <-  paste("(",as.character(n.class),")  ", " > ",round(fit$cutoffs.survival[n.class-1],2),sep="")
  }

  else{labels <- as.character(1:n.class)}

  ##win.metafile()
  plot(junk, col = 2:(2 + n.class - 1), xlab = "Time", ylab = "Probability of survival",
       main="Survival Strata Plot")
  legend(.01* max(fit$survival.time), 0.2, col = 2:(2 + n.class -
                                             1), lty = rep(1, n.class), legend = labels)
  text(0.1 * max(fit$survival.time), .25, paste("pvalue=",as.character(round(pv,4))))
  if (.Platform$OS.type == "windows") {
      savePlot("", type="wmf")
  }
  dev.off()
  return(TRUE)
}

pamr.xl.test.get.soft.classes  <- function(fit, survival.times, censoring.status) {
  proby <-  pamr.surv.to.class2(survival.times, censoring.status,
                                n.class=fit$ngroup.survival)$prob
  soft.classes  <- apply(proby,1,which.is.max)
  return(list(classes=soft.classes, probs = t(proby)))
}

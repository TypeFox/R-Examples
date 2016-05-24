####
## All functions that adjust a generic function for objects of classes PK, PKtest and halflife 
####

PKNews <- function() file.show(system.file("NEWS", package="PK"))


plot.halflife <- function (x, xlab="Time", ylab="Concentration", main="Half-life Estimation", xlim=NULL, ylim=NULL, add=FALSE, ...) {

    b1 <- x$parms[2,1]
    a1 <- x$parms[3,1]
    b2 <- x$parms[2,2]
    a2 <- x$parms[3,2]

    if (is.null(xlim)) {xlim <- c(min(x$time), max(x$time))}
    if (is.null(ylim)) {ylim <- c(min(x$conc), max(x$conc))}

    if(add==FALSE){plot(x=x$time, y=x$conc, xlim=c(xlim[1], xlim[2]), ylim=c(ylim[1], ylim[2]), xlab=xlab, ylab=ylab, main=main, ...)}
    if(add==TRUE){points(x=x$time, y=x$conc, ...)}

    switch(x$method, lee = {
        if (is.na(x$chgpt)) {x$chgpt <- min(x$time)}
        curve(10^(b1 * x + a1), from=xlim[1], to=x$chgpt, add=TRUE, ...)
        curve(10^(b2 * x + a2), from=x$chgpt, to=xlim[2], add=TRUE, ...)
      }, biexp = {
        if (b1 == b2) {curve(a1 * exp(-b1 * x), from=xlim[1], to=xlim[2], add=TRUE, ...)}
        if (b1 != b2) {
           curve(a1 * exp(-b1 * x) + a2 * exp(-b2 * x), from=xlim[1], to=xlim[2], add=TRUE, ...)
        }
    }, )
}


###
# Prints point estimate and se's for class PK objects
###
print.PK <- function (x, digits=max(3, getOption("digits") - 4), ...) {

  if (x$design == "complete") {
    cat("Estimation for a complete data design\n\n")
  }else {
    if (x$design == "batch") {
      cat("Estimation for a batch design\n\n")
    } else {
      cat("Estimation for a serial sampling design\n\n")
    }
  }

  for(m in unique(x$CIs[,'method'])){
  if(is.na(x$CIs[1,"stderr"])) {
    res <- data.frame(format(round(estimator(x, se = FALSE),digits=digits), digits=digits, 
                      scientific=FALSE, nsmall=min(digits,2)), se=rep(NA,nrow(x$CIs)),
             ci = paste("(",format(round(x$CIs[x$CIs[,'method']==m, "lower"],digits=digits), digits = digits,
                  scientific=FALSE, nsmall=min(digits,2)), ";",
             format(round(x$CIs[x$CIs[,'method']==m, "upper"],digits=digits), digits = digits, scientific=FALSE, 
             nsmall=min(digits,2)), ")", sep = ""))
  }else{
    res <- data.frame(format(round(estimator(x, se = TRUE),digits=digits), digits=digits, 
                      scientific=FALSE, nsmall=min(digits,2)), 
             ci = paste("(",format(round(x$CIs[x$CIs[,'method']==m, "lower"],digits=digits), digits = digits,
                  scientific=FALSE, nsmall=min(digits,2)), ";",
             format(round(x$CIs[x$CIs[,'method']==m, "upper"],digits=digits), digits = digits, scientific=FALSE, 
             nsmall=min(digits,2)), ")", sep = ""))

  }
    colnames(res) <- c("Estimate", "SE", paste(x$conf.level * 100,  "% ", m, "-CI", sep = ""))
    print(res); cat("\n")
  }   
}


###
# Gives point estimates and confidence intervals for class PK objects
###
summary.PK<-function(object,...){

  if(object$design=="complete"){
    cat("Confidence intervals for a complete data design\n\n")
  }else{
    if(object$design=="batch"){
      cat("Confidence intervals for a batch design\n\n")
    }else{
      cat("Confidence intervals for a serial sampling design\n\n")
    }
  }

  if(dim(estimator(object))[1]==1){ # only one estimator
    cat(" Point estimate\n\n")
  }else{
    cat(" Point estimates\n\n")
  }
  print(estimator(object,se=TRUE))

  if(dim(estimator(object))[1]==1 & dim(ci(object))[1]==1){ # only one estimator and one ci
    cat("\n Confidence Interval\n\n")
  }else{
    cat("\n Confidence Intervals\n\n")
  }

  print(ci(object))
  if(any(object$CIs[,'method']=='t') ){
    cat("\n")
    cat("degrees of freedom: ")
    cat(object$CIs[object$CIs[,'method']=='t','df'][1])
    cat("\n")
  }
  if(any(object$CIs[,'method']=='fieller') ){
    cat("\n")
    cat("degrees of freedom: ")
    cat(object$CIs[object$CIs[,'method']=='fieller','df'][1])
    cat("\n")
  }
}

plot.PK <- function (x, bygroup=FALSE, col=NULL, pch=NULL, main=NULL, xlab="Time", ylab="Concentration", ylim=NULL, xlim=NULL, add=FALSE, ...) {

  tmax <- ifelse(is.null(xlim), max(unlist(x$time)), xlim[2])
  tmin <- ifelse(is.null(xlim), min(unlist(x$time)), xlim[1])
  cmax <- ifelse(is.null(ylim), max(unlist(x$conc)), ylim[2])
  cmin <- ifelse(is.null(ylim), min(unlist(x$conc)), ylim[1])

  if(bygroup & !nlevels(as.factor(unlist(x$group)))==2){
    warning('No grouping variable available.')
    bygroup <- FALSE
  }

  if(bygroup & add){
    stop('Adding plots by group is currently not implemented.')
  }

  if(!bygroup){
    if (x$design == "batch") {
      B <- length(x$time)
      if (is.null(main)) {main = "Concentration versus time plot (Batch Design)"}
      if (is.null(col)) {col <- 1:B}
      if (is.null(pch)) {pch <- 1:B}
      pch <- rep(pch, length = B)
      col <- rep(col, length = B)
      if(!add){
        plot(x = x$time[[1]], y = x$conc[[1]], xlim = c(tmin, tmax), ylim = c(cmin, cmax), xlab = xlab, ylab = ylab, main = main, pch = pch[1], col = col[1], ...)
        for (i in 2:B) {
          points(x = x$time[[i]], y = x$conc[[i]], pch = pch[i], col = col[i])
        }
        labs <- 1:B
        legend(x = tmax, y = cmax, xjust=1, legend = paste("Batch",  labs), pch = pch, col = col)
      }
      if(add){
        for (i in 1:B) {
          points(x$time[[i]], x$conc[[i]], pch = pch[i], col = col[i], ...)
        }
      }
    } else {
      if (is.null(col)) {col <- 1}
      if (is.null(pch)) {pch <- 1}
      if (x$design == "ssd") {
        if (is.null(main)) {main = "Concentration versus time plot (Serial Sampling Design)"}
        if(!add){plot(x = x$time, y = x$conc, xlim = c(tmin, tmax), ylim = c(cmin, cmax), main = main, xlab = xlab, ylab = ylab, pch = pch, col = col, ...)}
        if(add){points(x = x$time, y = x$conc, col = col, pch = pch, ...)}
      } else {
        if (is.null(main)) {main = "Concentration versus time plot (Complete Data Design)"}
        if(!add){plot(x = x$time[[1]], y = x$conc[[1]], xlim = c(tmin, tmax), ylim = c(cmin, cmax), main = main, xlab = xlab, ylab = ylab, pch = pch, col = col, ...)}
        if(add){points(x = x$time[[1]], y = x$conc[[1]], col = col, pch = pch, ...)}
      }
    }

  }else{ ## graphs separated by group
    g1 <- sort(unique(unlist(x$group)))[1]
    g2 <- sort(unique(unlist(x$group)))[2]
    if (x$design == "batch") {
      B <- length(x$time)
      if (is.null(main)) {main = "Concentration versus time plot by group (Batch Design)"}
      if (is.null(col)) {col <- 1}
      if (is.null(pch)) {pch <- 1:B}
      pch <- rep(pch, length = B)
      col <- rep(col, length = 1)
      if (is.numeric(col)) {
        col1 <- col; col2 <-col+2
      }else{
        col1 <- col; col2 <-which(palette()==col)+1
        if(is.na(col2)){col2<-1}
      } 
      if(!add){
        plot(x = x$time[[1]][x$group[[1]]==g1], y = x$conc[[1]][x$group[[1]]==g1], xlim = c(tmin, tmax), ylim = c(cmin, cmax), xlab = xlab, ylab = ylab, main = main, pch = pch[1], col = col1[1], ...)
        for (i in 2:B) {
          points(x = x$time[[i]][x$group[[i]]==g1], y = x$conc[[i]][x$group[[i]]==g1], pch = pch[i], col = col1)
        }
        for (i in 1:B) {
          points(x$time[[i]][x$group[[i]]==g2], x$conc[[i]][x$group[[i]]==g2], pch = pch[i], col = col2, ...)
        }
        labs <- 1:B
        legend(x = tmax, y = cmax, xjust=1, legend = paste("Batch",  labs), pch = pch, col = 1)
        legend(x = tmax*4/5, y = cmax, xjust=1, legend = c("Group 1", "Group 2"), pch = 1, col = c(col1,col2))
      }
      if(add){
        for (i in 1:B) {
          points(x$time[[i]], x$conc[[i]], pch = pch[i], col = col[i], ...)
        }
      }
    } else {
      if (is.null(col)) {col <- 1}
      if (is.numeric(col)) {
        col1 <- col; col2 <-col+1
      }else{
        col1 <- col; col2 <-which(palette()==col)+1
        if(is.na(col2)){col2<-1}
      }
      if (is.null(pch)) {pch <- 1}
      if (x$design == "ssd") {
        if (is.null(main)) {main = "Concentration versus time plot by group (Serial Sampling Design)"}
        if(!add){
          plot(x = x$time[x$group==g1], y = x$conc[x$group==g1], xlim = c(tmin, tmax), ylim = c(cmin, cmax), main = main, xlab = xlab, ylab = ylab, pch = pch, col = col1, ...)
          points(x = x$time[x$group==g2], y = x$conc[x$group==g2], col = col2, pch = pch+1, ...)
          legend(x = tmax, y = cmax, xjust=1, legend = c("Group 1", "Group 2"), pch = c(pch,pch+1), col = c(col1,col2))
        }
        if(add){points(x = x$time, y = x$conc, col = col, pch = pch, ...)}
      } else {
        if (is.null(main)) {main = "Concentration versus time plot by group (Complete Data Design)"}
        if(!add){
          plot(x = x$time[[1]][x$group[[1]]==g1], y = x$conc[[1]][x$group[[1]]==g1], xlim = c(tmin, tmax), ylim = c(cmin, cmax), main = main, xlab = xlab, ylab = ylab, pch = pch, col = col1, ...)
          points(x = x$time[[1]][x$group[[1]]==g2], y = x$conc[[1]][x$group[[1]]==g2], col = col2, pch = pch + 1, ...)
          legend(x = tmax, y = cmax, xjust=1, legend = c("Group 1", "Group 2"), pch = c(pch,pch+1), col = c(col1,col2))
        }
        if(add){points(x = x$time[[1]], y = x$conc[[1]], col = col, pch = pch, ...)}
      }
    }
  }
}

print.PKtest<-function(x,hyp=FALSE,...){

  if(dim(x$stat)[1]>1){
    title <- ' Hypothesis tests for a'
  }else{
    title <- ' Hypothesis test for a'
  }
  if(x$design=="complete"){
    cat(paste(title,' complete data design\n\n',sep=''))
  }else{
    if(x$design=="batch"){
      cat(paste(title,' batch design\n\n',sep=''))
    }else{
      cat(paste(title,' serial sampling design\n\n',sep=''))
    }
  }

  if(hyp){
    cat("Hypothesis:\n\n")
    if(x$alternative=='less') {
      hyps<-paste('     ',c('H0','H1'),': ',rep(rownames(x$stat),each=2),' ',c('>=','< '),rep(x$theta,each=2),c('\n','\n\n'))
    }else{
      if(x$alternative=='greater') {
        hyps<-paste('     ',c('H0','H1'),': ',rep(rownames(x$stat),each=2),' ',c('<=','> '),rep(x$theta,each=2),c('\n','\n\n'))
      }else{
        hyps<-paste('     ',c('H0','H1'),': ',rep(rownames(x$stat),each=2),' ',c('==','<>'),rep(x$theta,each=2),c('\n','\n\n'))
      }
    }
  }else{
    if(x$alternative=='less') {
      hyp0<-paste('     H0: Parameter >= theta\n')
      hyp1<-paste('     H1: Parameter <  theta\n\n')
    }else{
      if(x$alternative=='greater') {
        hyp0<-paste('     H0: Parameter <= theta\n')
        hyp1<-paste('     H1: Parameter >  theta\n\n')
      }else{
        hyp0<-paste('     H0: Parameter == theta\n')
        hyp1<-paste('     H1: Parameter <> theta\n\n')
      }
    }
    hyps <- rbind(hyp0,hyp1)
  }

  cat(hyps)

  reject <- ifelse(x$p.value<=1-x$conf.level,'*','')
  pval<-matrix(NA,dim(x$p.value))
  pval[x$p.value<0.0001]<-'<0.0001'
  pval[x$p.value>=0.0001]<-format(x$p.value[x$p.value>=0.0001],digits=1,nsmall=4)

  if(x$method=='t'){
    cat("t-test\n\n")  
    out <- data.frame(Statistic=round(x$stat,3),df=round(x$df,2),pval=pval,re=reject)
    colnames(out) <- c('Test statistic','df','p-value',' ')
  }else{
    if(x$method=='z'){
      cat("z-test\n\n")
    }else{
      if(x$method=='fieller'){
        cat("Test based on a Fieller interval\n\n")
      }else{
        cat("Resampling-test\n\n")
      }
    }
    out <- cbind(x$stat,pval,reject)
    out <- data.frame(Statistic=round(x$stat,3),pval=pval,re=reject)
    colnames(out) <- c('Test statistic','p-value',' ')
  }

  print(out)

  cat('\n\n* indicates a parameter less than ',1-x$conf.level,' and so one can reject\nthe corresponding null hypothesis at a significance level of ',1-x$conf.level,'.\n\n',sep='')
 
}

summary.PKtest <- function(object,...){
  print(object,hyp=TRUE)
}

########################################################################
## This program is Open Source Software: you can redistribute it      ##
## and/or modify it under the terms of the GNU General Public License ##
## as published by the Free Software Foundation, either version 3 of  ##
## the License, or (at your option) any later version.                ##
##                                                                    ##
## This program is distributed in the hope that it will be useful,    ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of     ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  ##
## General Public License for more details.                           ##
##                                                                    ##
## You should have received a copy of the GNU General Public License  ##
## along with this program. If not, see http://www.gnu.org/licenses/. ##
########################################################################
myround <- function(x, what=what, digits=digits){
  if(what=="sample size" & digits==0)
    return(ceiling(x))
  round(x,digits)
}


integer2numeric <- function(x){
  if(class(x)=="integer")
    return(as.numeric(x))
  x
}

rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol")) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

gsbCgrid <- function (Ngrd, bnds){
    ## Ngrd - number of grid points (the actual size of the
    ##        grid is the smallest glp that has more than Ngrd points)
    ## bnds - matrix containing bounds
    ## dim - dimensional
    glp <- c(3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377,
            610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657,
            46368, 75025)
    if (Ngrd > 75025) {
      N <- glp[22]
      k <- 1:N
      mat <- cbind((k - 0.5)/N, ((glp[ind - 1] * k - 0.5)/N)%%1)
      mat2 <- cbind(runif(Ngrd - N), runif(Ngrd - N))
      mat <- rbind(mat, mat2)
    } else if (Ngrd < 5) {
      i <- 1:Ngrd
      mat <- cbind((i * sqrt(2))%%1, (i * sqrt(3))%%1)
    } else {
      ind <- min((1:22)[glp >= Ngrd])
      N <- glp[ind]
      k <- 1:N
      mat <- cbind((k - 0.5)/N, ((glp[ind - 1] * k - 0.5)/N)%%1)
    }
    mat[, 1] <- mat[, 1] * (bnds[1, 2] - bnds[1, 1]) + bnds[1, 1]
    mat[, 2] <- mat[, 2] * (bnds[2, 2] - bnds[2, 1]) + bnds[2, 1]
    mat
}


#as methods for classes: "OCsimulation", "OCtrial", "OCprior", "OCcriteria" <<<-----------------

as.data.frame.gsbCriteria <- function(x,...){
    if (class(x)[2]=="array"){
    value <- c(as.vector(x[1,1,,]),as.vector(x[1,2,,]))
    probability <- c(as.vector(x[2,1,,]),as.vector(x[2,2,,]))
    type <-factor(as.vector(t ( array(dim=c(2,dim(x)[3]*dim(x)[4]),data=c("success","futility") ))))

   # int <- array(dim=c(dim(x)[4],dim(x)[3]),data=c(paste("stage",1:dim(x)[4])))
    int <- array(dim=c(dim(x)[4],dim(x)[3]),data=1:dim(x)[4])
    stage <- as.factor(c(as.vector(t(int)),as.vector(t(int))))

    number <- factor(rep(1:dim(x)[3],times=2*dim(x)[4]))

    criteria <- data.frame(type,stage, number,value,probability)
}
    else{
        criteria <- x
    }
    class(criteria) <- c("gsbCriteria","data.frame")
    return(criteria)
}

## #print methods for classes: "OCsimulation", "OCtrial", "OCprior", "OCcriteria" <<<--------------

print.gsbMainOut <- function(x,...){
  summary(x,...)
  cat("\nnames in output list: \n")
  print(names(x))
  invisible(NULL)
}








#plot methods for calsses: "OCoutput", "OCtresult", "OCdesign" <---------------------------------

## x <- x5.table
## what <- "all"
## range.delta="default"
## stages="default"


plot.gsbSimArm.result <- function(x,
                                  what,
                                  range.delta,
                                  stages,
                                  delta.grid,
                                  color,
                                  smooth,
                                  contour,...) {



  ## defaults
  i <- length(levels(x$OC$stage))
  
  ## defaults
  if(stages[1]=="default"){
    if(what=="sample size"){
      stages <- i
    }
    if(what!="sample size"){
      stages <- 1:i
    }
  }
  if(class(stages)=="numeric"){
    stages <- stages[!duplicated(stages)]
    ## stages <- stages[stages<=i]
    ## stages <- stages[stages>=0]
  }

  if (range.delta[1]=="default"){
    range.delta.control <- c(min(x$delta.grid[,1]),max(x$delta.grid[,1]))
    range.delta.treatment <- c(min(x$delta.grid[,2]),max(x$delta.grid[,2]))
  }else{
    range.delta.control <- c(range.delta[1],range.delta[2])
    range.delta.treatment <- c(range.delta[3],range.delta[4])
  }

  x2 <- x

  pan.function <- function(x, y, z, subscripts = TRUE, form = z ~ x * y, method = "loess",..., args = list(), n = smooth){
    if (length(subscripts) == 0)
      return()
    missing.x <- missing(x)
    if (!missing.x && inherits(x, "formula")) {
      form <- x
      missing.x <- TRUE
    }
    if (missing.x)
      x <- environment(form)$x
    if (missing(y))
      y <- environment(form)$y
    if (missing(z))
      z <- environment(form)$z
    x <- x[subscripts]
    y <- y[subscripts]
    z <- z[subscripts]
    ok <- is.finite(x) & is.finite(y) & is.finite(z)
    if (sum(ok) < 1)
      return()
    x <- as.numeric(x)[ok]
    y <- as.numeric(y)[ok]
    z <- as.numeric(z)[ok]
    mod <- do.call(method, c(alist(form, data = list(x = x, y = y,z = z)), args))
    lims <- current.panel.limits()
    xrange <- c(max(min(lims$x), min(x)), min(max(lims$x), max(x)))
    yrange <- c(max(min(lims$y), min(y)), min(max(lims$y), max(y)))
    xseq <- seq(xrange[1], xrange[2], length = n)
    yseq <- seq(yrange[1], yrange[2], length = n)
    zseq <- seq(min(z), max(z), length = n)
    grid <- expand.grid(x = xseq, y = yseq)
    fit <- predict(mod, grid)
    panel.levelplot(x = grid$x, y = grid$y, z = fit, subscripts = TRUE,
                    ...)
    if (delta.grid){
      grid.points(x2$delta.grid[,1],x2$delta.grid[,2], pch=3)
    }
  }

  col.l <- c("#F7FBFF", "#F4F9FE", "#F2F8FD", "#F0F7FD", "#EEF5FC", "#ECF4FB", 
             "#EAF3FB", "#E8F1FA", "#E6F0F9", "#E4EFF9", "#E2EEF8", "#E0ECF7", 
             "#DEEBF7", "#DCEAF6", "#DAE8F5", "#D8E7F5", "#D6E6F4", "#D5E5F4", 
             "#D3E3F3", "#D1E2F2", "#CFE1F2", "#CDDFF1", "#CBDEF0", "#C9DDF0", 
             "#C7DBEF", "#C5DAEE", "#C1D9ED", "#BED7EC", "#BBD6EB", "#B8D5EA", 
             "#B5D3E9", "#B1D2E7", "#AED1E6", "#ABCFE5", "#A8CEE4", "#A4CCE3", 
             "#A1CBE2", "#9ECAE1", "#9AC8E0", "#96C5DF", "#92C3DE", "#8EC1DD", 
             "#89BEDC", "#85BCDB", "#81BADA", "#7DB8DA", "#79B5D9", "#75B3D8", 
             "#71B1D7", "#6DAFD6", "#69ACD5", "#66AAD4", "#62A8D2", "#5FA6D1", 
             "#5CA3D0", "#58A1CE", "#559FCD", "#529DCC", "#4E9ACB", "#4B98C9", 
             "#4896C8", "#4493C7", "#4191C5", "#3E8EC4", "#3C8CC3", "#3989C1", 
             "#3686C0", "#3484BE", "#3181BD", "#2E7EBC", "#2C7CBA", "#2979B9", 
             "#2776B8", "#2474B6", "#2171B5", "#1F6FB3", "#1D6CB1", "#1B69AF", 
             "#1967AD", "#1764AB", "#1562A9", "#135FA7", "#115CA5", "#0F5AA3", 
             "#0D57A1", "#0B559F", "#09529D", "#084F9A", "#084D96", "#084A92", 
             "#08478E", "#08458A", "#084286", "#083F82", "#083D7E", "#083A7A", 
             "#083776", "#083572", "#08326E", "#08306B")
  
  if (!(what=="cumulative all"|| what=="all" || what=="sample size")){
    x <- x$OC
    y <- subset(x, as.numeric(x$stage)%in%stages & x$delta.control<=range.delta.control[2]& x$delta.control>=range.delta.control[1]& x$delta.treatment<=range.delta.treatment[2]& x$delta.treatment>=range.delta.treatment[1] & x$type==what)

    if (color){
      res <- levelplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=contour,main="Operating Characteristics" , col.regions=col.l)
    }else{
      res <- contourplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y,main="Operating Characteristics"  , contour=TRUE)
    }
  }
  if(what=="sample size"){
    x <- x$OC
    y <- subset(x, as.numeric(x$stage)%in%stages & x$delta.control<=range.delta.control[2]& x$delta.control>=range.delta.control[1]& x$delta.treatment<=range.delta.treatment[2]& x$delta.treatment>=range.delta.treatment[1] & x$type%in%c("sample size"))
    if (color){
      res <- levelplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=contour, col.regions=col.l, main="Expected Sample Size" )
    }else{
      res <- contourplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=TRUE, main="sample size")
    }
  }

  if(what=="all"){
    x <- x$OC
    y <- subset(x, as.numeric(x$stage)%in%stages & x$delta.control<=range.delta.control[2]& x$delta.control>=range.delta.control[1]& x$delta.treatment<=range.delta.treatment[2]& x$delta.treatment>=range.delta.treatment[1] & x$type%in%c("success","futility","success or futility"))
    if (color){
      res <- levelplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=contour, col.regions=col.l, main="Operating Characteristics" )
    }else{
      res <- contourplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=TRUE, main="Operating Characteristics" )
    }
  }

  if(what=="cumulative all"){
    x <- x$OC
    y <- subset(x, as.numeric(x$stage)%in%stages & x$delta.control<=range.delta.control[2]& x$delta.control>=range.delta.control[1]& x$delta.treatment<=range.delta.treatment[2]& x$delta.treatment>=range.delta.treatment[1] & x$type%in%c("cumulative success","cumulative futility","cumulative success or futility"))
    if (color){
      res <- levelplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=contour, col.regions=col.l, main="Operating Characteristics" )
    }else{
      res <- contourplot(value~delta.control*delta.treatment|stage+type,panel=pan.function, data=y, contour=TRUE, main="Operating Characteristics" )
    }
  }

  return(res)
}




   ## what="both"
   ## range.delta="default"
   ## stages="default"
   ## delta.grid=TRUE
   ## color=TRUE
   ## smooth=100
   ## contour=TRUE
   ## export=FALSE
   ## path=getwd()


plot.gsbSimulation <- function(x,...){
  if(class(x$truth)[2]!="matrix")
    return()
  x <- x$truth
  class(x) <- "matrix"
  plot(x, main="delta.grid")
}





## plot boundary
plot.boundary <- function(x, what){

  
  panel <- function(...) {
    panel.grid(h = -1, v = -1, lty = "dotted", lwd=2, col = "light grey")
    panel.xyplot(...)
  }

  key <- list(columns=2,space="bottom", points=list(col="black",fill=c("green4","red2"),cex=2, lwd=1.5, pch=c(24,21)), text=list(c("success / upper bound","futility / lower bound"),cex=1.1), border=FALSE)
  
  if(what=="boundary"){
    r <- xyplot(observed_difference~ stage,
                main=list("Decision criteria translated to bounds \n for the observed treatment effect",cex=1.5),
                group=x$boundary$type,
                data=x$boundary,
                type="o",
                ylab=list("Observed treatment effect",cex=1.1),
                xlab="",
                lty=2,
                lwd=2.5,
                cex=2,
                pch=c(21,24),
                col="black",
                fill=c("red2","green4"),
                panel=panel,
                key=key,
                scales=list(cex=1.1),
                aspect=.9)
  }
  if(what=="std.boundary"){
    r <- xyplot(std_observed_difference ~ stage,
                group=x$boundary$type,
                data=x$boundary,
                type="o",
                main=list("Decision criteria translated to bounds \n for the standardized observed treatment effect",cex=1.5),
                ylab=list("Standardized observed treatment effect",cex=1.1),
                xlab="",
                lty=2,
                lwd=2.5,
                cex=2,
                pch=c(21,24),
                col="black",
                fill=c("red2","green4"),
                panel=panel,
                key=key,
                scales=list(cex=1.1),
                aspect=.9)
  }
     
 return(r)
}

  




plot.gsbCalcDif.result <- function(...){
    r <- plot.gsbSimDif.result(...)
    return(r)
}
plot.gsbSimCalcDif.result <- function(...){
    plot.gsbSimDif.result(...)
}

plot.gsbDesign <- function(x,...){
    plot.gsbTrial(x$trial)
}

plot.gsbDesign <- function(x,...){
  x <- x$patients
  
    e <- x[,1]
    e2 <- x[,2]
    ne <- paste("control",names(e))
    ne2 <- paste("treatment",names(e2))
    t <- as.vector(t( cbind(e,e2)))
    n <- as.vector(t( cbind(ne,ne2)))
    names(t) <- n

    return( barplot(t,ylab="number of patients",main="Patients per Stage",...))
}


  
gsbBayesUpdate <- function(alpha,
                           beta,
                           meanData,
                           precisionData,
                           with.alpha=TRUE){

#alpha <- alpha[1,,]
#beta <- beta[1,,]
#percisionData <- b[1]
#meanData <- d2[1,,]

    bet <- beta + precisionData
    w <- beta/bet
    result <- list(beta=bet, weight=w)

    if(with.alpha){
        alp <- w*alpha + (1-w)*meanData
        result <- list(alpha=alp, beta=bet, weight=w)
    }
    return(result)
}

gsbCumulativeProbability <- function(list){
  i <- dim(list$success)[1]

  ##list <- li
  LS <- array(data=list$success,dim=c(dim(list$success),i))
  LF <- array(data=list$futility,dim=c(dim(list$success),i))
    LSF <- array(data=list$success_and_futility,dim=c(dim(list$success),i))

    llS <- array(data=NA,dim=c(i,length(LS[1,,1])))
    llF <- array(data=NA,dim=c(i,length(LS[1,,1])))
    llSF <- array(data=NA,dim=c(i,length(LS[1,,1])))

    for (de in 1:length(LS[1,,1])){
                                        #  de <- 1
        LS.t <- LS[,de,]
        LS.t[lower.tri(LS.t)] <- c(0)
        LS.t <- apply(LS.t, c(2),sum)
        llS[,de] <-LS.t

        LF.t <- LF[,de,]
        LF.t[lower.tri(LF.t)] <- c(0)
        LF.t <- apply(LF.t, c(2),sum)
        llF[,de] <-LF.t

        LSF.t <- LSF[,de,]
        LSF.t[lower.tri(LSF.t)] <- c(0)
        LSF.t <- apply(LSF.t, c(2),sum)
        llSF[,de] <-LSF.t

    }
    r <- list(cum.suc=llS,cum.fut=llF,cum.suc.fut=llSF)
    return(r)
}


gsbSampleSize <- function(design,cum.prob.suc.fut ){

    ## design <- design
    ## prob.suc.fut <- cum$cum.suc.fut

    cp <- cum.prob.suc.fut
    i <- design$nr.stages

    patients <- design$trial[,1]+design$trial[,3]

    ess.per.stage <- array(dim=dim(cp),data=NA)
    ess.per.stage[1,] <- patients[1]
    r <- ess.per.stage

    if (i>1){
        ess.per.stage[2:i,] <- patients[2:i]*(1-cp[1:(i-1),])
        ess.cum <- aperm(array(dim=c(i,dim(ess.per.stage)[2],i),data=ess.per.stage),c(1,3,2))

        for (k in 1:dim(ess.per.stage)[2]){
             ess.cum[,,k][lower.tri(ess.cum[,,k])] <- 0
        }
        r <- apply(ess.cum,c(2,3),sum)
        dimnames(r)[[1]] <- paste("stage",1:i)
    }
    return(r)
}

gsbCriteria <- function(criteria,
                        priorMean,
                        postPrecision,
                        weight){

  ## criteria=CR
  ## priorMean=prior[1]
  ## postPrecision=post$beta
  ## weight=post$weight
  
  i <- length(weight)
  nr.crit <- length(criteria[1,1,,1])
  
  if (i==1){
    ##success criteria
    if(sum(is.na(criteria[1,1,,]))==nr.crit){
      CS <- c(1e+09)
    }
    
    if(sum(is.na(criteria[1,1,,]))!=nr.crit){
      s <- criteria[1,1,,]
      ps <- criteria[2,1,,]
      qs <- qnorm(1-ps)
      
      cS <- (s-qs*postPrecision^(-0.5)-priorMean*weight)/(1-weight)
      CS <- max(cS, na.rm=TRUE)   
    }
    
    ##futility criteria
    if(sum(is.na(criteria[1,2,,]))==nr.crit){
      CF <- c(-1e+09)
    }
    if(sum(is.na(criteria[1,2,,]))!=nr.crit){
      f <- criteria[1,2,,]
      pf<- criteria[2,2,,]
      qf <- qnorm(pf)
      
      cF <- (f-qf*postPrecision^(-0.5)-priorMean*weight)/(1-weight)
      CF <- min(cF, na.rm=TRUE)
    }
  }

  if (i>1){
    
    ##boost for multiplication
    b.priorMean <- array(data=priorMean,dim=dim(t(criteria[1,1,,])))
    b.postPrecision <- array(data=postPrecision,dim=dim(t(criteria[1,1,,])))
    b.weight <- (array(data=weight,dim=dim((t(criteria[1,1,,])))))
    
    
    s <- t(criteria[1,1,,])
    ps <- criteria[2,1,,]
    qs <- t(qnorm(1-ps))
    cS <- (s-qs*b.postPrecision^(-0.5)-b.priorMean*b.weight)/(1-b.weight)
    
    f <- t(criteria[1,2,,])
    pf <- criteria[2,2,,]
    qf <- t(qnorm(pf))

    cF <- (f-qf*b.postPrecision^(-0.5)-b.priorMean*b.weight)/(1-b.weight)
    
    CS <- rep(NA,times = length(cS[,1]))
    CF <- rep(NA,times = length(cF[,1]))

    ## if one criteria per stage
    if (dim(cF)[1]==1){
      CS <- cS
      CS[is.na(CS)] <- c(1e+09)
      
      CF <- cF
      CF[is.na(CF)] <- c(-1e+09)
    }

    cs.max <- function(nr.crit, cs){
      if(sum(is.na(cs)) == nr.crit){
        CS <- c(1e+09)
      }
      if(sum(is.na(cs)) != nr.crit){
        CS <- max(cs, na.rm=TRUE)
      }      
      return(CS)
    }
    
    cf.max <- function(nr.crit, cf){
      if(sum(is.na(cf)) == nr.crit){
        CF <- c(-1e+09)
      }
      if(sum(is.na(cf)) != nr.crit){
        CF <- min(cf, na.rm=TRUE)
      }      
      return(CF)
    }
    
    
    ## if possibly more than one criteria per stage
    if (dim(cF)[1]>1){
      CS <- apply(cS, c(1), cs.max, nr.crit=nr.crit)
      CF <- apply(cF, c(1), cf.max, nr.crit=nr.crit)
    }
  }
  
  return(list(CS=CS,CF=CF))
}



gsbCalcDif <- function(design,
                       simulation){
  
  ## design <- design2
  ## simulation <- sim
  
  EX <- design$trial
  PRI <- design$prior.difference
  CR <- design$criteria
  DELTA <- simulation$truth
  i <- design$nr.stages
  
  ## PR.calc.1. Preparations
  id <- length(DELTA)
  
  ## PR.calc.1.1.1 generate b
  b <- matrix(nrow=i, ncol=1)
  b[] <- cumsum(EX[, 1]) * cumsum(EX[, 3])/(cumsum(EX[, 1]) * EX[, 4]^2 + cumsum(EX[,3])*EX[,2]^2)

  ## PR.calc.1.2 prior: genarate PRIOR with (alpha0, beta0 , nr control, nr treatment)
  if(design$prior.difference[1]=="non-informative")
    prior <- c(0,0,0,1,0,1)
  if(design$prior.difference[1]!="non-informative"){
    prior <- c(design$prior.difference[1],NA,design$prior.difference[2], design$sigma[1,1], design$prior.difference[3],design$sigma[1,2] )
    prior[2] <- (prior[4]^2/prior[3]+prior[6]^2/prior[5])^(-1)
  }
 
  ## PR.calc.2. generate all posteriors
  post <- gsbBayesUpdate(beta=prior[2],precisionData=b, with.alpha=FALSE)
  
  ## PR.calc.3. compute criteria
  cri <- gsbCriteria(criteria=CR, priorMean=prior[1], postPrecision=post$beta, weight=post$weight)
  
  test <- gsProbability(k=i, theta=(DELTA), n.I=as.vector(b), a=as.vector(cri$CF)*sqrt(b), b=as.vector(cri$CS)*sqrt(b))
  
  ## PR.calc.3. generate output frame
  if (i==1){
    value <- c(as.vector(t(test$upper$prob)),as.vector(t(test$lower$prob)),as.vector(t(test$upper$prob))+as.vector(t(test$lower$prob)),1-(as.vector(t(test$upper$prob))+as.vector(t(test$lower$prob))), rep( EX[1,1]+EX[1,3],times=length(DELTA)))
    delta <- rep(DELTA,times=i*5)
    type <- c(rep("success",times=length(DELTA)*i),rep("futility",times=length(DELTA)*i),rep("success or futility",times=length(DELTA)*i),rep("indeterminate",times=length(DELTA)*i),rep("sample size",times=length(DELTA)*i))
    stage <- rep(paste("stage",t(array(dim=c(i,length(DELTA)),data=1:i))),times=5)
    stage <- factor(stage, levels=paste("stage",1:i))
    method <- rep("numerical integration",times=length(DELTA)*i*5)
  }else{
    li <- list(success=test$upper$prob,futility=test$lower$prob,success_and_futility=test$upper$prob+test$lower$prob)
    
    cum <- gsbCumulativeProbability(list=li)
    ss <- gsbSampleSize(design=design, cum.prob.suc.fut=cum$cum.suc.fut)
    
    value <- c(as.vector(t(test$upper$prob)),as.vector(t(test$lower$prob)),as.vector(t(test$upper$prob))+as.vector(t(test$lower$prob)),1-(as.vector(t(test$upper$prob))+as.vector(t(test$lower$prob))),as.vector(t(cum$cum.suc[,])), as.vector(t(cum$cum.fut[,])), as.vector(t(cum$cum.suc.fut[,])), 1-as.vector(t(cum$cum.suc.fut[,])), as.vector(t(ss)))
    
    delta <- rep(DELTA,times=i*9)
    type <- c(rep("success",times=length(DELTA)*i),rep("futility",times=length(DELTA)*i),rep("success or futility",times=length(DELTA)*i),rep("indeterminate",times=length(DELTA)*i),rep("cumulative success",times=length(DELTA)*i),rep("cumulative futility",times=length(DELTA)*i),rep("cumulative success or futility",times=length(DELTA)*i),rep("cumulative indeterminate",times=length(DELTA)*i),rep("sample size",times=length(DELTA)*i))
    stage <- rep(paste("stage",t(array(dim=c(i,length(DELTA)),data=1:i))),times=9)
    stage <- factor(stage, levels=paste("stage",1:i))
    method <- rep("numerical integration",times=length(DELTA)*i*9)
  }
  r <- data.frame(value,type,delta,stage,method)
  
  class(r) <- c("gsbCalcDif.result","data.frame")
  
  ## prepar criteria for output
  criteria <- c(cri$CS,cri$CF)
  criteria[criteria > 1e+08] <- NA
  criteria[criteria < -1e+08] <- NA
  
  stand.criteria <- criteria*sqrt(rep(b,times=2))

  stage <- paste("stage",1:i)
  stage <- factor(stage)
  levels(stage) <- paste("stage",1:i)
  stage[1:i] <- paste("stage",1:i)
  stage <- rep(stage,2)
  
  names(criteria) <- NULL
  type <- rep(c("success / upper bound","futility / lower bound"), each=i)
  criteria <- data.frame(type,stage,observed_difference=criteria, std_observed_difference=stand.criteria)
  
  
  return(list(OC=r, boundary=criteria))
}



gsbSimDif <- function(design,
                      simulation){

    nr.of.sim <- simulation$nr.sim
    delta3 <- simulation$truth
    EX <- design$trial
    pri <- design$prior
    CR <- design$criteria

    i <- design$nr.stages

    nr.of.success.and.futility <- rep(NA,times=(i*length(delta3)))
    Psuccess.and.futility <- rep(NA,times=(i*length(delta3)))
    nr.of.success <- rep(NA,times=(i*length(delta3)))
    Psuccess <- rep(NA,times=(i*length(delta3)))
    nr.of.futility <- rep(NA,times=(i*length(delta3)))
    Pfutility <- rep(NA,times=(i*length(delta3)))

    ## prior prep.
    if(design$prior.difference[1]=="non-informative")
      prior <- c(0,0,0,1,0,1)
    if(design$prior.difference[1]!="non-informative"){
      prior <- c(design$prior.difference[1],NA,design$prior.difference[2], design$sigma[1,1], design$prior.difference[3],design$sigma[1,2] )
      prior[2] <- (prior[4]^2/prior[3]+prior[6]^2/prior[5])^(-1)
    }
    

    warning.index <- rep(0,times=length(delta3))

    for (dd in 1:length(delta3)){
  #   dd <- 1
#        print(dd)

                                        #nr os used simulated values
        nd <- rep(NA,times=i+1)
        nd[1] <- nr.of.sim
                                        #Experiment tau
        b <- matrix(nrow=i, ncol=1)
        b[] <- cumsum(EX[,1])*cumsum(EX[,3])/(cumsum(EX[,1])*EX[,4]^2+cumsum(EX[,3])*EX[,2]^2)
        tau <- 1/sqrt(b)

        ## sim z
        covmat <- diag(i)
        ut <-  sqrt(outer(c(b),c(b),"/"))
       covmat[upper.tri(covmat)] <- ut[upper.tri(ut)]
       covmat[lower.tri(covmat)] <- t(covmat)[lower.tri(covmat)]
        #covmat[upper.tri(covmat)] <- ut[upper.tri(ut)]
       

       
        z <- rmvnorm(nd[1],delta3[dd]*sqrt(b),covmat,method="svd")
        D <- sweep(z,2,sqrt(b), "/")

        tt <- gsbBayesUpdate(alpha=rep(prior[1],times=nd[1]),beta=rep(prior[2],nd[1]),precisionData=b[1],meanData=D[,1],with.alpha=TRUE)

        NS <- length(CR[1,1,,1][!is.na(CR[1,1,,1])])
        if(NS==0){
            index.SS  <- array(dim=c(1,nr.of.sim),data=FALSE)
            index.S  <- array(dim=c(1,nr.of.sim),data=FALSE)
        }
        if(!(NS==0)){
            qnorm.S <- array(dim=c(NS,nr.of.sim),data=NA)
            index.S  <- array(dim=c(NS,nr.of.sim),data=NA)

            for (jk in 1:NS){
                qnorm.S[jk,] <- qnorm(1-CR[2,1,jk,1],tt$alpha,1/sqrt(tt$beta))
                index.S[jk,] <- qnorm.S[jk,]>CR[1,1,jk,1]
            }

            index.SS <- Reduce('&',as.data.frame(t(index.S)))
        }


        NF <- length(CR[1,2,,1][!is.na(CR[1,2,,1])])
        if(NF==0){
            index.FF  <- array(dim=c(1,nr.of.sim),data=FALSE)
            index.F  <- array(dim=c(1,nr.of.sim),data=FALSE)
        }
        if(!(NF==0)){
            qnorm.F <- array(dim=c(NF,nr.of.sim),data=NA)
            index.F  <- array(dim=c(NF,nr.of.sim),data=NA)

            for (jk in 1:NF){
                qnorm.F[jk,] <- qnorm(CR[2,2,jk,1],tt$alpha,1/sqrt(tt$beta))
                index.F[jk,] <- qnorm.F[jk,]<CR[1,2,jk,1]
            }
            index.FF <- Reduce('&',as.data.frame(t(index.F)))
        }

        indeX <- index.SS | index.FF

        nr.of.success.and.futility[1+i*(dd-1)] <- sum(indeX)
           Psuccess.and.futility[1+i*(dd-1)] <- sum(indeX)/nd[1]

        nr.of.success[1+i*(dd-1)] <- sum(index.SS)
        Psuccess[1+i*(dd-1)] <- sum(index.SS)/nd[1]

        nr.of.futility[1+i*(dd-1)] <- sum(index.FF)
        Pfutility[1+i*(dd-1)] <- sum(index.FF)/nd[1]

        nd[2] <- sum(!indeX)

        if (i >1){
            for (jn in 2:i){
#                                         jn <- 1
                                        #warning message
                if (nd[jn] < simulation$warnings.sensitivity){
                    if(warning.index[dd]==0){
                        warning.index[dd] <- jn
                    }
                }


                                        #prevent from error
                if(nd[jn]<2){
                    Psuccess.and.futility[jn+i*(dd-1)] <- 0
                    Psuccess[jn+i*(dd-1)] <- 0
                    Pfutility[jn+i*(dd-1)] <- 0

                    nd[jn+1] <- 0

                }else{
                                        #select nr of used simulations
                    D <- D[!indeX,]
                                        #count the nr. of silumated values that fullfill criteria
                    tt <- gsbBayesUpdate(alpha=rep(prior[1],times=nd[jn]),beta=rep(prior[2],nd[jn]),precisionData=b[jn],meanData=D[,jn],with.alpha=TRUE)
                    #tt <- gsbBayesUpdate(alpha=tt$alpha[!indeX],beta=tt$beta[!indeX],precisionData=b[jn],meanData=d[jn,],with.alpha=TRUE)

                    NS <- length(CR[1,1,,jn][!is.na(CR[1,1,,jn])])
                    if(NS==0){
                        index.SS  <- array(dim=c(1,nd[jn]),data=FALSE)
                    }
                    if(!(NS==0)){
                        qnorm.S <- array(dim=c(NS,nd[jn]),data=NA)
                        index.S  <- array(dim=c(NS,nd[jn]),data=NA)

                        for (jk in 1:NS){
                            qnorm.S[jk,] <- qnorm(1-CR[2,1,jk,jn],tt$alpha,1/sqrt(tt$beta))

                            index.S[jk,] <- qnorm.S[jk,]>CR[1,1,jk,jn]
                    }

                        index.SS <- Reduce('&',as.data.frame(t(index.S)))
                    }

                    NF <- length(CR[1,2,,jn][!is.na(CR[1,2,,jn])])
                    if(NF==0){
                        index.FF  <- array(dim=c(1,nd[jn]),data=FALSE)
                    }
                    if(!(NF==0)){
                        qnorm.F <- array(dim=c(NF,nd[jn]),data=NA)
                        index.F  <- array(dim=c(NF,nd[jn]),data=NA)

                        for (jk in 1:NF){
                            qnorm.F[jk,] <- qnorm(CR[2,2,jk,jn],tt$alpha,1/sqrt(tt$beta))
                            index.F[jk,] <- qnorm.F[jk,]<CR[1,2,jk,jn]
                        }
                        index.FF <- Reduce('&',as.data.frame(t(index.F)))
                    }

                    indeX <- index.SS | index.FF

                    nr.of.success.and.futility[jn+i*(dd-1)] <- sum(indeX)
                    ## Psuccess.and.futility[jn+i*(dd-1)] <- mean(indeX)
                    Psuccess.and.futility[jn+i*(dd-1)] <- sum(indeX)/nd[1]

                    nr.of.success[jn+i*(dd-1)] <- sum(index.SS)
                    ## Psuccess[jn+i*(dd-1)] <- mean(index.SS)
                    Psuccess[jn+i*(dd-1)] <- sum(index.SS)/nd[1]

                    nr.of.futility[jn+i*(dd-1)] <- sum(index.FF)
                    ## Pfutility[jn+i*(dd-1)] <- mean(index.FF)
                    Pfutility[jn+i*(dd-1)] <- sum(index.FF)/nd[1]

                    nd[jn+1] <- sum(!indeX)
                }
            }
        }
    }


    liS <- array(dim=c(i,length(delta3)), data= Psuccess)
    liF <- array(dim=c(i,length(delta3)), data=Pfutility)
    liSF <- array(dim=c(i,length(delta3)), data=Psuccess.and.futility)

    li <- list(success=liS,futility=liF,success_and_futility=liSF)

                                           #output-data.frame
    if (i==1){
        sample.size <- gsbSampleSize(design=design,cum.prob.suc.fut=liSF)

        method <- factor(rep("simulation",times=i*11*length(delta3)))

        value <- c(Psuccess,Pfutility,Psuccess.and.futility,nr.of.success,nr.of.futility,nr.of.success.and.futility, as.vector((li$success[,])), as.vector((li$futility[,])), as.vector((li$success_and_futility[,])), 1-as.vector((li$success_and_futility[,])),as.vector(sample.size))

        type <- factor(c(rep("success",times=(i*length(delta3)) ),rep("futility",times=(i*length(delta3)) ),rep("success or futility",times=(i*length(delta3)) ),rep("nr of success",times=(i*length(delta3)) ),rep("nr of futility",times=(i*length(delta3)) ),rep("nr of success or futility",times=(i*length(delta3))),rep("success",times=(i*length(delta3))),rep("futility",times=(i*length(delta3))),rep("success or futility",times=(i*length(delta3))),rep("indeterminate",times=(i*length(delta3))),rep("sample size",times=length(delta3)) ))

        stage <- c(rep(paste("stage",1:i),times=length(delta3)*11 ))
        stage <- factor(stage, levels=paste("stage",1:i))
        
        delta <- c(as.vector( t(array(dim=c(length(delta3)*11,i),data=delta3))))

        r <- data.frame(value,type,delta,stage,method)
    }

    if (i>1){
        cum <- gsbCumulativeProbability(list=li)
        sample.size <- gsbSampleSize(design=design,cum.prob.suc.fut=cum$cum.suc.fut)

        method <- factor(rep("simulation",times=i*12*length(delta3)))

        value <- c(Psuccess,Pfutility,Psuccess.and.futility,1-Psuccess.and.futility,
                   nr.of.success,nr.of.futility,nr.of.success.and.futility,
                   as.vector((cum$cum.suc)),as.vector((cum$cum.fut)),as.vector((cum$cum.suc.fut)),1-as.vector((cum$cum.suc.fut)),
                   as.vector(sample.size) )

        type <- factor(c(rep("success",times=(i*length(delta3)) ),rep("futility",times=(i*length(delta3)) ),rep("success or futility",times=(i*length(delta3)) ),rep("indeterminate",times=(i*length(delta3))),
                         rep("nr of success",times=(i*length(delta3)) ),rep("nr of futility",times=(i*length(delta3)) ),rep("nr of success or futility",times=(i*length(delta3))),
                         rep("cumulative success",times=(i*length(delta3))),rep("cumulative futility",times=(i*length(delta3))),rep("cumulative success or futility",times=(i*length(delta3))),rep("cumulative indeterminate",times=(i*length(delta3))),
                         rep("sample size",times=(i*length(delta3)))))
        stage <- c(rep(paste("stage",1:i),times=length(delta3)*12 ) )
        stage <- factor(stage, levels=paste("stage",1:i))
        
        delta <- c(as.vector( t(array(dim=c(length(delta3)*12,i),data=delta3))))

        r <- data.frame(value,type,delta,stage,method)
    }

    class(r) <- c("gsbSimDif.result","data.frame")

    if(sum(warning.index)!=0){
        war <- cbind(delta3[as.logical(warning.index)],warning.index[as.logical(warning.index)])
        dimnames(war) <- list(paste("#",1:sum(as.logical(warning.index))),c("delta","incorrect in stages >= "))
        class(war) <- c(class(war),"tr")
    }else{
        war <- 0
        class(war) <- c(class(war),"fa")
    }

    return(list(OC=r,warnings=war))
}


gsbSimArm <- function(design,
                      simulation){

    ## design <- design4
    ## simulation <- simulation4

    EX <- design$trial
    CR <- design$criteria
   
    i <- design$nr.stages

    nr.of.sim  <- simulation$nr.sim
    delta.grid <- simulation$truth
    delta.p <- delta.grid[,1]
    delta.t <- delta.grid[,2]

    ld.p <- length(delta.p)
    ld.t <- length(delta.t)

    ##precision dicributions prior and trial
    if(design$prior.control[1]=="non-informative")
      design$prior.control <- c(0,0)
    if(design$prior.treatment[1]=="non-informative")
      design$prior.treatment <- c(0,0)

    prior <- c(design$prior.control[1],design$prior.treatment[1],design$prior.control[2], design$sigma[1],design$prior.treatment[2], design$sigma[2])
    sigmaDp <- prior[4]^2/prior[3]      #sigma^2 / n
    sigmaDt <- prior[6]^2/prior[5]
    names(sigmaDp) <- "sigma^2 control distr"
    names(sigmaDt) <- c("sigma^2 treatment distr")
    prior <- c(prior,sigmaDp,sigmaDt)

    tri <- array(dim=c(i,8),data=c(EX,rep(NA,times=4*i)))#
    tri[,5] <- EX[,2]^2/EX[,1]
    tri[,6] <- EX[,4]^2/EX[,3]
    tri[,7] <- 1/tri[,5]
    tri[,8] <- 1/tri[,6]

    dimnames(tri) <- list(paste("stage",1:i),c("nr.control","sigma control","nr treatment","sigma treatment","sigma^2 distr. plasebo","sigma^2 distr. treatment","precision distr. control","precision distr. treatment"))


 #test if patients =0 in stages
    test.p <- (tri[,1]&tri[,2])!=0
    test.t <- (tri[,3]&tri[,4])!=0

    Psucfut <- array(data=NA,dim=c(i,length(delta.p)))
    Psuc <- array(data=NA,dim=c(i,length(delta.p)))
    Pfut <- array(data=NA,dim=c(i,length(delta.p)))


    # index for warnings
    warning.index <- rep(0,times=ld.p)


    for (dd.t in 1:length(delta.t)){
     #   dd.t <- 1

                                        #nr of simulations
        nd <- rep(0,times=i+1)
        nd[1] <- nr.of.sim

                                        #simulate control
        dt.p <- matrix(nrow=i,ncol=nd)
        dt.t <- matrix(nrow=i,ncol=nd)


        if (i>1){
            for (k in (1:i)[test.p]){
                dt.p[k,] <- rnorm(n=nd[1],mean=delta.p[dd.t],sd=sqrt(tri[k,5]))#simulate !
            }
                d.p <- dt.p[,1:nd[1]]
            }

        if ((i==1)&test.p[1]){
            dt.p <-  rnorm(n=nd[1],mean=delta.p[dd.t],sd=sqrt(tri[1,5]))
            d.p <- t(as.matrix(dt.p))
        }


        if (i>1){
            for (k in (1:i)[test.t]){
                dt.t[k,] <- rnorm(n=nd[1],mean=delta.t[dd.t],sd=sqrt(tri[k,6]))#simulate !
            }
            d.t <- dt.t[,1:nd[1]]
        }

        if ((i==1)&test.t[1]){
            dt.t <-  rnorm(n=nd[1],mean=delta.t[dd.t],sd=sqrt(tri[1,6]))
            d.t <- t(as.matrix(dt.t))
        }

                                        #test if there are patients in stage, update if there are
        if (test.p[1]){
            post.p <- gsbBayesUpdate(alpha=rep(prior[1],times=nd[1]),beta=rep(1/prior[7],nd[1]),precisionData=tri[1,7],meanData=d.p[1,],with.alpha=TRUE)##
        }else{
            post.p <- list(alpha=rep(prior[1],times=nd[1]),beta=rep(1/prior[7],times=nd[1]))
        }

        if (test.t[1]){
            post.t <- gsbBayesUpdate(alpha=rep(prior[2],times=nd[1]),beta=rep(1/prior[8],nd[1]),precisionData=tri[1,8],meanData=d.t[1,],with.alpha=TRUE)##
        }else{
            post.t <- list(alpha=rep(prior[2],times=nd[1]),beta=rep(1/prior[8],times=nd[1]))
        }
        tt <- list(alpha=post.t$alpha-post.p$alpha,beta=(post.p$beta*post.t$beta)/(post.p$beta+post.t$beta))

        NS <- length(CR[1,1,,1][!is.na(CR[1,1,,1])])
        if(NS==0){
            index.SS  <- array(dim=c(1,nr.of.sim),data=FALSE)
            index.S  <- array(dim=c(1,nr.of.sim),data=FALSE)
        }
        if(!(NS==0)){

            qnorm.S <- array(dim=c(NS,nr.of.sim),data=NA)
            index.S  <- array(dim=c(NS,nr.of.sim),data=NA)

            for (jk in 1:NS){
                qnorm.S[jk,] <- qnorm(1-CR[2,1,jk,1],tt$alpha,1/sqrt(tt$beta))
                index.S[jk,] <- qnorm.S[jk,]>CR[1,1,jk,1]
            }

            index.SS <- Reduce('&',as.data.frame(t(index.S)))
        }


        NF <- length(CR[1,2,,1][!is.na(CR[1,2,,1])])
        if(NF==0){
            index.FF  <- array(dim=c(1,nr.of.sim),data=FALSE)
            index.F  <- array(dim=c(1,nr.of.sim),data=FALSE)
        }
        if(!(NF==0)){
            qnorm.F <- array(dim=c(NF,nr.of.sim),data=NA)
            index.F  <- array(dim=c(NF,nr.of.sim),data=NA)

            for (jk in 1:NF){
                qnorm.F[jk,] <- qnorm(CR[2,2,jk,1],tt$alpha,1/sqrt(tt$beta))
                index.F[jk,] <- qnorm.F[jk,]<CR[1,2,jk,1]
            }
            index.FF <- Reduce('&',as.data.frame(t(index.F)))
        }

        indeX <- index.SS | index.FF


        Psucfut[1,][dd.t] <- sum(indeX)/nd[1]
        Psuc[1,][dd.t] <- sum(index.SS)/nd[1]
        Pfut[1,][dd.t] <- sum(index.FF)/nd[1]

        nd[2] <- sum(!indeX)

        if (i >1){
            for (jn in 2:i){
               #  jn <- 2

                #for warning message....
                if (nd[jn]<simulation$warnings.sensitivity){
                    if(warning.index[dd.t]==0){
                        warning.index[dd.t] <- jn
                    }
                }

                                        #prevent d[jn,] incorrect dimension error
                if(nd[jn]<2){
                    Psucfut[jn,][dd.t] <- 0
                    Psuc[jn,][dd.t] <- 0
                    Pfut[jn,][dd.t] <- 0
                    nd[jn+1] <- 0
                }else{

                    d.p <- d.p[,!indeX]
                    d.t <- d.t[,!indeX]
                    if (test.p[jn]){
                        post.p <- gsbBayesUpdate(alpha=post.p$alpha[!indeX],beta=post.p$beta[!indeX],precisionData=tri[jn,7],meanData=d.p[jn,],with.alpha=TRUE)
                    }else{
                        post.p <- list(alpha=post.p$alpha[!indeX],beta=post.p$beta[!indeX])
                    }
                    if (test.t[jn]){
                        post.t <- gsbBayesUpdate(alpha=post.t$alpha[!indeX],beta=post.t$beta[!indeX],precisionData=tri[jn,8],meanData=d.t[jn,],with.alpha=TRUE)
                    }else{
                        post.t <- list(alpha=post.t$alpha[!indeX],beta=post.t$beta[!indeX])
                    }

                    tt <- list(alpha=post.t$alpha-post.p$alpha,beta=(post.p$beta*post.t$beta)/(post.p$beta+post.t$beta))

                    NS <- length(CR[1,1,,jn][!is.na(CR[1,1,,jn])])
                    if(NS==0){
                        index.SS  <- array(dim=c(1,nd[jn]),data=FALSE)
                    }
                    if(!(NS==0)){
                        qnorm.S <- array(dim=c(NS,nd[jn]),data=NA)
                        index.S  <- array(dim=c(NS,nd[jn]),data=NA)

                        for (jk in 1:NS){
                            qnorm.S[jk,] <- qnorm(1-CR[2,1,jk,jn],tt$alpha,1/sqrt(tt$beta))
                            index.S[jk,] <- qnorm.S[jk,]>CR[1,1,jk,jn]
                        }

                        index.SS <- Reduce('&',as.data.frame(t(index.S)))
                    }


                    NF <- length(CR[1,2,,jn][!is.na(CR[1,2,,jn])])
                    if(NF==0){
                        index.FF  <- array(dim=c(1,nd[jn]),data=FALSE)
                    }
                    if(!(NF==0)){
                        qnorm.F <- array(dim=c(NF,nd[jn]),data=NA)
                        index.F  <- array(dim=c(NF,nd[jn]),data=NA)

                        for (jk in 1:NF){
                            qnorm.F[jk,] <- qnorm(CR[2,2,jk,jn],tt$alpha,1/sqrt(tt$beta))
                            index.F[jk,] <- qnorm.F[jk,]<CR[1,2,jk,jn]
                        }
                        index.FF <- Reduce('&',as.data.frame(t(index.F)))
                    }

                    indeX <- index.SS | index.FF

                    nd[jn+1] <- sum(!indeX)

                    Psucfut[jn,][dd.t] <- sum(indeX)/nd[1]

                    Psuc[jn,][dd.t] <- sum(index.SS)/nd[1]

                    Pfut[jn,][dd.t] <- sum(index.FF)/nd[1]
                }
            }
        }

    }

    if (i==1){
        value <- c(as.vector(t(Psuc)),as.vector(t(Pfut)),as.vector(t(Psucfut)),1-as.vector(t(Psucfut)), rep( EX[1,1]+EX[1,3],times=length(delta.grid[,1])))

        type <- c(rep("success",times=ld.p*i),rep("futility",times=ld.p*i),rep("success or futility",times=ld.p*i),rep("indeterminate",times=ld.p*i),rep("sample size",times=ld.p*i))

        delta.control <- rep(delta.p,times=i*5)
        delta.treatment <- rep(delta.t,times=i*5)
        delta <- delta.treatment-delta.control

        stage <- as.vector(t(array(data=paste("stage",1:i),dim=c(i,ld.t))))
        stage <- factor(stage, levels=paste("stage",1:i))
        
        method <-  as.vector(rep("simArm",times=5*ld.t))

        r <- data.frame(value,delta.control,delta.treatment,delta,type,stage,method)
    }else{

        li <- list(success=Psuc,futility=Pfut,success_and_futility=Psucfut)
        cum <- gsbCumulativeProbability(list=li)
        ss <- gsbSampleSize(design=design,cum.prob.suc.fut=cum$cum.suc.fut)

        value <- c(as.vector(t(Psuc)),as.vector(t(Pfut)),as.vector(t(Psucfut)),1-as.vector(t(Psucfut)), as.vector(t(cum$cum.suc)),as.vector(t(cum$cum.fut)), as.vector(t(cum$cum.suc.fut)),1-as.vector(t(cum$cum.suc.fut)), as.vector(t(ss)))

        type <- c(rep("success",times=ld.p*i),rep("futility",times=ld.p*i),rep("success or futility",times=ld.p*i),rep("indeterminate",times=ld.p*i),rep("cumulative success",times=ld.p*i),rep("cumulative futility",times=ld.p*i),rep("cumulative success or futility",times=ld.p*i),rep("cumulative indeterminate",times=ld.p*i),rep("sample size",times=ld.p*i))

        delta.control <- rep(delta.p,times=i*9)
        delta.treatment <- rep(delta.t,times=i*9)
        delta <- delta.treatment-delta.control

        stage <- rep(as.vector(aperm(array(data=paste("stage",1:i),dim=c(i,ld.t)),c(2,1))),times=9)
        stage <- factor(stage, levels=paste("stage",1:i))
        
        method <-  rep("simArm",length(value))

        r <- data.frame(value,delta.control,delta.treatment,delta,type,stage,method)

    }

    class(r) <- c("gsbSimArm.result","data.frame")#noDif

    if(sum(warning.index)!=0){
        war <- cbind(delta.t[as.logical(warning.index)],delta.p[as.logical(warning.index)],warning.index[as.logical(warning.index)])
        dimnames(war) <- list(paste("# ",1:sum(as.logical(warning.index))),c("delta.treatment","delta.control","incorrect in stages >="))
        class(war) <- c(class(war),"tr")
    }else{
        war <- 0
        class(war) <- c(class(war),"fa")
    }

    return(list(OC=r,warnings=war,delta.grid=delta.grid))
}




gsb <- function(design=NULL,
                simulation=NULL){

  ## set.seed
  if(!is.null(simulation$seed))
    set.seed(simulation$seed)

  ## check NULL
  if(is.null(design))
    stop("Please specify the argument \"design\". The design object can be created with function gsbDesign().")

  if(is.null(simulation))
    stop("Please specify the argument \"simulation\". The simulation object can be created with function gsbSimulation.")

  ## check class
  if(class(design)[1]!= "gsbDesign")
    stop("The argument \"design\" has to be of class \"gsbDesign\". The design object can be created with function gsbDesign().")

  if(class(simulation)[1] != "gsbSimulation")
    stop("The argument \"design\" has to be of class \"gsbSimulation\". The simulation object can be created with function gsbSimulation().")
    
  ## check if design is compatible with Simulation
  if(simulation$type.update[1]=="treatment effect" & (design$prior.control[1]!="non-informative" |design$prior.treatment[1]!="non-informative")  )
    stop("If the argument \"simulation\" is specified for a update on \"treatment effect\", the prior in the design object has to be specified on difference too. (use argument \"prior.difference\" of function \"gsbDesign()\".")

  if(simulation$type.update=="per arm" & design$prior.difference[1]!="non-informative")
    stop("If the argument \"simulation\" is specified for a update \"per arm\", the prior in the design object has to be specified per arm too. (Use arguments \"prior.control\" and \"prior.treatment\" of function \"gsbDesign()\".")
  

  ## convert to old format (end complete)
  ## ---------------------------------------------------------------------------
  designInput <- design


  design$sigma <- matrix(design$sigma, nrow=design$nr.stages, ncol=2, byrow=TRUE)
  design$trial <- cbind(design$patients[,1],design$sigma[,1],design$patients[,2], design$sigma[,2])

 
  if(class(simulation)[2]=="gsbSimArm"){
    
    s <- system.time(r <- gsbSimArm(design = design,
                                    simulation = simulation))
    if(class(r$warnings)[2]=="tr"){
      message(simulation$nr.sim," simlations per (delta.control,delta.treatment) is smaller than the specified threshold = ",simulation$warnings.sensitivity ," for the following delta.control and delta.treatment. Further information in help of \"gsb()\": \n" )
      w <- as.matrix(r$warnings)
      class(w) <- "matrix"
      print(w, digits=3)
    }
    res <- list(OC = r$OC, design=designInput, simulation=simulation, delta.grid=r$delta.grid, warnings=r$warnings,  system.time=s)
    
  }
  
  if(class(simulation)[2]=="gsbSimDif"){
   
    s <- system.time(r <- gsbSimDif(design = design,
                                    simulation = simulation
                                    )
                     )
    if(class(r$warnings)[2]=="tr")
      message(simulation$nr.sim," simlations per (delta.control,delta.treatment) is smaller than the specified threshold = ",simulation$warnings.sensitivity ," for the following delta. Further information in help of \"gsb()\": \n" )
    w <- as.matrix(r$warnings)
    class(w) <- "matrix"
    print(w,digits=3)
    res <- list(OC = r$OC, design=designInput, simulation=simulation, warnings=r$warnings,  system.time=s)
  }
  
  if(class(simulation)[2]=="gsbCalcDif"){
    s <- system.time(r <- gsbCalcDif(design = design,
                                     simulation = simulation
                                     )
                     )
    res <- list(OC = r$OC, boundary=r$boundary, design=designInput, simulation=simulation,  system.time=s)
  }
  
  if(class(simulation)[2]=="gsbSimCalcDif"){
   
    both <- function(design = design,
                     simulation = simulation){
      r1 <- gsbSimDif(design = design,
                      simulation = simulation)
      
      r2 <- gsbCalcDif(design = design,
                       simulation = simulation)
      r <- list(r1=r1,r2=r2)
      return(r)
    }
    
    s <- system.time(r <- both(design = design,
                               simulation = simulation))
    
    r1 <- r$r1
    r2 <- r$r2
    r <- rbind(r1$OC,r2$OC)
   
    class(r) <- c("gsbSimCalcDif.result","data.frame")
    
    if(class(r1$warnings)[2]=="tr"){
      cat(simulation$nr.sim," simlations per (delta.control,delta.treatment) is smaller than the specified threshold = ",simulation$warnings.sensitivity ," for the following delta. Further information in help of \"gsb()\": \n" )
      w <- as.matrix(r1$warnings)
      class(w) <- "matrix"
      print(w, digits=3)
    }
    res <- list(OC = r, boundary=r2$boundary, design=designInput, simulation=simulation, warnings=r1$warnings,  system.time=s)
  }
  
  class(res) <- c("gsbMainOut" ,class(res))
  return(res)
}


  ## nr.stages=1
  ## patients= t(c(1,2))
  ## sigma=1
  ## criteria.success=c(0,.8)
  ## criteria.futility=c(0.8)
  ## prior.difference=c(1,2,3,4)
  ## prior.control="non-informative"
  ## prior.treatment="non-informative"



gsbDesign <- function(nr.stages=NULL,
                      patients=NULL,
                      sigma=1,
                      criteria.success=NULL,
                      criteria.futility=NULL,
                      prior.difference="non-informative",
                      prior.control="non-informative",
                      prior.treatment="non-informative")
{
  
  ## integrer2numeric
  nr.stages <- integer2numeric(nr.stages)
  patients <- integer2numeric(patients)
  sigma <- integer2numeric(sigma)
  criteria.success <- integer2numeric(criteria.success)
  criteria.futility <- integer2numeric(criteria.futility)
  prior.difference <- integer2numeric(prior.difference)
  prior.control <- integer2numeric(prior.control)
  prior.treatment <- integer2numeric(prior.treatment)
  
  ## nr.stages
  if(is.null(nr.stages))
    stop("Please specify the argument \"nr.stges\".")

  if(class(nr.stages)!="numeric")
    stop("Argument \"nr.stages\" has to be of class \"numeric\". \n\n")
  
  if(length(nr.stages)!=1)
    stop("Argument \"nr.stages\" has to be a numeric of length 1. \n\n")
  
  if(round(nr.stages, digits=0)!=nr.stages)
    stop("Argument \"nr.stages\" has to be a numeric containing an integer value. \n\n")
  
  if(nr.stages<0)
    stop("Argument \"nr.stages\" has to be positive. \n\n")

  names(nr.stages) <- "nr of stages"
  
  ## patients
  if(is.null(patients))
    stop("Please specify the argument \"patients\".")

 
  if(class(patients)!="numeric" & class(patients)!="matrix")
    stop("Argument \"patients\" has to be of class \"numeric\" or \"matrix\".\n")
  
  if(class(patients)=="numeric"){
    
    if(length(patients)>2)
      stop("Argument \"patients\" has to be a vector of length 1 or 2 or a matrix of size \"nr.stages x 2\"  or size \"nr.stages x 1\".")
    
    if(length(patients)==1)
      patients <- c(patients,patients)
    
    if(patients[1]<0 | patients[2]<0)
      stop("Argument \"patients\" has to be non-negative.")
    
    patients <- matrix(patients,nrow=nr.stages, ncol=2, byrow=TRUE)  
    
  }else{
    
    if(dim(patients)[2]==1)
      patients <- cbind(patients,patients)
    
    if(dim(patients)[2]>2 | dim(patients)[1] != nr.stages)    
      stop("Argument \"patients\" has to be a vector of length 1 or 2 or a matrix of size \"nr.stages x 2\" or size \"nr.stages x 1\".")
  }
  
  dimnames(patients) <- list(paste("stage",1:nr.stages),c("control", "treatment"))
  
  

  
  ## sigma
  if(is.null(sigma))
    stop("Please specify the argument \"sigma\".")
  
  
  if(class(sigma)!="numeric" & class(sigma)!="matrix")
    stop("Argument \"sigma\" has to be of class \"numeric\" or \"matrix\".\n")
  
  if(class(sigma)=="numeric"){
    
    if(length(sigma)>2)
      stop("Argument \"sigma\" has to be a vector of length 1 or 2. It also can be a matrix of size \"nr.stages x 2\" or size \"nr.stages x 1\"")
    
    if(length(sigma)==1)
      sigma <- c(sigma,sigma)
    
    if(sigma[1]<0 | sigma[2]<0)
      stop("Argument \"sigma\" has non-negative.")
    
    names(sigma) <- c("control", "treatment")
  }
  
  if(class(sigma)=="matrix"){
    if(dim(sigma)[2]==1)
      sigma <- cbind(sigma,sigma)
    
    if(dim(sigma)[2]>2 | dim(sigma)[1] != nr.stages)    
      stop("Argument \"sigma\" has to be a vector of length 1 or 2. It also can be a matrix of size \"nr.stages x 2\" or size \"nr.stages x 1\"")
    
    dimnames(sigma) <- list(paste("stage",1:nr.stages),c("control", "treatment"))
  }
  
  ## prior
  on.diff <- NULL
  if(is.null(prior.difference) & is.null(prior.control) & is.null(prior.treatment))
    stop("Please specify the argument \"prior.difference\" or the arguments \"prior.control\" and \"prior.treatment\"")
  
  if((prior.difference[1] !="non-informative" & prior.control[1] !="non-informative") | (prior.difference[1] !="non-informative" & prior.treatment[1] !="non-informative"))
    stop("Please specify the argument \"prior.difference\" OR the arguments \"prior.control\" and \"prior.treatment\"")

  if(!(prior.difference[1] =="non-informative" & prior.control[1] =="non-informative"  & prior.treatment[1] =="non-informative") ){
    
    
    if(prior.difference[1] !="non-informative"){
      
      if(is.null(prior.difference))
        stop("Please specify the argument \"prior.difference\".")
      
      if(class(prior.difference) != "numeric")
        stop("the argument \"prior.difference\" has to be a \"numeric\" of length 3.")
      
      if(length(prior.difference) != 3)
        stop("the argument \"prior.difference\" has to be a \"numeric\" of length 3.")
      
      if(sum(prior.difference[2:3] > 0)!=2)
        stop("the entries 2 and 3 of the argument \"prior.difference\" have to strictly positive.")
      
      names(prior.difference) <- c("difference", "n_control" , "n_treatment")
      
    }
    
    
    if(prior.control[1] !="non-informative"){
      
      if(class(prior.control)!="numeric")
        stop("the argument \"prior.control\" has to a \"numeric\" of length 2.")
      if(length(prior.control)!=2)
        stop("the argument \"prior.control\" has to a \"numeric\" of length 2.")        
      if(prior.control[2] < 0)
        stop("the entry 2 of the argument \"prior.control\" has to > 0.")
      
      names(prior.control) <- c("mu_contol", "n_control")
    }
    
    if(prior.treatment[1] !="non-informative"){
      if(class(prior.treatment)!="numeric")
        stop("the argument \"prior.treatment\" has to be a \"numeric\" of length 2.")
      if(length(prior.treatment)!=2)
        stop("the argument \"prior.treatment\" has to be a \"numeric\" of length 2.")        
      if(prior.treatment[2] < 0)
        stop("the entry 2 of the argument \"prior.treatment\" has to > 0.")
      
      names(prior.treatment) <- c("mu_contol", "n_treatment")
    }
    
  }
  

  
  ## criteria
  if(is.null(criteria.success))
    stop("Please specify the argument \"criteria.success\".")

  if(is.null(criteria.futility))
    stop("Please specify the argument \"criteria.futility\".")

  if(length(criteria.success)==1 & is.na(criteria.success[1]))
    criteria.success <- rep(NA,2)

  if(length(criteria.futility)==1 & is.na(criteria.futility[1]))
    criteria.futility <- rep(NA,2)

  if(!(class(criteria.futility)=="numeric" | class(criteria.futility)=="matrix") & !is.na(criteria.futility[1])) 
    stop("The argument \"criteria.futility\" has to be of class \"numeric\" or \"matrix\".")

  if(!(class(criteria.success)=="numeric" | class(criteria.success)=="matrix") & !is.na(criteria.success[1])) 
    stop("The argument \"criteria.success\" has to be of class \"numeric\" or \"matrix\".")

  
  ## check maximal nr of criteria
  if(class(criteria.success)=="numeric"){
    if(length(criteria.success)%%2) ## check if length is odd.
      stop(" Argument \"criteria.success\" has to contain pairs, i.e. length(criteria.success) %% 2 has to be 0.")
    if(sum(criteria.success[ (1:(length(criteria.success)/2)) *2]>1, na.rm=TRUE)!=0 |sum(criteria.success[ (1:(length(criteria.success)/2))*2]<0, na.rm=TRUE)!=0)
      stop("Every second value in argument \"criteria.success\" is a probability and has to be in (0,1).")
    
    names(criteria.success) <- rep(c("value","probability"), length(criteria.success)/2)
  }
  if(class(criteria.success)=="matrix"){
    if(dim(criteria.success)[2]%%2 ||dim(criteria.success)[1] != nr.stages )
      stop("If \"criteria.success\" is inputed as matrix, the dimension has to be \"c(nr.stages, 'odd number')\"")
    if(sum(criteria.success[ ,(1:(length(criteria.success[1,])/2)) *2]>1, na.rm=TRUE)!=0 |sum(criteria.success[ ,(1:(length(criteria.success[1,])/2))*2]<0, na.rm=TRUE)!=0)
      stop("Every second colum in argument \"criteria.success\" is a probability and has to be in (0,1).")
    dimnames(criteria.success)[[2]] <- rep(c("value","probability"), dim(criteria.success)[2]/2)
    dimnames(criteria.success)[[1]] <- paste("stage", 1:nr.stages)
  }
  if(class(criteria.futility)=="numeric"){
    if(length(criteria.futility)%%2) ## check if length is odd.
      stop("\"criteria.futility\" has to contain pairs, i.e. length(criteria.futility) %% 2 has to be 0. \n\n")
    if(sum(criteria.futility[ (1:(length(criteria.futility)/2)) *2]>1, na.rm=TRUE)!=0 |sum(criteria.futility[ (1:(length(criteria.futility)/2))*2]<0, na.rm=TRUE)!=0)
      stop("Every second value in argument \"criteria.futility\" is a probability and has to be in (0,1).")

    names(criteria.futility) <- rep(c("value","probability"), length(criteria.futility)/2)
  }
  if(class(criteria.futility)=="matrix"){
    if(dim(criteria.futility)[2]%%2 ||dim(criteria.futility)[1] != nr.stages )
      stop("If \"criteria.futility\" is inputed as matrix, the dimension has to be \"c(nr.stages, 'odd number')\" \n\n")
    if(sum(criteria.futility[ ,(1:(length(criteria.futility[1,])/2)) *2]>1, na.rm=TRUE)!=0 |sum(criteria.futility[, (1:(length(criteria.futility[1,])/2))*2]<0, na.rm=TRUE)!=0)
      stop("Every second colum in argument \"criteria.futility\" is a probability and has to be in (0,1).")

    dimnames(criteria.futility)[[2]] <- rep(c("value","probability"), dim(criteria.futility)[2]/2)
    dimnames(criteria.futility)[[1]] <- paste("stage", 1:nr.stages)
  }

  if(class(criteria.success)[1]=="numeric")
    ns <- length(criteria.success)/2
  if(class(criteria.success)[1]=="matrix")
    ns <- dim(criteria.success)[2]/2
  if(class(criteria.success)[1]=="logical")
    ns <- 1

  if(class(criteria.futility)[1]=="numeric")
    nf <- length(criteria.futility)/2 
  if(class(criteria.futility)[1]=="matrix")
    nf <- dim(criteria.futility)[2]/2
  if(class(criteria.futility)[1]=="logical")
    nf <- 1
  

  m <- max(nf,ns) # maximal length of success and fut crit
  
  ## array and labels
  lab.criteria.1 <- c("value", "prob")
  lab.criteria.2 <- c("success","futility")
  lab.criteria.3 <- paste("number",1:max(ns,nf))
  lab.criteria.4 <- paste("stage",1:nr.stages)
  criteria <- array(dim=c(2,2,m,nr.stages), dimnames=list(lab.criteria.1,lab.criteria.2,lab.criteria.3,lab.criteria.4))
  
  ## fill criteria array
  if(class(criteria.success)[1]!="matrix"){
    if(ns < m){
      criteria.success <- c(criteria.success,rep(NA,times=2*(m-ns)))
    }
    criteria[,1,,] <- t(matrix(data=criteria.success, ncol= m ,nrow=2,byrow=TRUE))
  }else{
    if(ns < m)
      criteria.success <- cbind(criteria.success, matrix(nrow=nr.stages,ncol=2*(m-ns), data=NA) )
    criteria[,1,,] <- t(criteria.success)
  }
  if(class(criteria.futility)[1]!="matrix"){
    if(nf < m){
      criteria.futility <- c(criteria.futility,rep(NA,times=2*(m-nf)))
    }
    criteria[,2,,] <- t(matrix(data=criteria.futility, ncol= m ,nrow=2,byrow=TRUE))
  }else{
    if(nf < m)
      criteria.futility <- cbind(criteria.futility, matrix(nrow=nr.stages,ncol=2*(m-nf), data=NA) )
    criteria[,2,,] <- t(criteria.futility)
  }
    
  class(patients)=c("gsbPatients",class(patients))
  class(sigma)=c("gsbSigma",class(sigma))
  class(criteria)=c("gsbCriteria",class(criteria))
   class(prior.difference)=c("gsbPrior",class(prior.difference))
  class(prior.control)=c("gsbPrior",class(prior.control))
  class(prior.treatment)=c("gsbPrior",class(prior.treatment))

  if (is.null(on.diff))
    design <- list(nr.stages=nr.stages, patients=patients, sigma=sigma, criteria=criteria, prior.difference=prior.difference,prior.control=prior.control, prior.treatment=prior.treatment)
  if (!is.null(on.diff)){
    if(on.diff)
      design <- list(nr.stages=nr.stages, patients=patients, sigma=sigma, criteria=criteria, prior.difference=prior.difference,prior.control=NULL, prior.treatment=NULL)
    if(!on.diff)
      design <- list(nr.stages=nr.stages, patients=patients, sigma=sigma, criteria=criteria, prior.difference=NULL, prior.control=prior.control, prior.treatment=prior.treatment)
  }
  class(design) <- c("gsbDesign","list")
 
   
  return(design)
}

print.gsbDesign <- function(x, ...){
  cat("\n*** Trial Design ***\n\n")
  
  cat("number of stages: ",x$nr.stages,"\n\n")
   
  if(x$prior.difference[1]!="non-informative"){
    cat("prior difference:\t ", x$prior.difference[1],"\n")    
    cat("prior patients control:\t ", x$prior.difference[2],"\n")
    cat("prior patients treatment:", x$prior.difference[3],"\n")
  }
  if(x$prior.control[1]!="non-informative"){
    if(x$prior.control[1]=="non-informative")
        cat("prior control:\t\t ", x$prior.control[1],"\n")    
    if(x$prior.control[1]!="non-informative"){
      cat("prior mean control: \t ", x$prior.control[1],"\n")
      cat("prior patients control:  ", x$prior.control[2],"\n")
    }
  }
  if(x$prior.treatment[1]!="non-informative"){
    if(x$prior.treatment[1]=="non-informative")
      cat("prior treatment:\t ", x$prior.treatment[1],"\n")    
    if(x$prior.treatment[1]!="non-informative"){
      cat("prior mean treatment:\t ", x$prior.treatment[1],"\n")
      cat("prior patients treatment:", x$prior.treatment[2],"\n")
    }
  }
  if(x$prior.difference[1]=="non-informative" & x$prior.control[1]=="non-informative" & x$prior.treatment[1]=="non-informative") 
     cat("prior:  non-informative\n")
 
  cat("\npatients: \n")
  class(x$patients) <- "matrix"
  print(x$patients)
  cat("\nsigma control:", x$sigma[1], ";\tsigma treatment:", x$sigma[2],"\n")
  cat("\ncriteria: \n")
  print(na.omit(as.data.frame(x$criteria)),row.names=FALSE)
  invisible(NULL)
}


gsbSimulation <- function(truth=NULL,
                          type.update = c("treatment effect","per arm"),
                          method = c("numerical integration", "simulation", "both"),
                          grid.type= c("table","plot","sliced","manually"),
                          nr.sim= 50000,
                          warnings.sensitivity = 100,
                          seed=NULL){

  type.update <- match.arg(type.update)
  method <- match.arg(method)
  grid.type <- match.arg(grid.type)

  
  if(is.null(truth))
    stop("Please specify the argument \"truth\".")
  truth <- integer2numeric(truth)
  nr.sim <- integer2numeric(nr.sim)
  warnings.sensitivity <- integer2numeric(warnings.sensitivity)

  

  if(class(warnings.sensitivity)!="numeric")
    stop("\"warnings.sensitivity\" has to be a vector of length 1 containing a positive integer value.")
  if( length(warnings.sensitivity)>1)
    stop("\"warnings.sensitivity\" has to be a vector of length 1 containing a positive integer value.")
  if( round(warnings.sensitivity,digits=0)!=warnings.sensitivity | warnings.sensitivity<1)
    stop("\"warnings.sensitivity\" has to be a vector of length 1 containing a positive integer value.")



  if(class(nr.sim)!="numeric" | length(nr.sim)>1)
    stop("\"nr.sim\" has to be a numeric of length 1 containing a positive integer value.")
  if(round(nr.sim,digits=0)!=nr.sim | nr.sim<warnings.sensitivity)
    stop("\"nr.sim\" has to be a numeric of length 1 containing a integer value > \"warnings.sensitivity\".")

 
  ## seed
  if(!is.null(seed)){
    if(seed[1] == "generate"){
      seed <- round(runif(1,-1000,10000), digits=0)
    }
    if(seed[1] != "generate"){
      if(length(seed)>1 | class(seed) != "numeric" )
        stop("seed has to be a vector with one integer value.")
      seed <- round(seed, digits=1)
    }
  }
  

  ## if (method=="calc treatment effect"){
  if (type.update[1]=="treatment effect" && method == "numerical integration"){
    if (grid.type=="manually"){
      if(class(truth) !="numeric")
        stop("\"truth\" has to be a numeric if type.update = \"treatment effect\" and grid.type = \"manually\".")
      r <- list(truth=truth,type.update=type.update, method=method, grid.type=grid.type)
    }else{
      if(class(truth) !="numeric" | length(truth)!=3 | truth[2]<truth[1])
        stop("\"truth\" has to be a numeric of length 3 specifying a sequence as in seq() if type.update = \"treatment effect\" and grid.type != \"manually\".")
      r <- list(truth=seq(truth[1],truth[2],length.out=truth[3]),type.update=type.update,method=method, grid.type=grid.type, seed=seed)
    }
    class(r) <- c("gsbCalcDif","list")
  }
  
  
  ##if(type=="sim treatment effect"){
  if (type.update=="treatment effect" && method == "simulation"){
    
    if (grid.type=="manually"){
      if(class(truth) !="numeric")
        stop("\"truth\" has to be a numeric if type.update = \"treatment effect\" and grid.type = \"manually\".")
    }
    
    if (grid.type!="manually"){
      if(class(truth) !="numeric" | length(truth)!=3 | truth[2]<truth[1])
        stop("\"truth\" has to be a numeric of length 3 specifying a sequence as in seq() if type.update = \"treatment effect\" and grid.type != \"manually\".")
      truth <- seq(truth[1],truth[2],length.out=truth[3])
    }
    
    r <- list(truth=truth,
              type.update=type.update,
              method=method,
              grid.type=grid.type,
              nr.sim=nr.sim,
              seed=seed,
              warnings.sensitivity=warnings.sensitivity)
    class(r) <- c("gsbSimDif","list")
  }
  
  ##if(type=="sim and calc treatment effect"){
  
  if (type.update=="treatment effect" && method == "both"){

    if (grid.type=="manually"){
      if(class(truth) !="numeric")
      stop("\"truth\" has to be a numeric if type.update = \"treatment effect\" and grid.type = \"manually\".")
     }

    if (grid.type!="manually"){
      if(class(truth) !="numeric" | length(truth)!=3 | truth[2]<truth[1])
        stop("\"truth\" has to be a numeric of length 3 specifying a sequence as in seq() if type.update = \"treatment effect\" and grid.type != \"manually\".")
      truth <- seq(truth[1],truth[2],length.out=truth[3])
    }
    r <- list(truth=truth,
              type.update=type.update,
              method=method,
              grid.type=grid.type,
              nr.sim=nr.sim,
              seed=seed,warnings.sensitivity=warnings.sensitivity)
    class(r) <- c("gsbSimCalcDif","list")
  }
  
  
  if (type.update=="per arm"){
    if (grid.type=="manually"){
      if(class(truth)[1] != "matrix")
        stop("\"truth\" has to be a n x 2 matrix if type.update = \"per arm\" and grid.type = \"manually\".")
      if(dim(truth)[2] != 2)
        stop("\"truth\" has to be a n x 2 matrix if type.update = \"per arm\" and grid.type = \"manually\".")
    }

    if (grid.type=="table"){
      if(class(truth) != "list" )
        stop("If type.update = \"per arm\" and grid.type = \"table\", \"truth\" has to be a list where the first obect is a vector containing the control values and in the second object is a vector containing the treatment values. These are use to build a regular grid.")
      if(length(truth) != 2)
        stop("If type.update = \"per arm\" and grid.type = \"table\", \"truth\" has to be a list where the first obect is a vector containing the control values and in the second object is a vector containing the treatment values. These are use to build a regular grid.")
      if(class(integer2numeric(truth[[1]])) != "numeric" | class(integer2numeric(truth[[2]])) != "numeric" )
        stop("If type.update = \"per arm\" and grid.type = \"table\", \"truth\" has to be a list where the first object is a numeric containing the control values and the second object is a numeric containing the treatment values. These are use to build a regular grid.")

      truth <- as.matrix(expand.grid(truth[[1]],truth[[2]]))
    }

    if (grid.type=="plot"){
      if(class(truth) != "numeric" | length(truth)!=5)
        stop("If type.update = \"per arm\" and grid.type = \"plot\", \"truth\" has to be a numeric of length 5 specifying c(control.range.min, control.range.max, treatment.range.min, treatment.range.max, nr.of.points). ")
   
      bnds <- rbind(c(truth[1],truth[2]), c(truth[3],truth[4]))
      truth <- rbind(gsbCgrid(truth[5]-4,bnds),
                     c(truth[1],truth[3]),
                     c(truth[1],truth[4]),
                     c(truth[2],truth[3]),
                     c(truth[2],truth[4]))                     
    }

    if (grid.type=="sliced"){
      if(class(truth) != "list" | length(truth)!=2)
        stop("If type.update = \"per arm\" and grid.type = \"sliced\", \"truth\" has to be a list of length 2 containing two vectors. The first specifies the control values. The second the differences (treatment - control) of interest. ")

      treat <- truth[[1]] + rep(truth[[2]], each=length(truth[[1]]))
      truth <- cbind(truth[[1]], treat)
      dimnames(truth)[[2]] <- c("control","treatment")
                         
    }

   
    dimnames(truth) <- NULL
    dimnames(truth)[[2]] <- c("control","treatment")
    class(truth) <- c("gsbDelta.grid",class(truth))
    
    r <- list(truth=truth,
              type.update=type.update,
              method="simulation",
              grid.type=grid.type,
              nr.sim=nr.sim,
              seed=seed,
              warnings.sensitivity=warnings.sensitivity)

    class(r) <- c("gsbSimArm","list")
    class(r) <- c("gsbSimulation",class(r))

    return(r)
  }
  
  class(r) <- c("gsbSimulation",class(r))
  return(r)
}


print.gsbSimulation <- function(x,...){
  cat("\n*** Simulation Settings ***\n\n")

  cat("type.update: \t\t",x$type.update,"\n")  
  if(x$method=="both")
    cat("method: \t\t numerical integration & simulation\n")
  if(x$method!="both")
    cat("method: \t\t",x$method, "\n")  

  if(x$method!="numerical integration"){
    cat("nr of simulations: \t",x$nr.sim, "\n")  
    cat("warnings.sensitivity: \t",x$warnings.sensitivity, "\n")  
  }
  if(x$type.update=="per arm")
    cat("grid.type: \t\t",x$grid.type,"\n")  

  if(!is.null(x$seed))
    cat("seed: \t\t\t",x$seed,"\n")  

  if(length(class(x$truth))==2){ ## i.e. per arm
    if(x$grid.type=="table"){
      cat("\ngrid of true control arm and treatment arm values:\n")
      cat("- ",length(x$truth[,1][!duplicated(x$truth[,1])])," distinct control arm values from ", min(x$truth[,1]), " to ", max(x$truth[,1]) , ".\n",sep="")
      
      cat("- ",length(x$truth[,2][!duplicated(x$truth[,2])])," distinct treatment arm values from ", min(x$truth[,2]), " to ", max(x$truth[,2]) , ".\n",sep="") 
    }

    if(x$grid.type=="manually"){
      cat("\ngrid of true control arm and treatment arm values:\n")
      cat("- control arm values from ", min(x$truth[,1]), " to ", max(x$truth[,1]) , ".\n",sep="")
      cat("- tretment arm from ", min(x$truth[,2]), " to ", max(x$truth[,2]) , ".\n",sep="") 
      cat("- in total ",length(x$truth[,1])," points.\n", sep="")
    }

    if(x$grid.type=="plot"){
      cat("\ngrid of true control arm and treatment arm values:\n")
      cat("- ",length(x$truth[,1])," distinct points.\n", sep="")
      cat(" control arm values from ", round(min(x$truth[,1]),3), " to ", round(max(x$truth[,1]),3) , ".\n",sep="")
      cat("- treatment arm values from ", round(min(x$truth[,2]),3), " to ", round(max(x$truth[,2]),3) , ".\n",sep="") 
    }
    if(x$grid.type=="sliced"){
      cat("\ngrid of true control arm and treatment arm values:\n")
      cat("- ",length(x$truth[,1])," distinct points.\n", sep="")
      cat("- control arm values from ", min(x$truth[,1]), " to ", max(x$truth[,1]) , ".\n",sep="")
      cat("- delta values from ", round(min(x$truth[,2]-x$truth[,1]),3), " to ", max(x$truth[,2]-x$truth[,1]) , ".\n",sep="") 
    }
  }
  
  if(length(class(x$truth))==1){
    cat("\ngrid of true differences (= treatment - control):\n")
    cat(length(x$truth[!duplicated(x$truth)])," distinct values from ", min(x$truth), " to ", max(x$truth) , ".\n",sep="")
    
  }
  
  invisible(NULL)
}





plot.gsbMainOut <- function(x,
                            what=c("all", "cumulative all", "both", "cumulative both","sample size", "success", "futility", "success or futility", "indeterminate", "cumulative success", "cumulative futility", "cumulative success or futility", "cumulative indeterminate", "boundary", "std.boundary", "delta.grid", "patients"),
                            range.delta="default",
                            stages="default",
                            delta.grid=TRUE,
                            color=TRUE,
                            smooth=100,
                            contour=TRUE,
                            export=FALSE,
                            path=getwd(),
                            sliced=FALSE,
                            range.control="default",
                            ...){


  ## global lattice settings --------------
  fontsize <- trellis.par.get("fontsize")
  fontsize$text <- 11
  trellis.par.set("fontsize", fontsize)
  # ---------------------------------------

  
  what <- match.arg(what)

  range.delta <- integer2numeric(range.delta)
  smooth <- integer2numeric(smooth)
  stages <- integer2numeric(stages)

  
  if(range.delta[1]!="default")
    if(class(range.delta)!="numeric")
      stop("\"range.delta\" hat to be a numeric.")
  
  if(stages[1]!="default"){
    if(class(stages)!="numeric")
      stop("\"stages\" has to b a numeric.")
    if(x$design$nr.stages < max(stages))
      stop("\"stages\" has to be <= nr of stages.")
  }
  if(class(delta.grid) !="logical")
    stop("\"delta.grid\" has to be of class logical.")

  if(class(color) !="logical")
    stop("\"color\" has to be of class logical.")

  if(class(smooth)!="numeric" | length(smooth)!=1)
    stop("\"smooth\" has a numeric of length 1.")    

  if(class(export) !="logical")
    stop("\"export\" has to be of class logical.")

  if(class(path) !="character")
    stop("\"path\" has to be of class character .")

  if((what=="both" | what=="cumulative both") & class(x$OC)[1]!= "gsbSimCalcDif.result")
    stop("what = \"both\" is only a valid input if method = \"both\".")
  
  if(x$design$nr.stages==1){
    if(what=="cumulative all")
      what <- "all"
    if(what=="cumulative both")
      what <- "both"
    if(what=="cumulative success")
      what <- "success"
    if(what=="cumulative futility")
      what <- "futility"
    if(what=="cumulative success or futility")
      what <- "success or futility"
  }


  ## boundary
  if(what=="boundary" || what=="std.boundary"){
    if(!(class(x$simulation)[2]=="gsbSimCalcDif" || class(x$simulation)[2]=="gsbCalcDif")){
      stop("\n To plot criteria, \"type.update\" has to be \"treatment effect\" and \"method\" has to be \"numerical integration\" or \"both\"\n")
    }
    r <- plot.boundary(x, what = what)
  }

  
  if(sliced){
    if(class(x$OC)[1] != "gsbSimArm.result")
      stop("if \"sliced = TRUE\" the \"type.update\" has to be \"per arm\".")
    if(x$simulation$grid.type != "sliced")
      stop("if \"sliced = TRUE\" the \"grid.type\" has to be \"sliced\".")
    
    
    r <- plot.gsbSimArm.sliced(x, what=what, stages = stages, range.control=range.control,range.delta=range.delta)
    
    if(export){
      png(filename=paste(what,".png",sep=""),width=700, height=700)
      print(r)
      dev.off()
    }
    return(r)
  }
  
     
  if(what=="patients"){
    f <- x$design$patients
    if(export){
      png(filename=paste(what,".png",sep=""),width=700, height=700)
      plot(f)
      dev.off()
    }
    return(plot(f))
  }

  if(what=="delta.grid"){
    if(class(x$simulation)[2]!="gsbSimArm"){
      return(cat("the argument \"type.update\" of the function \"gsbSimulation\" has to be set to \"per arm\" to use this plot"))
    }
    if(export){
      png(filename=paste(what,".png",sep=""),width=700, height=700)
      plot(x$delta.grid, main="delta.grid")
      dev.off()
      return(NULL)
    }
    if(!export){
      return(plot(x$delta.grid, main="delta.grid"))
    }
  }

  if(!(what=="patients" || what== "delta.grid" || what== "boundary" || what=="std.boundary")){
    if(class(x$OC)[1]=="gsbSimArm.result"){
      r <- plot.gsbSimArm.result(x=x, what=what, range.delta = range.delta, stages=stages,color=color,delta.grid=delta.grid,smooth=smooth, contour=contour)
    }else{
      r <- plot.gsbSimDif.result(x,what=what,stages=stages, range.delta=range.delta )
    }
  }
  
  if(export){
    if(substr(path,nchar(path),nchar(path))!="/")
      path <- paste(path,"/",sep="")
    
    if(!(what=="cumulative all"|| what=="all" || what=="sample size")){
       png(filename=paste(path,what,".png",sep=""), width=800, height=400)
     }else{
       png(filename=paste(path,what,".png",sep=""), width=800, height=800)
     }
    print(r)
      dev.off()
  }
  ## dev.new( "apropriate size" and return invisible)
  
  return(r)
}

stage <- NULL
plot.gsbSimDif.result <- function(x,
                                  what,
                                  range.delta,
                                  stages,...) {


  i <- length(levels(x$OC$stage))
  
  ## defaults

  if(stages[1]=="default"){
    if(what=="sample size"){
      stages <- i
    }
    if(what!="sample size"){
      stages <- 1:i
    }
  }
  if(class(stages)=="numeric"){
    stages <- stages[!duplicated(stages)]
    ## stages <- stages[stages<=i]
    ## stages <- stages[stages>=0]
  }

  if (range.delta[1]=="default"){
    range.delta <- c(min(x$OC$delta),max(x$OC$delta))
  }
  

  panel <- function(...) {
    panel.grid(h = -1, v = -1, lty = "dotted", lwd=2, col = "light grey")
    panel.xyplot(...)
  }

  if(what=="all"){
    d1 <- subset(x$OC,
                 x$OC$type%in%c("success","futility","indeterminate") & as.numeric(x$OC$stage)%in%stages & x$OC$delta >= range.delta[1] & x$OC$delta<= range.delta[2])
    if(length(levels(d1$method))==2){
      d1 <- subset(d1, d1$method=="numerical integration")
    }
    d1$type <- factor(d1$type, levels=c("success","indeterminate","futility"))
    r <- xyplot(value~delta|stage,
                group=d1$type,
                data=d1,
                col=c("green4","darkorange1","red2"),
                type="l",
                lwd=2.5,
                panel=panel,
                layout=c(length(stages),1),
                xlim=range.delta,
                ylim=c(-.03,1.03),
                ylab="probability",
                xlab="delta = treatment - control\n",
                main="Operating Characteristics",
                scales=list(alternating=1),
                key=list(columns=3,space="bottom",
                  lines=list(col=c("green4","darkorange1","red2"),lwd=2.5),
                  text=list(levels(d1$type),col=1), border=FALSE)
                )
  }
  
  
  if(what=="cumulative all"){
    d1 <- subset(x$OC,
                 x$OC$type%in%c("cumulative success","cumulative futility","cumulative indeterminate") & as.numeric(x$OC$stage)%in%stages & x$OC$delta>= range.delta[1] & x$OC$delta<= range.delta[2])
    if(length(levels(d1$method))==2){
      d1 <- subset(d1, d1$method=="numerical integration")
    }
    d1$type <- factor(d1$type, levels=c("cumulative success","cumulative indeterminate","cumulative futility"))
    
    r <- xyplot(value~delta|stage,
                group=d1$type,
                data=d1,
                col=c("green4","darkorange1","red2"),
                type="l",
                lwd=2.5,
                panel=panel,
                layout=c(length(stages),1),      
                xlim=range.delta,
                ylim=c(-.03,1.03),
                ylab="probability",
                xlab="delta = treatment - control\n",
                main="Cumulative Operating Characteristics",
                scales=list(alternating=1),
                key=list(columns=3,space="bottom",
                            lines=list(col=c("green4","darkorange1","red2"),lwd=2.5),
                  text=list(levels(d1$type),col=1), border=FALSE)
                )
 }


  if(what=="both"){
    d1 <- subset(x$OC,
                 x$OC$type%in%c("success","futility","success or futility") & as.numeric(x$OC$stage)%in%stages & x$OC$delta>= range.delta[1] & x$OC$delta<= range.delta[2])
    d1$type <- factor(d1$type, levels=c("success","futility","success or futility"))
    r <- xyplot(value~delta|stage+type,
                group=d1$method,
                data=d1,
                col=c("blue1","red1"),
                type="l",
                xlim=range.delta,
                ylim=c(-.03,1.03),
                lwd=2.5,
                lty=c(1,2),
                panel=panel,
                ylab="probability of stopping",
                xlab="delta = treatment - control\n",
                main="Operating Characteristics",
                 scales=list(alternating=1),
                key=list(columns=2,space="bottom",
                            lines=list(col=c("blue1","red1"),lwd=2.5, lty=c(1,2)),
                  text=list(as.character(d1$method[!duplicated(d1$method)]),col=1), border=FALSE))
  }

  if(what=="cumulative both"){
    d1 <- subset(x$OC,
                 x$OC$type%in%c("cumulative success","cumulative futility","cumulative success or futility") & as.numeric(x$OC$stage)%in%stages & x$OC$delta>= range.delta[1] & x$OC$delta<= range.delta[2])
    d1$type <- factor(d1$type, levels=c("cumulative success","cumulative futility","cumulative success or futility"))
    r <- xyplot(value~delta|stage+type,
                group=d1$method,
                data=d1,
                col=c("blue1","red1"),
                type="l",
                xlim=range.delta,
                ylim=c(-.03,1.03),
                lwd=2.5,
                lty=c(1,2),
                panel=panel,
                ylab="cumulative probability of stopping",
                xlab="delta = treatment - control\n",
                main="Operating Characteristics",
                scales=list(alternating=1),
                key=list(columns=2,space="bottom",
                            lines=list(col=c("blue1","red1"),lwd=2.5, lty=c(1,2)),
                  text=list(as.character(d1$method[!duplicated(d1$method)]),col=1)))
  }

  if(what=="sample size"){
    d1 <- subset(x$OC,
                 x$OC$type%in%c("sample size") & as.numeric(x$OC$stage)%in%stages & x$OC$delta>= range.delta[1] & x$OC$delta<= range.delta[2])

    if(length(levels(d1$method))==2){
      d1 <- subset(d1, d1$method=="numerical integration")
    }

    if(length(stages)==1)
      lty <- 1
    else
      lty <- rev(1:length(unique(d1$stage)))
    
    
    r <- xyplot(value~delta,
                group=d1$stage,
                data=d1,
                type="l",
                xlim=range.delta,
                ylim=c(0,max(d1$value)+4),
                lwd=2.5,
                col=1,
                lty=lty,
                panel=panel,
                ylab="Patients",
                xlab="delta = treatment - control\n",
                main="Expected Sample Size",
                scales=list(alternating=1),
                key=list(columns=ifelse(length(stages)>4,4,length(stages)), rows= ceiling(length(stages)/4),
                  col=1, space="bottom",
                  lines=list(lwd=2.5, lty=lty),
                  text=list(as.character(d1$stage[!duplicated(d1$stage)]),col=1), border=FALSE))


  }

  if(!(what=="all" || what=="cumulative all" || what=="sample size" || what== "both" || what == "cumulative both")){
    d1 <- subset(x$OC,
                 x$OC$type== what &  as.numeric(x$OC$stage)%in%stages & x$OC$delta >= range.delta[1] & x$OC$delta <= range.delta[2])
    if(length(levels(d1$method))==2){
      d1 <- subset(d1, d1$method=="numerical integration")
    }

    if(length(stages)==1)
      lty <- 1
    else
      lty <- rev(1:length(unique(d1$stage)))
    
    if(what == "success" | what == "cumulative success")
      col <- "green4"
    if(what == "indeterminate" | what == "cumulative indeterminate")
      col <-  "darkorange1"
    if(what == "futility" | what == "cumulative futility")
      col <-  "red2"
    if(what == "success or futility" | what == "cumulative success or futility")
      col <-  1
    
    r <- xyplot(value~delta,
                group=stage, 
                data=d1,
                type="l",
                col=col,
                xlim=range.delta,
                ylim=c(-.03,1.03),
                lwd=2.5,
                lty=lty,
                panel=panel,
                ylab="probability of stopping",
                xlab="delta = treatment - control\n",
                scales=list(alternating=1),
                main=paste("Operating Characteristics: \"",what,"\"",sep=""),
                key=list(columns=ifelse(length(stages)>4,4,length(stages)), rows= ceiling(length(stages)/4),
                  col=1, space="bottom",
                  lines=list(lwd=2.5, lty=lty, col=col),
                  text=list(as.character(d1$stage[!duplicated(d1$stage)]),col=1), border=FALSE))
  }
  return(r)
}

tab <- function(x,
                what = c("all", "cumulative all", "success", "futility",
                  "indeterminate", "success or futility", 
                  "cumulative success", "cumulative futility", "cumulative indeterminate",
                  "cumulative success or futility", "sample size"),
                atDelta = "default",
                wide=FALSE,
                digits = 3,
                export = FALSE, 
                sep = ",",
                path = getwd()) 
{
  
  ## control / prepare arguments
  if (class(x)[1] != "gsbMainOut") 
    stop("\"x\" has to be of class gsbMainOut (output object of function 'gsb()').")
  OC <- x$OC
  i <- length(levels(OC$stage))
  grid.type <- x$simulation$grid.type
  type.update <- x$simulation$type.update
  
  what <- match.arg(what)
  if (i == 1){
    if(what == "cumulative success")
      what <- "success"
    if(what == "cumulative futility")
      what <- "futility"
    if(what == "cumulative indeterminate")
      what <- "indeterminate"
    if(what == "cumulative success or futility")
      what <- "success or futility"
    if(what == "cumulative all")
      what <- "all"
  }
  
  atDelta <- integer2numeric(atDelta)
  if (atDelta[1] != "default" & class(atDelta)[1] != "numeric") 
    stop("\"atDelta\" has to be of class numeric.")
  
  digits <- integer2numeric(digits)
  if (class(digits) != "numeric" | length(digits) != 1 | round(digits, 0) != digits) 
    stop("\"digits\" has to be an integer value of class numeric or intger with length 1.")
  
 
  if (class(export) != "logical") 
    stop("\"export\" has to be of class logical.")
  
  if (class(sep) != "character") 
    stop("\"sep\" has to be of class character.")
  
  if (class(path) != "character") 
    stop("\"path\" has to be of class character.")
  
  
  ## check if valid imput combinations
  
  if (type.update == "treatment effect") 
    result <- tabDif(OC=OC, what=what, atDelta= atDelta, digits=digits)
  
  if (type.update == "per arm"){
    if(wide){
      if(grid.type != "table")
        stop("wide format is only availble for grid.type = \"table\".")
      result <- tabArmWide(OC=OC, what=what, digits=digits)
    }
    else
      result <- tabArmLong(OC=OC, what=what, digits=digits)
  }
  
  if (export) {
    if (substr(path, nchar(path), nchar(path)) != "/") 
      path <- paste(path, "/", sep = "")
    write.table(format((result), scientific = FALSE, trim = TRUE), 
                file = paste(path, what, ".csv", sep = ""), qmethod = c("escape"), 
                sep = sep, col.names = NA)
  }
  return(result)
}




tabDif <- function(OC, what, atDelta, digits){
  i <- length(levels(OC$stage))
  OC <- OC[order(OC$delta),]; OC <- OC[order(OC$stage),]; OC <- OC[order(OC$type),]
  delta <- unique(OC$delta)
  if (length(levels(OC$method)) == 2) OC <- subset(OC, OC$method == "numerical integration")
    
  if(what != "all" & what != "cumulative all"){        
    value <- matrix(nrow = length(delta), ncol = i, 
                    data = subset(OC, OC$type == what)$value, byrow = FALSE)
    dimnames(value) <- list(1:length(delta), paste("stage", 1:i))
    
    if (atDelta[1] != "default") {
      value <- matrix(apply(value,2,function(x){approx(delta,x , xout = atDelta, method = "linear")$y}),ncol=i,nrow=length(atDelta))
      dimnames(value) <- list(1:length(atDelta), paste("stage", 1:i))
      delta <- atDelta
    }
    
    result <- cbind(round(delta,digits=digits),myround(value,digits=digits, what=what))
  }
  
  if(what == "all"){        
    success <- matrix(OC[OC$type=="success","value"], ncol=i)
    futility <- matrix(OC[OC$type=="futility","value"], ncol=i)
    ind <- matrix(OC[OC$type=="indeterminate","value"], ncol=i)
    value <- cbind(success, futility, ind)
    dimnames(value) <- list(1:length(delta),
                            c(paste(paste("stage", 1:i, sep=""), ".suc", sep=""),
                              paste(paste("stage", 1:i, sep=""), ".fut", sep=""),
                              paste(paste("stage", 1:i, sep=""), ".ind", sep="")))
    
    
    if (atDelta[1] != "default") {
      value <- matrix(apply(value,2,function(x){approx(delta,x , xout = atDelta, method = "linear")$y}),ncol=i*3,nrow=length(atDelta))
      dimnames(value) <- list(1:length(atDelta),
                              c(paste(paste("stage", 1:i, sep=""), ".suc", sep=""),
                                paste(paste("stage", 1:i, sep=""), ".fut", sep=""),
                                paste(paste("stage", 1:i, sep=""), ".ind", sep="")))
      
      delta <- atDelta
    }
    
    index <- rep(seq(1, 3*i, i), i) + rep(0:(i-1), each=3)

    if(length(delta)==1){
      result <- round(matrix(c(delta, value[,index]),nrow=1),digits=digits)
      dimnames(result) <- list(1,c("delta", dimnames(value)[[2]][index]))
    }
    else
      result <- round(cbind(delta,value[,index]), digits=digits)
  }
  
  if(what == "cumulative all"){        
    success <- matrix(OC[OC$type=="cumulative success","value"], ncol=i)
    futility <- matrix(OC[OC$type=="cumulative futility","value"], ncol=i)
    ind <- matrix(OC[OC$type=="cumulative indeterminate","value"], ncol=i)
    value <- cbind(success, futility, ind)
    dimnames(value) <- list(1:length(delta),
                            c(paste(paste("stage", 1:i, sep=""), ".cum.suc", sep=""),
                              paste(paste("stage", 1:i, sep=""), ".cum.fut", sep=""),
                              paste(paste("stage", 1:i, sep=""), ".cum.ind", sep="")))
   
    if (atDelta[1] != "default") {
      value <- matrix(apply(value,2,function(x){approx(delta,x , xout = atDelta, method = "linear")$y}),ncol=3*i,nrow=length(atDelta))
      dimnames(value) <- list(1:length(atDelta),
                              c(paste(paste("stage", 1:i, sep=""), ".cum.suc", sep=""),
                                paste(paste("stage", 1:i, sep=""), ".cum.fut", sep=""),
                                paste(paste("stage", 1:i, sep=""), ".cum.ind", sep="")))
      delta <- atDelta
    }
    
    index <- rep(seq(1, 3*i, i), i) + rep(0:(i-1), each=3)

    
    if(length(delta)==1){
      result <- round(matrix(c(delta, value[,index]),nrow=1),digits=digits)
      dimnames(result) <- list(1,c("delta", dimnames(value)[[2]][index]))
    }
    else
      result <- round(cbind(delta,value[,index]), digits=digits)
  }

  result
}

tabArmLong <- function(OC, what, digits){
  i <- length(levels(OC$stage))
  OC <- OC[order(OC$delta.control), ]; OC <- OC[order(OC$delta.treatment), ]
  OC <- OC[order(OC$stage), ]; OC <- OC[order(OC$type), ]
  delta.treatment <- OC$delta.treatment[!duplicated(OC$delta.treatment)]
  delta.control <- OC$delta.control[!duplicated(OC$delta.control)]
  delta.extend <- round(as.matrix(OC[OC$stage %in% "stage 1" & OC$type == "success",c("delta.control","delta.treatment","delta")]), digits)
  dimnames(delta.extend) <- list(paste(1:dim(delta.extend)[1]),c("control", "treatment", "delta"))
  
  if(what == "all" ){
    success <- matrix(OC[OC$type == "success","value"], ncol=i)
    dimnames(success) <- list(rep("",nrow(success)),paste("stage",1:i,".suc", sep=""))
    futility <- matrix(OC[OC$type == "futility","value"], ncol=i)
    dimnames(futility) <- list(rep("",nrow(futility)),paste("stage",1:i,".fut", sep=""))
    ind <- matrix(OC[OC$type == "indeterminate","value"], ncol=i)
    dimnames(ind) <- list(rep("",nrow(ind)),paste("stage",1:i,".ind", sep=""))
    value <- round(cbind(success, futility, ind),digits)
    dimnames(value)[[1]] <- 1:dim(value)[1]
    index <- rep(seq(1, 3*i, i), i) + rep(0:(i-1), each=3)
    result <- cbind(delta.extend,value[,index])
    return(result)
  }
  
  if(what == "cumulative all" ){
    success <- matrix(OC[OC$type == "cumulative success","value"], ncol=i)
    dimnames(success) <- list(rep("",nrow(success)),paste("stage",1:i,".suc", sep=""))
    futility <- matrix(OC[OC$type == "cumulative futility","value"], ncol=i)
    dimnames(futility) <- list(rep("",nrow(success)),paste("stage",1:i,".fut", sep=""))
    ind <- matrix(OC[OC$type == "cumulative indeterminate","value"], ncol=i)
    dimnames(ind) <- list(rep("",nrow(success)),paste("stage",1:i,".ind", sep=""))
    value <- round(cbind(success, futility, ind),digits)
    dimnames(value)[[1]] <- 1:dim(value)[1]
    index <- rep(seq(1, 3*i, i), i) + rep(0:(i-1), each=3)
    result <- cbind(delta.extend,value[,index])
    return(result)
  }
 
  value <- myround(matrix(OC[OC$type == what,"value"], ncol=i), what=what, digits=digits)
  dimnames(value) <- list(1:dim(value)[1],paste("stage",1:i, sep=""))
  result <- round(cbind(delta.extend,value),digits=digits)
  return(result)
}


tabArmWide <- function(OC, what, digits){

  if(what %in% c("all","cumulative all"))
    stop("wide forlmat not available for what = \"all\" and \"cumulative all\"")
 
  OC <- OC[order(OC$delta.control), ]
  OC <- OC[order(OC$delta.treatment), ]
  OC <- OC[order(OC$stage), ]
  OC <- OC[order(OC$type), ]
  delta.treatment <- OC$delta.treatment[!duplicated(OC$delta.treatment)]
  delta.control <- OC$delta.control[!duplicated(OC$delta.control)]
  i <- length(levels(OC$stage))
  value <- array(dim = c(length(delta.control), length(delta.treatment), i), data = subset(OC, OC$type == what)$value)
  value <- round(value, digits = digits)
  dimnames(value) <- list(c(paste("contr", round(delta.control, digits))), c(paste("treat", round(delta.treatment, digits))), paste("stage", 1:i))
  result <- value
  result
}



summaryDif <- function(x, atDelta){
  method <- x$simulation$method
  i <- x$design$nr.stages

  
  if(method== "simulation"){
    dat <- data.frame(matrix(nrow=x$design$nr.stages + 1, ncol= 3))
    colnames(dat) <- c("Analysis","N1","N2")
    dat[,1] <- c("Prior",1:(nrow(dat)-1))
    if(x$design$prior.difference[1]=="non-informative"){
      dat[1,2] <- 0
      dat[1,3] <- 0
    }else{
      dat[1,2] <- x$design$prior.difference[2]
      dat[1,3] <- x$design$prior.difference[3]
    }
    dat[-1,2] <- as.numeric(x$design$patients[,1])
    dat[-1,3] <- as.numeric(x$design$patients[,2])
  }
  else{
    dat <- data.frame(matrix(nrow=x$design$nr.stages + 1, ncol= 7))
    colnames(dat) <- c("Analysis","N1","N2","S","F","std.S", "std.F")
    dat[,1] <- c("Prior",1:(nrow(dat)-1))
    if(x$design$prior.difference[1]=="non-informative"){
      dat[1,2] <- 0
      dat[1,3] <- 0
    }else{
      dat[1,2] <- x$design$prior.difference[2]
      dat[1,3] <- x$design$prior.difference[3]
    }
    dat[-1,2] <- x$design$patients[,1]
    dat[-1,3] <- x$design$patients[,2]
    dat[-1,4] <- subset(x$boundary,x$boundary$type=="success / upper bound")$observed_difference
    dat[-1,5] <- subset(x$boundary,x$boundary$type!="success / upper bound")$observed_difference
    dat[-1,6] <- subset(x$boundary,x$boundary$type=="success / upper bound")$std_observed_difference
    dat[-1,7] <- subset(x$boundary,x$boundary$type!="success / upper bound")$std_observed_difference
  }
  print(dat, digits=3, row.names=FALSE)
  
  cat("\nsigma treatment:", x$design$sigma[1], "\tsigma control:", x$design$sigma[2],"\n")
  
  
  ## table success
  if(atDelta[1]=="default"){
    if(i>1){
      val.suc <- as.data.frame(x$design$criteria)
      val.suc <- subset(val.suc, val.suc$type=="success")$value
      val.suc <- val.suc[min(x$simulation$truth) <= val.suc & val.suc <= max(x$simulation$truth)] 
      logic.suc <- sum(!is.na(val.suc))>0
      
      if(logic.suc){
        val.suc <- val.suc[(!duplicated(val.suc) & !is.na(val.suc))]
        
        suc <- tab(x, what="success", atDelta=val.suc, digits=4)
        dim.suc <- dim(suc)
        total <- tab(x, what="cumulative success", atDelta=val.suc, digits=4)[,dim.suc[2]]
        en <- tab(x, what="sample size", atDelta=val.suc, digits=1)[,dim.suc[2]]
        
        table.s <- cbind(suc,total,en)
        dimnames(table.s)[[2]][dim(suc)[2]+2] <- "E{N}"
        dimnames(table.s)[[2]][1] <- "delta"
        cat("\nstopping for success:\n")
        table.s <- as.data.frame(table.s)
        table.s <- table.s[order(table.s$delta),]
        print(table.s, row.names=FALSE)
      }
      
      
      ## table futility
      val.fut <- as.data.frame(x$design$criteria)
      val.fut <- subset(val.fut, val.fut$type=="futility")$value
      val.fut <- val.fut[min(x$simulation$truth) <= val.fut & val.fut <= max(x$simulation$truth)] 
      logic.fut <- (sum(!is.na(val.fut))>0)
      if(logic.fut){
        val.fut <- val.fut[(!duplicated(val.fut) & !is.na(val.fut))]
        
        fut <- tab(x, what="futility", atDelta=val.fut, digits=4)
        dim.fut <- dim(fut)
        total <- tab(x, what="cumulative futility", atDelta=val.fut, digits=4)[,dim.fut[2]]
        table.f <- cbind(fut,total)
        dimnames(table.f)[[2]][1] <- "delta"
        
        cat("\nstopping for futility:\n")
        table.f <- as.data.frame(table.f)
        table.f <- table.f[order(table.f$delta),]
        print(table.f, row.names=FALSE)
      }
    }
    if(i==1){
      val.suc <- as.data.frame(x$design$criteria)
      val.suc <- subset(val.suc, val.suc$type=="success")$value
      val.suc <- val.suc[min(x$simulation$truth) <= val.suc & val.suc <= max(x$simulation$truth)] 
      logic.suc <- sum(!is.na(val.suc))>0
      
      if(logic.suc){
        val.suc <- val.suc[(!duplicated(val.suc) & !is.na(val.suc))]
        
        suc <- tab(x, what="success", atDelta=val.suc, digits=4)
        en <- tab(x, what="sample size", atDelta=val.suc, digits=1)[,2]
        
        table.s <- cbind(suc,en)
        dimnames(table.s)[[2]][3] <- "E{N}"
        dimnames(table.s)[[2]][1] <- "delta"
        
        cat("\nstopping for success:\n")
        table.s <- as.data.frame(table.s)
        table.s <- table.s[order(table.s$delta),]
        print(table.s, row.names=FALSE)
      }
      ## table futility
      val.fut <- as.data.frame(x$design$criteria)
      val.fut <- subset(val.fut, val.fut$type=="futility")$value
      val.fut <- val.fut[min(x$simulation$truth) <= val.fut & val.fut <= max(x$simulation$truth)] 
      logic.fut <- (sum(!is.na(val.fut))>0)
      if(logic.fut){
        val.fut <- val.fut[(!duplicated(val.fut) & !is.na(val.fut))]
        
        table.f <- tab(x, what="futility", atDelta=val.fut, digits=4)
        dimnames(table.f)[[2]][1] <- "delta"
        
        cat("\nstopping for futility:\n")
        table.f <- as.data.frame(table.f)
        table.f <- table.f[order(table.f$delta),]
        print(table.f, row.names=FALSE)
      }
    }
    ## return invisible
    if(logic.suc){
      if(logic.fut){
        l <- list(design=dat, table.success=table.s, table.futility=table.f)
      }
      if(!logic.fut){
        l <- list(design=dat, table.success=table.s)
      }
    }
    if(!logic.suc){
      if(logic.fut){
        l <- list(design=dat, table.futility=table.f)
      }
      if(!logic.fut){
        l <- list(design=dat)
      }
    }
  }
  if(atDelta[1]!="default"){
    if(i>1){
      suc <- tab(x, what="success", atDelta=atDelta, digits=4)
      dim.suc <- dim(suc)
      total <- tab(x, what="cumulative success", atDelta=atDelta, digits=4)[,dim.suc[2]]
      en <- tab(x, what="sample size", atDelta=atDelta, digits=1)[,dim.suc[2]]
      
      table.s <- cbind(suc,total,en)
      dimnames(table.s)[[2]][dim(suc)[2]+2] <- "E{N}"
      dimnames(table.s)[[2]][1] <- "delta"
      cat("\nstopping for success:\n")
      table.s <- as.data.frame(table.s)
      table.s <- table.s[order(table.s$delta),]
      print(table.s, row.names=FALSE)
      
      ## table futility
      fut <- tab(x, what="futility", atDelta=atDelta, digits=4)
      total <- tab(x, what="cumulative futility", atDelta=atDelta, digits=4)[,dim.suc[2]]
      table.f <- cbind(fut,total)
      dimnames(table.f)[[2]][1] <- "delta"
      
      cat("\nstopping for futility:\n")
      table.f <- as.data.frame(table.f)
      table.f <- table.f[order(table.f$delta),]
      print(table.f, row.names=FALSE)
      
    }
    if(i==1){
      suc <- tab(x, what="success", atDelta=atDelta, digits=4)
      en <- tab(x, what="sample size", atDelta=atDelta, digits=1)[,2]
      
      table.s <- cbind(suc,en)
      dimnames(table.s)[[2]][3] <- "E{N}"
      dimnames(table.s)[[2]][1] <- "delta"
      
      cat("\nstopping for success:\n")
      table.s <- as.data.frame(table.s)
      table.s <- table.s[order(table.s$delta),]
      print(table.s, row.names=FALSE)
      
      ## table futility  
      table.f <- tab(x, what="futility", atDelta=atDelta, digits=4)
      dimnames(table.f)[[2]][1] <- "delta"
      
      cat("\nstopping for futility:\n")
      table.f <- as.data.frame(table.f)
      table.f <- table.f[order(table.f$delta),]
      print(table.f, row.names=FALSE)
    }
    l <- list(design=dat, table.success=table.s, table.futility=table.f)
  }
  return(l)
}


summaryArm <- function(x){
        
  dat <- data.frame(matrix(nrow=x$design$nr.stages + 1, ncol= 3))
  colnames(dat) <- c("Analysis","N1","N2")
  dat[,1] <- c("Prior",1:(nrow(dat)-1))
  if(x$design$prior.control[1]=="non-informative"){
    dat[1,2] <- 0
  }else{
    dat[1,2] <- x$design$prior.control[2]
  }
  if(x$design$prior.treatment[1]=="non-informative"){
    dat[1,3] <- 0
  }else{
      dat[1,3] <- x$design$prior.treatment[2]
    }
  
  dat[-1,2] <- x$design$patients[,1]
  dat[-1,3] <- x$design$patients[,2]
  
  print(dat,digits=3, row.names=FALSE)
  cat("\nsigma treatment:", x$design$sigma[1], "\tsigma control:", x$design$sigma[2],"\n\n")
  cat("access the operating characteristics via the data.frame \"OC\" in the output of \"gsb()\"\nor the functions \"tab()\" and \"plot()\".\n")
  l <- list(design=dat)
  
 return(l) 
}


summary.gsbMainOut <- function(object, atDelta="default", ...){

  x <- object
  i <- x$design$nr.stages
  grid.type <- x$simulation$grid.type
  type.update <- x$simulation$type.update
  
  
  atDelta <- integer2numeric(atDelta)

  if(atDelta[1]!="default" & class(atDelta)[1]!="numeric")
    stop("\"atDelta\" has to be of class numeric.") 
  
  cat("\n*** Group Sequential Bayesian Design ***\n\n")

  if (type.update=="treatment effect")
    l <- summaryDif(x=x, atDelta=atDelta)
  else
    l <- summaryArm(x=x)
  
  invisible(l)
}

plot.boundary <- function(x, what){
  
  panel <- function(...) {
    panel.grid(h = -1, v = -1, lty = "dotted", lwd=2, col = "light grey")
    panel.xyplot(...)
  }

  key <- list(columns=2,space="bottom", points=list(col="black",fill=c("green4","red2"),cex=2, lwd=1.5, pch=c(24,21)), text=list(c("success / upper bound","futility / lower bound"),cex=1.1), border=FALSE)

  x$boundary$observed_difference <- round(x$boundary$observed_difference,5)
  x$boundary$std_observed_difference <- round(x$boundary$std_observed_difference,5)
  
  if(what=="boundary"){
    r <- xyplot(observed_difference~ stage,
                main=list("Decision criteria translated to bounds \n for the observed treatment effect",cex=1.5),
                group=x$boundary$type,
                data=x$boundary,
                type="o",
                ylab=list("Observed treatment effect",cex=1.1),
                xlab="",
                lty=2,
                lwd=2.5,
                cex=2,
                pch=c(21,24),
                col="black",
                fill=c("red2","green4"),
                panel=panel,
                key=key,
                scales=list(cex=1.1),
                aspect=.9)
  }
  if(what=="std.boundary"){
   
    r <- xyplot(as.numeric(std_observed_difference) ~ stage,
                group=x$boundary$type,
                data=x$boundary,
                type="o",
                main=list("Decision criteria translated to bounds \n for the standardized observed treatment effect",cex=1.5),
                ylab=list("Standardized observed treatment effect",cex=1.1),
                xlab="",
                lty=2,
                lwd=2.5,
                cex=2,
                pch=c(21,24),
                col="black",
                fill=c("red2","green4"),
                panel=panel,
                key=key,
                scales=list(cex=1.1),
                aspect=.9)
  }
     
 return(r)
}



type <- NULL
plot.gsbSimArm.sliced <- function(x,
                                  what,
                                  range.control="default",
                                  range.delta="default",
                                  stages="default",
                                  ...) {

  
  
  ## defaults
  i <- length(levels(x$OC$stage))
  
  ## defaults
  if(stages[1]=="default"){
    if(what=="sample size"){
      stages <- i
    }
    if(what!="sample size"){
      stages <- 1:i
    }
  }
  if(class(stages)=="numeric"){
    stages <- stages[!duplicated(stages)]
  }

  if (range.control[1]=="default"){
    range.control <- c(min(x$delta.grid[,1]),max(x$delta.grid[,1]))
  }
  else{
    if(length(range.control) != 2 )
      stop("\"range.control\" has to be numeric of length 2.")
  }
  
  if (range.delta[1]=="default"){
    range.delta <- c(min(x$OC$delta),max(x$OC$delta))
  }
  else{
    if(length(range.delta) != 2 )
      stop("\"range.delta\" has to be numeric of length 2.")
  }
  
  ## global lattice settings --------------
  fontsize <- trellis.par.get("fontsize")
  fontsize$text <- 10
  trellis.par.set("fontsize", fontsize)
  # ---------------------------------------

  
  panel <- function(...) {
    panel.grid(h = -1, v = -1, lty = "dotted", lwd=1.5, col = "light grey")
    panel.xyplot(...)
  }
  
  tmp <- x$OC
  tmp <- subset(tmp, tmp$delta.control <= max(range.control) & 
                tmp$delta.control >= min(range.control) &
                tmp$delta <= max(range.delta) &
                tmp$delta >= min(range.delta) &
                as.numeric(tmp$stage)%in%stages)

## -----------------------------------------------------------------------  

  if(what=="all"){
    d1 <- subset(tmp, tmp$type %in% c("success","futility","indeterminate"))
    d1$type <- factor(d1$type, levels=c("success","indeterminate","futility"))

    tt <- unique(d1$delta.control)[order(unique(d1$delta.control))]
    d1$delta.control <- factor(d1$delta.control, levels=tt, labels=paste("control =",tt) )
    range.delta <- c(min(tmp$delta), max(tmp$delta))
    
    r <- xyplot(value~delta|stage + delta.control,
                group=type,
                data=d1,
                col=c("green4","darkorange1","red2"),
                type="l",
                lwd=2.5,
                panel=panel,
                xlim=range.delta,
                ylim=c(-.03,1.03),
                ylab="probability",
                xlab="delta = treatment - control\n",
                main="Operating Characteristics",
                scales=list(alternating=1),
                key=list(columns=3,space="bottom",
                  lines=list(col=c("green4","darkorange1","red2"),lwd=2.5),
                  text=list(levels(d1$type),col=1), border=FALSE))
    
    return(r)
  }
  
  if(what=="cumulative all"){
    d1 <- subset(tmp, tmp$type %in% c("cumulative success","cumulative futility","cumulative indeterminate"))
    d1$type <- factor(d1$type, levels=c("cumulative success","cumulative indeterminate","cumulative futility"))
    tt <- unique(d1$delta.control)[order(unique(d1$delta.control))]
    d1$delta.control <- factor(d1$delta.control, levels=tt, labels=paste("control =",tt) )
    
    range.delta <- c(min(tmp$delta), max(tmp$delta))
    
    r <- xyplot(value~delta|stage + delta.control,
                group=d1$type,
                data=d1,
                col=c("green4","darkorange1","red2"),
                type="l",
                lwd=2.5,
                panel=panel,
                xlim=range.delta,
                ylim=c(-.03,1.03),
                ylab="probability",
                xlab="delta = treatment - control\n",
                main="Operating Characteristics",
                scales=list(alternating=1),
                key=list(columns=3,space="bottom",
                  lines=list(col=c("green4","darkorange1","red2"),lwd=2.5),
                  text=list(levels(d1$type),col=1), border=FALSE))
    return(r)
  }

  if(what=="sample size"){
    d1 <- subset(tmp,tmp$type%in%"sample size")
    
    if(length(stages)==1)
      lty <- 1
    else
      lty <- rev(1:length(unique(d1$stage)))
    
    tt <- unique(d1$delta.control)[order(unique(d1$delta.control))]
    d1$delta.control <- factor(d1$delta.control, levels=tt, labels=paste("control =",tt) )
    
    r <- xyplot(value~delta | delta.control ,
                group=d1$stage,
                data=d1,
                type="l",
                xlim=range.delta,
                ylim=c(0,max(d1$value)+4),
                lwd=2.5,
                col=1,
                lty=lty,
                panel=panel,
                ylab="Patients",
                xlab="delta = treatment - control\n",
                main="Expected Sample Size",
                scales=list(alternating=1),
                key=list(columns=ifelse(length(stages)>4,4,length(stages)), rows= ceiling(length(stages)/4),
                  col=1, space="bottom",
                  lines=list(lwd=2.5, lty=lty),
                  text=list(as.character(d1$stage[!duplicated(d1$stage)]),col=1), border=FALSE))
  
    return(r)
  }

  
  if(!(what=="all" || what=="cumulative all" || what=="sample size")){
    d1 <- subset(tmp, tmp$type%in%what)
    
    d1$type <- factor(d1$type)
    tt <- unique(d1$delta.control)[order(unique(d1$delta.control))]
    d1$delta.control <- factor(d1$delta.control, levels=tt, labels=paste("control =",tt) )
        
    range.delta <- c(min(tmp$delta), max(tmp$delta))

    if(what == "success" | what == "cumulative success")
      col <- "green4"
    if(what == "indeterminate" | what == "cumulative indeterminate")
      col <-  "darkorange1"
    if(what == "futility" | what == "cumulative futility")
      col <-  "red2"
    if(what == "success or futility" | what == "cumulative success or futility")
      col <-  1

    r <- xyplot(value~delta| delta.control,
                group=d1$stage,
                data=d1,
                type="l",
                lwd=2.5,
                col=col,
                lty=rev(1:length(unique(d1$stage))),
                panel=panel,
                xlim=range.delta,
                ylim=c(-.03,1.03),
                ylab="probability",
                xlab="delta = treatment - control\n",
                main=paste("Operating Characteristics: \"",what,"\"",sep=""),
                scales=list(alternating=1),
                key=list(columns=3,space="bottom",
                  lines=list(col=col,lwd=2.5, lty=rev(1:length(unique(d1$stage)))),
                  text=list(as.character(unique(d1$stage)),col=1), border=FALSE))
    return(r)
  }
}

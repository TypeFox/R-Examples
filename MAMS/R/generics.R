MAMSNews <- function() file.show(system.file("NEWS", package="MAMS"))

print.MAMS <- function (x, digits=max(3, getOption("digits") - 4), ...) {

  cat(paste("Design parameters for a ", x$J, " stage trial with ", x$K, " treatments\n\n",sep=""))

  if(!is.na(x$power)){
    res <- matrix(NA,nrow=2,ncol=x$J)
    colnames(res)<-paste("Stage",1:x$J)
    rownames(res) <- c("Cumulative sample size per stage (control):", "Cumulative sample size per stage (active):")

    res[1,] <- ceiling(x$n*x$rMat[1,])
    res[2,] <- ceiling(x$n*x$rMat[2,])

    print(res)
  
    cat(paste("\nMaximum total sample size: ", x$N,"\n\n"))

  }

  res <- matrix(NA,nrow=2,ncol=x$J)
  colnames(res)<-paste("Stage",1:x$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1,] <- round(x$u,digits)
  res[2,] <- round(x$l,digits)
  
  print(res)

}


summary.MAMS<-function(object, digits=max(3, getOption("digits") - 4), ...){

  cat(paste("Design parameters for a ", object$J, " stage trial with ", object$K, " treatments\n\n",sep=""))

  if(!is.null(object$n)){
    res <- matrix(NA,nrow=2,ncol=object$J)
    colnames(res)<-paste("Stage",1:object$J)
    rownames(res) <- c("Cumulative sample size per stage (control):", "Cumulative sample size per stage (active):")

    res[1,] <- ceiling(object$n*object$rMat[1,])
    res[2,] <- ceiling(object$n*object$rMat[2,])

    print(res)
  
    cat(paste("\nMaximum total sample size: ", object$N,"\n\n"))
  }

  res <- matrix(NA,nrow=2,ncol=object$J)
  colnames(res)<-paste("Stage",1:object$J)
  rownames(res) <- c("Upper bound:", "Lower bound:")
  res[1,] <- round(object$u,digits)
  res[2,] <- round(object$l,digits)
  
  print(res)
}

plot.MAMS <- function (x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", ylab="Test statistic", ylim=NULL, type=NULL, las=1, ...) {

  if(is.null(type))type<-"p"
  if(is.null(pch))pch<-1
  if(is.null(col))col<-1
  if(is.null(lty))lty<-2
  if(is.null(las))las<-1
  if(is.null(ylim)){
    r<-range(x$l,x$u)
    ylim <- c(r[1]-diff(r)/6,r[2]+diff(r)/6)
  }

  matplot(1:x$J, cbind(x$l,x$u), type=type, pch=pch, col=col, ylab=ylab, xlab=xlab, ylim=ylim, main=main, axes=FALSE, las=las, ...)
  mtext(1:x$J,side=1,at=1:x$J)
#  axis(side=2)
  axis(side=2,at=seq(-10,10,1), las=las)
  lines(x$u,lty=lty)
  lines(x$l[1:(x$J)],lty=lty)

}


print.MAMS.sim <- function (x, digits=max(3, getOption("digits") - 4), ...) {

  cat(paste("Simulated error rates based on ", x$nsim," simulations\n\n",sep=""))

  res <- matrix(NA,nrow=4,ncol=1)

  res[1,1] <- round(x$typeI,digits)
  res[2,1] <- round(x$power,digits)
  res[3,1] <- round(x$prop.rej,digits)
  res[4,1] <- round(x$exss,digits)

  if(length(x$ptest)==1){
  rownames(res) <- c("Prop. rejecting at least 1 hypothesis:", "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                      paste("Prop. rejecting hypothesis ",x$ptest,":",sep=""),"Expected sample size:")
  }else{
  rownames(res) <- c("Prop. rejecting at least 1 hypothesis:", "Prop. rejecting first hypothesis (Z_1>Z_2,...,Z_K)",
                      paste("Prop. rejecting hypotheses ",paste(as.character(x$ptest),collapse=" or "),":",sep=""),"Expected sample size:")
  }
  colnames(res)<-""
  
  print(res)

}

summary.MAMS.sim<-function(object, digits=max(3, getOption("digits") - 4), ...){

  print(object)

}



print.MAMS.step_down<- function (x, digits=max(3, getOption("digits") - 4), ...) {

    get.hyp <- function(n){ # find the nth intersection hypothesis (positions of 1s in binary n)
        indlength = ceiling(log(n)/log(2)+.0000001)
        ind = rep(0,indlength)
        newn=n
        
        for (h in seq(1,indlength)){
            ind[h] = (newn/(2^(h-1))) %% 2
            newn = newn - ind[h]*2^(h-1)
        }
        seq(1,indlength)[ind==1]
    }
    
    cat(paste("Design parameters for a ", x$J, " stage trial with ", x$K, " treatments\n\n",sep=""))
    res <- t(x$sample_sizes)
    colnames(res)<-paste("Stage",1:x$J)
    rownames(res) <- c("Cumulative sample size  (control):", paste("Cumulative sample size per stage (treatment ", 1:x$K, "):"))

    print(res)
    cat(paste("\nMaximum total sample size: ", sum(x$sample_sizes[x$J,]),"\n\n"))

    for (i in 1:length(x$l)){

        cat(paste("\nIntersection hypothesis H_{", paste(get.hyp(i), collapse = " "), "}:","\n\n"))
        
        res <- matrix(NA,nrow=3,ncol=x$J)
        colnames(res)<-paste("Stage",1:x$J)
        rownames(res) <- c("Conditional error", "Upper boundary", "Lower boundary")
        res[1,] <- x$alpha_star[[i]]
        res[2,] <- x$u[[i]]
        res[3,] <- x$l[[i]]
  
        print(res)

    }
}


summary.MAMS.step_down<-function(object, digits=max(3, getOption("digits") - 4), ...){

  print(object)

}
             

plot.MAMS.step_down <- function (x, col=NULL, pch=NULL, lty=NULL, main=NULL, xlab="Analysis", ylab="Test statistic", ylim=NULL, type=NULL, bty="n", las=1, ...) {

    get.hyp <- function(n){ # find the nth intersection hypothesis (positions of 1s in binary n)
        indlength = ceiling(log(n)/log(2)+.0000001)
        ind = rep(0,indlength)
        newn=n
        
        for (h in seq(1,indlength)){
            ind[h] = (newn/(2^(h-1))) %% 2
            newn = newn - ind[h]*2^(h-1)
        }
        seq(1,indlength)[ind==1]
    }
    
    if(is.null(type))type<-"p"
    if(is.null(bty))bty<-"n"
    if(is.null(pch))pch<-1
    if(is.null(las))las<-1
    if(is.null(col))col<-1:length(x$l)
    if(length(col) != length(x$l)) stop("There must be as many colours as hypotheses.")
    if(is.null(lty))lty<-2
    if(is.null(ylim)){

        l_min <- min(unlist(lapply(x$l, function(a) min(a[(a!=Inf)&(a!=-Inf)]))))
        if (!is.null(x$z_scores)) l_min <- min(l_min, min(unlist(x$z_scores)[unlist(x$z_scores) != -Inf]))
        u_max <- max(unlist(lapply(x$u, function(a) max(a[(a!=Inf)&(a!=-Inf)]))))
        r <- u_max - l_min
        ylim <- c(l_min - r/6, u_max + r/6)
        
    }
    
    matplot(1:x$J, cbind(x$l[[1]], x$u[[1]]), type=type, pch=pch, main=main, col=0, ylab=ylab, xlab=xlab, ylim=ylim, axes=FALSE, las=las, ...)
    mtext(1:x$J,side=1,at=1:x$J)
    #  axis(side=2)
    axis(side=2,at=seq(-10,10,1), las=las)
    lines(x$u[[1]],lty=lty)
    lines(x$l[[1]][1:(x$J)],lty=lty)

    completed_stages <- length(x$z_scores)
    if (completed_stages > 0){
        for (i in 1:completed_stages){
            for (k in 1:x$K){
                points(i, x$z_scores[[i]][k], col = 2 ^ (k - 1), pch = 3)
            }
        }
    }
        

    
    legend_text <- NULL
    
    #if (length(col) < length(x$l)) col <- rep(col, length(x$l))
        
    for (i in 1:length(x$l)){
        legend_text <- c(legend_text, paste("H_{", paste(get.hyp(i), collapse = " "), "}"))
        legend_col <- c(col, i)
        if ((x$alpha_star[[i]][x$J] > 0) && (x$alpha_star[[i]][x$J] < 1)){
            
            matpoints(1:x$J, cbind(x$l[[i]], x$u[[i]]), type=type, pch=pch, col=col[i], ylab=ylab, xlab=xlab, ylim=ylim, axes=FALSE, ...)
           
            lines(x$u[[i]],lty=lty, col = col[i])
            lines(x$l[[i]][1:(x$J)],lty=lty, col=col[i])
        }
        
      
        
    }

    legend("bottomright", legend=legend_text, bty=bty, lty=lty, col=col)
}

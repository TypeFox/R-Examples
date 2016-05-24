plot.ndlClassify <- function(x, values="weights", ...)
{ 
  if(values=="weights")
    plot.ndlWeights(x, ...)
  if(values=="probabilities")
    plot.ndlProbabilities(x, ...)
}

plot.ndlWeights <- function(x, type="density", predictors=NULL, outcomes=NULL, panes="single", lty=NULL, col=NULL, mfrow=NULL,  main=NULL, legend.position="topright", ...)
{ 
  weightMatrix <- x$weightMatrix

  if(!is.null(predictors))  
    { weightMatrix <- weightMatrix[grep(predictors,rownames(weightMatrix)),]
      if(length(weightMatrix)==0)
        stop("Predictors not present in NDL model.")
    }
  if(!is.null(outcomes))
     { if(!all(outcomes %in% colnames(weightMatrix)))
         stop("Outcomes not present in NDL model.")
       n.outcomes <- length(outcomes)
       weightMatrix <- weightMatrix[,which(colnames(weightMatrix) %in% outcomes)]
     }
  if(is.null(predictors) & is.null(outcomes) & panes=="single")
    { weightMatrix <- as.vector(weightMatrix)
      n.outcomes=1
    }
  if(is.null(predictors) & is.null(outcomes) & panes=="multiple")
    { outcomes = colnames(weightMatrix)
      n.outcomes = length(outcomes)
    }
  if(is.null(outcomes))
    { outcomes=colnames(weightMatrix)
      n.outcomes <- length(outcomes)
    }

  if(panes=="single" & is.null(mfrow))
    { par(mfrow=c(1,1))
      mfrow=NULL
      if(type=="hist")
        col=0
    }
  if(panes=="multiple" & is.null(mfrow))
    mfrow=c(ceiling(n.outcomes/sqrt(n.outcomes)), ceiling(sqrt(n.outcomes)))
  if(is.null(lty) & !is.null(mfrow))
    lty=1
  if(is.null(col) & !is.null(mfrow))
    if(type=="hist") 
      col=0
    else
      col=1

  if(is.null(lty))
    lty=seq(1:n.outcomes)
  if(is.null(col))
    col=seq(1:n.outcomes)

  if(length(lty)<n.outcomes)
    lty <- rep(lty,ceiling(n.outcomes/length(lty)))[1:n.outcomes]
  if(length(col)<n.outcomes)
    col <- rep(col,ceiling(n.outcomes/length(col)))[1:n.outcomes]

  if(is.null(outcomes))
    mains=colnames(weightMatrix)
  else
    mains=outcomes

  if(type=="hist")
    { if(is.null(main)) main="Strenghts of outcome-specific predictor weights"
      if(is.null(mfrow))
        hist(weightMatrix, col=col[1], main=main, xlab="Weights", ...)
      else
        { par(mfrow=mfrow)
          for(i in 1:n.outcomes)
               hist(weightMatrix[,outcomes[i]], col=col[i], main=mains[i], xlab="Weights", ...)
        }
    }

  if(type=="density")
    { if(n.outcomes>1)
        w.densities <- apply(weightMatrix,2,density)
      else
        w.densities <- list(density(weightMatrix))
      if(is.null(main)) main="Densities of outcome-specific predictor weights"
      if(is.null(mfrow))
        { x.range <- range(sapply(w.densities, function(z) range(z$x)))
          y.range <- range(sapply(w.densities,function(z) range(z$y)))
          plot(w.densities[[1]], lty=lty[1], col=col[1], xlim=x.range, ylim=y.range, main=main, xlab="Weight", ...)
          if(n.outcomes>1)
            for(i in 2:n.outcomes)
               lines(w.densities[[i]], lty=lty[i], col=col[i], ...)
          if(!is.null(outcomes))
            legend(legend.position, legend=mains, lty=lty, col=col, ...)
        }
      else
        { par(mfrow=mfrow)
          plot(w.densities[[1]], lty=lty[1], col=col[1], main=mains[1], ...)
          if(n.outcomes>1)
            for(i in 2:n.outcomes)
               plot(w.densities[[i]], lty=lty[i], col=col[i], main=mains[i], ...)
        }
    }
}


plot.ndlProbabilities <- function(x, type="density", select="all", panes="single", lty=NULL, col=NULL, pch=NULL, mfrow=NULL, main=NULL, legend.position="topright", ...)
{ activationMatrix <- x$activationMatrix
  outcomes <- colnames(activationMatrix)
  n.outcomes <- length(outcomes)
  n.instances <- nrow(activationMatrix)
  probabilities <- apply(acts2probs(activationMatrix)$p, 1, function(z) sort(z,decreasing=TRUE))

  if(length(select)>1)
    if(all(is.numeric(select)) & min(select)>=1 & max(select)<=n.outcomes)
      select.index=select
    else
      stop(paste(c("Indices out of bounds [1,",n.outcomes,"] in parameter 'select': ",paste(select,collapse=",")),collapse=""))
  else
    { if(select=="max")
        select.index=1
      if(select=="min")
        select.index=n.outcomes
      if(select=="maxmin" | select=="minmax")
        select.index=c(1,n.outcomes)
      if(select=="all")
        select.index=seq(1:n.outcomes)
    }

  ranks <- as.character(select.index)
  ranks[which(ranks=="1")]="Maximum"
  ranks[which(ranks==as.character(n.outcomes))]="Minimum"
  mains <- sapply(ranks, function(z) paste(c("Rank: ",z),collapse=""))

  if(panes=="single" & is.null(mfrow))
    { par(mfrow=c(1,1))
      mfrow=NULL
    }
  if(panes=="multiple" & is.null(mfrow))
    mfrow=c(ceiling(length(select.index)/ceiling(sqrt(length(select.index)))), ceiling(sqrt(length(select.index))))
  if(is.null(lty) & !is.null(mfrow))
    lty=1
  if(is.null(col) & !is.null(mfrow))
    col=1
  if(is.null(pch) & !is.null(mfrow))
    pch=1

  if(is.null(lty))
    lty=seq(1:length(select.index))
  if(is.null(col))
    col=seq(1:length(select.index))
  if(is.null(pch))
    pch=seq(1:length(select.index))

  if(length(lty)<length(select.index))
    lty <- rep(lty,ceiling(length(select.index)/length(lty)))[1:length(select.index)]
  if(length(col)<length(select.index))
    col <- rep(col,ceiling(length(select.index)/length(col)))[1:length(select.index)]
  if(length(pch)<length(select.index))
    pch <- rep(pch,ceiling(length(select.index)/length(pch)))[1:length(select.index)]

  if(type=="density")
    { p.densities <- apply(probabilities,1,density)[select.index]
      if(is.null(main)) main="Densities of instance-wise ranked probabilities"
      if(is.null(mfrow))
        { x.range <- range(sapply(p.densities, function(z) range(z$x)))
          y.range <- range(sapply(p.densities,function(z) range(z$y)))
          plot(p.densities[[1]], lty=lty[1], col=col[1], xlim=x.range, ylim=y.range, main=main, xlab=paste("Probability (N = ",n.instances,")",sep=""), ...)
          if(length(select.index)>1)
            for(i in 2:length(select.index))
               lines(p.densities[[i]], lty=lty[i], col=col[i], ...)
          legend(legend.position, legend=mains, lty=lty, col=col, ...)
        }
      else
        { par(mfrow=mfrow)
          plot(p.densities[[1]], lty=lty[1], col=col[1], main=mains[1], ...)
          if(length(select.index)>1)
            for(i in 2:length(select.index))
               plot(p.densities[[i]], lty=lty[i], col=col[i], main=mains[i], ...)
        }
    }

  if(type=="sort")
    { p.values <- matrix(apply(probabilities,1,sort)[,select.index],,length(select.index))
      if(is.null(main)) main="Sorted values of instance-wise ranked probabilities"
      if(is.null(mfrow))
        { y.range <- range(p.values)
          plot(p.values[,1], pch=pch[1], col=col[1], main=main, xlab=paste("N = ",n.instances,sep=""), ylab="Probability", ylim=y.range, ...)
          if(length(select.index)>1)
            for(i in 2:length(select.index))
               points(p.values[,i], pch=pch[i], col=col[i], ...)
          legend(legend.position, legend=mains, pch=pch, col=col, ...)
        }
      else
        { par(mfrow=mfrow)
          plot(p.values[,1], pch=pch[1], col=col[1], main=mains[1], xlab=paste("N = ",n.instances,sep=""), ylab="Probability", ...)
          if(length(select.index)>1)
            for(i in 2:length(select.index))
               plot(p.values[,i], pch=pch[i], col=col[i], main=mains[i], xlab=paste("N = ",n.instances,sep=""), ylab="Probability", ...)
        }
    }

}

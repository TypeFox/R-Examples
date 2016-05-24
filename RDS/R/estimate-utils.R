#TODO what is this function???
#Is it the empirical likelihood standard error? What assumptions?
#was EL_se
#  @export 
EL.se<-function(weights,outcome,N=NULL, use.second.order=TRUE, homophily=1){
  nas <- is.na(weights)|is.na(outcome)
  nas <- nas | is.infinite(weights)|is.infinite(outcome)
  if(use.second.order){
      wij <- outer(weights[!nas],weights[!nas],"*")
      wij <- wij * (1+(homophily-1)*outer(outcome[!nas],outcome[!nas],"!="))
#
      wij[wij > 10*stats::median(wij)] <- 10*stats::median(wij)
#     Compute the EL standard error (if available)
      iLj <- row(wij)>col(wij)
      oij <- 0.5*outer(outcome[!nas],outcome[!nas],"+")
      estij <- sum(oij*wij*iLj)/sum(wij*iLj)
      Gi <- apply((wij*wij*(oij-estij))*iLj,1,sum)
      wvi <- apply((wij*wij)*iLj,1,sum)
      varest <- 4*sum(Gi^2)/(sum(wvi)^2)
  }else{
      wi <- weights[!nas]
      oi <- outcome[!nas]
#
      wi[wi > 10*stats::median(wi)] <- 10*stats::median(wi)
#     Compute the EL standard error (if available)
      esti  <- sum(oi*wi)/sum(wi)
      varest <- sum(((oi-esti)*wi*wi)^2)/(sum(wi*wi)^2)
  }
  if(is.null(N) | !is.numeric(N)){
    sqrt( varest )
  }else{
    sqrt( (N-length(outcome))*varest/(N-1) )
  }
}

#  @export 
EL.se.new<-function(weights,outcome,N=NULL, use.second.order=FALSE, homophily=1){
  nas <- is.na(weights)|is.na(outcome)
  nas <- nas | is.infinite(weights)|is.infinite(outcome)
  if(use.second.order){
      wi <- weights[!nas]
      wi[wi > 10*median(wi)] <- 10*median(wi)
      wij <- outer(weights[!nas],weights[!nas],"*")
#     wij <- wij * (1+(homophily-1)*outer(outcome[!nas],outcome[!nas],"!="))
      wij[wij > 10*median(wij)] <- 10*median(wij)
#     Compute the EL standard error (if available)
      iLj <-  row(wij) >col(wij)
      iLjn <- row(wij)!=col(wij)
      oij <- 0.5*outer(outcome[!nas],outcome[!nas],"+")
      diag(oij) <- 0
      estij <- sum(oij*wij*iLj)/sum(wij*iLj)
      #
      qij <- 2*(oij-estij)*wij*wij*iLjn
      q <- sum(wij*wij*iLjn)
      qi <- apply(qij,2,sum)
      cross.mat=1-wij/outer(wi,wi,"*")
      sgn=sign(sum((cross.mat/wij)[iLjn]))
      varest <- (sum(qij*qij)+sgn*sum(cross.mat*outer(qi,qi,"*")*iLjn))/(q*q)
  }else{
      wij <- outer(weights[!nas],weights[!nas],"*")
#     wij <- wij * (1+(homophily-1)*outer(outcome[!nas],outcome[!nas],"!="))
      wij[wij > 10*stats::median(wij)] <- 10*stats::median(wij)
#     Compute the EL standard error (if available)
      iLj <-  row(wij) >col(wij)
      iLjn <- row(wij)!=col(wij)
      oij <- 0.5*outer(outcome[!nas],outcome[!nas],"+")
      diag(oij) <- 0
      estij <- sum(oij*wij*iLj)/sum(wij*iLj)
      qij <- 2*(oij-estij)*wij*wij*iLjn
      q <- sum(wij*wij*iLjn)
      varest <- sum(qij*qij)/(q*q)
  }
# if(is.null(N) | !is.numeric(N)){
    sqrt( varest )
# }else{
#   sqrt( (N-length(outcome))*varest/(N-1) )
# }
}




EL.est<-function(weights,outcome,N=NULL, use.second.order=TRUE, homophily=1.35){
  nas <- is.na(weights)|is.na(outcome)
  nas <- nas | is.infinite(weights)|is.infinite(outcome)
  if(use.second.order){
      wij <- outer(weights[!nas],weights[!nas],"*")
      wij <- wij / (1+(homophily-1)*outer(outcome[!nas],outcome[!nas],"=="))
#
#     wij[wij > 10*stats::median(wij)] <- 10*stats::median(wij)
#     Compute the EL standard error (if available)
      iLj <- row(wij)>col(wij)
      oij <- 0.5*outer(outcome[!nas],outcome[!nas],"+")
      estij <- sum(oij*wij*iLj)/sum(wij*iLj)
      Gi <- apply((wij*wij*(oij-estij))*iLj,1,sum)
      wvi <- apply((wij*wij)*iLj,1,sum)
      varest <- 4*sum(Gi^2)/(sum(wvi)^2)
      est <- estij
  }else{
      wi <- weights[!nas]
      oi <- outcome[!nas]
#
      wi[wi > 10*stats::median(wi)] <- 10*stats::median(wi)
#     Compute the EL standard error (if available)
      esti  <- sum(oi*wi)/sum(wi)
      varest <- sum(((oi-esti)*wi*wi)^2)/(sum(wi*wi)^2)
      est <- esti
  }
  if(is.null(N) | !is.numeric(N)){
    attr(est,"EL.se") <- sqrt( varest )
  }else{
    attr(est,"EL.se") <- sqrt( (N-length(outcome))*varest/(N-1) )
  }
  est
}

#' Get Horvitz-Thompson estimator assuming inclusion probability proportional
#' to the inverse of network.var (i.e. degree).
#' @param rds.data An rds.data.from
#' @param group.variable The grouping variable.
#' @param network.var The network.size variable.
#' @export
get.h.hat <- function(rds.data,
		group.variable,
		network.var=attr(rds.data,"network.size")){
	
	if(is(rds.data,"rds.data.frame")){
		stopifnot(group.variable %in% names(rds.data))
		stopifnot(network.var %in% names(rds.data))
	}
	else{
		stop("rds.data must be of type rds.data.frame")
	}
	
	harmonic.mean <- function(x){
		x <- x[!is.na(x)]
		length(x)/sum(1/x) 	
	}
	
	
	degrees <- rds.data[[network.var]]
	temp <- function(x){  
		deg <- degrees[as.character(rds.data[[group.variable]]) == x]
		return(harmonic.mean(deg))
	}	
	
	group.names <- sort(unique(as.character(rds.data[[group.variable]])))
	group.names <- group.names[!is.na(group.names)]
	
	h.hat <- sapply(group.names,temp)
	names(h.hat) <- group.names
	return(h.hat)    	
}


# Compute the horvitz thompson estimate
HT.estimate<-function(weights,outcome){
	nas <- is.na(weights)|is.na(outcome)
	if(is.factor(outcome)){
		num <- rep(0,length=length(levels(outcome[!nas])))
		a <- tapply(weights[!nas],as.numeric(outcome[!nas]),sum)
		num[as.numeric(names(a))] <- a
	}else{
		num<-sum(outcome[!nas]*weights[!nas])
	}
	den<-sum(weights[!nas])
	num/den
}




#' Counts the number or recruiter->recruitee transitions between different levels
#' of the grouping variable.
#' @param rds.data An rds.data.frame
#' @param group.variable The name of a categorical variable in rds.data
#' @export
#' @examples 
#' data(faux)
#' count.transitions(faux,"X")
count.transitions <- function(rds.data, group.variable)
{
	stopifnot(is.rds.data.frame(rds.data))
	
	#extract data
	id <- get.id(rds.data)
	recruiter.id <- get.rid(rds.data)
	grp <- as.factor(rds.data[[group.variable]])
#       Make sure the factor labels are alphabetic!
	grp=factor(grp,levels=levels(grp)[order(levels(grp))])
	
	#count transitions
	ri <- match(recruiter.id,id)
	rgrp <- grp[ri]
	res <- table(rgrp,grp)
	for(i in 1:nrow(res)){
	  if(sum(res[i,])==0){res[i,] <- 0.1}
	}
	class(res) <- "matrix"
	res
}



# This function assumes that transition counts can be turned into MLEs.
# This will not be possible if any of the rows in the transition counts
# matrix have sum zero.

# IF 7/12/13: not sure what the point of this function is, I simplified the
# code, but why not just use prop.table

#' calculates the mle. i.e. the row proportions of the transition matrix
#' @param transition.counts a matrix or table of transition counts
#' @details depreicated. just use prop.table(transition.counts,1)
transition.counts.to.Markov.mle <- function(transition.counts){
	#return(prop.table(tij,1))
	tij <- transition.counts + .Machine$double.eps
	ti <- margin.table(tij, 1)
	stopifnot(ti>0)
	sweep(tij, 1, ti, "/", check.margin = FALSE)
}

#' Markov chain statistionary distribution
#' @param mle The transition probabilities
#' @export
#' @return A vector of proportions representing the proportion in each group at the
#' stationary distribution of the Markov chain.
get.stationary.distribution <- function(mle){
	stat.dist <- NULL	
	eigen.analysis <- eigen(t(mle))	
	
	# If -1 is an eigenvalue then we have a periodic chain ... which is bad.
	values.close.to.minus.one <- abs(eigen(t(mle))$values + 1.0) < 
			sqrt(.Machine$double.eps)
	if(sum(values.close.to.minus.one) == 1){
		cat(" Markov chain apears to be periodic.\n")
		stat.dist <- rep(NA,nrow(mle))
		attr(stat.dist,'status') <- "Periodic Markov MLE"
		return(stat.dist)
	}
	
	# The stationary distribution will be given by
	# the eigenvector with eigenvalue equal to one.    
	values.close.to.one <- abs(eigen(t(mle))$values - 1.0) <
			sqrt(.Machine$double.eps)
	
	# If there is exactly one of these, then we can take it as eigenvalue corresponding to the
	# stationary distribution.
	if(sum(values.close.to.one) == 1){
		stationary.dist.eigen.vector <- eigen.analysis$vectors[,which(values.close.to.one)]
		stationary.dist.eigen.vector <- Re(stationary.dist.eigen.vector)
		stat.dist <- stationary.dist.eigen.vector/sum(stationary.dist.eigen.vector)
		attr(stat.dist,'status') <- "MLE gives unique stationary distribution."
		names(stat.dist) <- colnames(mle)	
	}
	return(stat.dist)	
}






#' Estimate LR exceedance probabilities (with importance sampling)
#' 
#' @param t numeric (vector), threshold
#' @param dists list of per-locus probability distributions of a likelihood ratio
#' @param N integer, number of samples
#' @param dists.sample if dists.sample is not equal to dists, then importance sampling is applied. Sampling is done according to dists.sample, while exceedance probabilities are estimated for dists.
#' @return numeric (vector) with estimated probabilities
#' @details For a combined likelihood ratio \deqn{LR=LR_1 LR_2 \times LR_m,}
#' define \eqn{q_{t|H}} as the probability that the LR exceeds \eqn{t} under hypothesis \eqn{H}, i.e.:
#' \deqn{q_{t|H} := P(LR>t|H).}
#' The hypothesis \eqn{H} can be \eqn{H_p}, \eqn{H_d} or even another hypothesis. The current function estimates \eqn{q_{t|H}} by taking \eqn{N} samples from the distributions specified by the \code{dists} parameter and computing the empirical fraction of the product of the samples that exceeds \eqn{t}.
#' 
#' Importance sampling can be used by supplying different distributions to sample from and to estimate the exceedance probabilities for. For instance, the exceedance probability for \eqn{H_d} can be estimated by sampling from \eqn{H_p} and an appropriate weighting of the samples. See the paper and examples for details.
#' 
#' @examples
#'
#' data(freqsNLngm)
#'
#'# dist of PI for true parent/offspring pairs
#'hp <- ki.dist(hyp.1="PO",hyp.2="UN",hyp.true="PO",freqs.ki=freqsNLngm)
#'
#'# dist of PI for unrelated pairs
#'hd <- ki.dist(hyp.1="PO",hyp.2="UN",hyp.true="UN",freqs.ki=freqsNLngm)
#'
#'set.seed(100)
#'
#'# estimate P(PI>1e6) for true PO
#'sim.q(t=1e6,dists=hp)
#'
#'# estimate P(PI>1e6) for unrelated pairs
#'sim.q(t=1e6,dists=hd) # small probability, so no samples exceed t=1e6
#'
#'# importance sampling can estimate the small probability reliably
#'# by sampling from H_p and weighting the samples appropriately
#'sim.q(t=1e6,dists=hd,dists.sample=hp)
#'
sim.q <- function(t,dists,N=1e5,dists.sample=dists){
  if (any(sapply(dists,function(d0) anyDuplicated(d0$x)))) stop(
    "dists contains duplicated events, apply dist.unique.events() first")
  if (any(sapply(dists.sample,function(d0) anyDuplicated(d0$x)))) stop(
    "dists.sample contains duplicated events, apply dist.unique.events() first")
  
  if (identical(dists,dists.sample)){ # simple monte carlo
    samples <- lapply(dists,function(d0)
      sample(x=d0$x,size=N,replace=TRUE,prob=d0$fx)
    )  
    # take product of samples
    samples <- Reduce("*",samples)
    # simple monte carlo estimate
    return(  sapply(t,function(t0) mean(samples>t0)))    
  }else{ # importance sampling
    # sample from dists.sample, but compute exc. prob for dists
    
    # align the outcomes of the two distributions
    for(i in seq_along(dists)){
      ds <- dists.sample[[i]]
      d <- dists[[i]]
      ds.pos <- ds$fx>0
      map.ds.to.d <- match(ds$x[ds.pos],d$x,NA)
      
      # check if all outcomes of d are in ds
      map.d.to.ds <- match(d$x[d$fx>0],ds$x,NA)
#       if (any(is.na(map.d.to.ds))||any(ds$fx[map.d.to.ds]<=0)) warning("Not every outcome with positive probability under dists also has positive probability under dists.sample.")
      
      dists.sample[[i]] <- list(x=ds$x[ds.pos],fx=ds$fx[ds.pos])
      dists[[i]] <- list(x=d$x[map.ds.to.d],fx=d$fx[map.ds.to.d])      
      if (any(is.na(dists[[i]]$fx))) stop("For importance sampling, every outcome with positive probability under dists.sample should have positive probability under dists.")      
      if (any(dists[[i]]$fx==0)) stop("For importance sampling, every outcome with positive probability under dists.sample should have positive probability under dists.")      
    }
    
    if (!identical(lapply(dists,function(y) y$x),
                   lapply(dists.sample,function(y) y$x))){
      stop("For importance sampling, dists and dists.sample should have same support")
    }
    
    inds <- lapply(dists.sample,function(d0)
      sample.int(n=length(d0$x),size=N,replace=TRUE,prob=d0$fx)
    ) 
    samples <- lapply(seq_along(inds),function(i)
      dists.sample[[i]][["x"]][inds[[i]]])
    samples <- Reduce("*",samples)
    weights <- lapply(seq_along(inds),function(i)
      dists.sample[[i]][["fx"]][inds[[i]]]/dists[[i]][["fx"]][inds[[i]]])
    weights <- Reduce("*",weights)
    # importance sampling estimate
    return(sapply(t, function(t0) mean(as.numeric(samples>t0)/weights)))
  }
}
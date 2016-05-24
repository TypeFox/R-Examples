findRecruiter <- function(active.coupons, coupon){
  #     coupon <- 71114
  #      active.coupons <- list()
  #      active.coupons[[1]] <- c(TRUE, TRUE, TRUE)
  #      names(active.coupons[[1]]) <- c(61004,61005,61006)
  #      active.coupons[[2]] <- c(TRUE, TRUE, TRUE)
  #      names(active.coupons[[2]]) <- c(71114,71115,71116)
  
  result <- NA
  recruiter.ind <- sapply(active.coupons, function(x) coupon %in% names(x))
  if(any(recruiter.ind)){
    # Report coupon number:
    recruiter.num <- which(recruiter.ind)
    coupon.ind <- names(active.coupons[[recruiter.num]]) %in% coupon
    coupon.num <- which(coupon.ind)
    
    result <- list(
      recruiter=recruiter.num,
      coupon=coupon.num
    ) 
  }
  return(result)
}
## Testing:
# coupon <- 71114s
# active.coupons <- list()
# active.coupons[[1]] <- c(TRUE, TRUE, TRUE)
# names(active.coupons[[1]]) <- c(61004,61005,61006)
# active.coupons[[2]] <- c(TRUE, TRUE, TRUE)
# names(active.coupons[[2]]) <- c(71114,71115,71116)
# findRecruiter(active.coupons, coupon)






makeSnowBall <- function(rds.sample, seeds){
  coupon.inds <- grepl('coup[0-9]*', names(rds.sample))
#   sample.length <- ncol(rds.sample)
  
  I.t <- rep(NA, nrow(rds.sample))
  I.t[1] <- seeds
  degree.in <- rep(0, nrow(rds.sample))
degree.in[1] <- rds.sample[1, 'NS1']
  degree.out <- rep(0, nrow(rds.sample))
  active.coupons <- list()
  for(period in 2:length(I.t)){
    ### Sketch:
    ## if recruiter in data: 
    # remove incoming coupn
    # update I.t if coupons depleted
    ## if recruiter not in data:
    # update I.t(?)
    ## if coupons handed:
    # add distributed coupons
    # update I.t 
    
    # period <- 2
    # period <- period+1
    I.t.running <- I.t[period-1] 
    
    reference.coupon <- rds.sample[period,'refCoupNum']
    # reference.coupon <- 71114
    recruiter <- findRecruiter(active.coupons, reference.coupon )
    
    if(length(recruiter)>1){
      # revoke coupon
      active.coupons[[recruiter$recruiter]][recruiter$coupon] <- FALSE
      
      # if coupons depleted:
      if(all(!active.coupons[[recruiter$recruiter]])) {
        I.t.running <- I.t.running-1
        degree.out[period] <- attributes(active.coupons[[recruiter$recruiter]])$degree
      }
    }
    
    
    # Increase snowball if coupons handed:
    new.coupons <- unlist(rds.sample[period, coupon.inds])
    if(isTRUE(any(as.logical(new.coupons)))){
      new.guy <- rep(TRUE, sum(!is.na(new.coupons)))
      attributes(new.guy) <-list(
        names=new.coupons, 
        degree= rds.sample[period,'NS1'])
      
      active.coupons <- c(active.coupons, list(new.guy))
      
      I.t.running <- I.t.running+1
      degree.in[period] <- rds.sample[period, 'NS1']
    }
    
    I.t[period] <- I.t.running
    
  }
  
  return(list(I.t=I.t,
              degree.in=degree.in,
              degree.out=degree.out))
}
# ## Testing:
# test.snowball <- makeSnowBall(rds.sample)
# str(test.snowball)
# table(rds.sample$NS1)
# table(test.snowball$degree.in)
# table(test.snowball$degree.out)

rdsObjectConstructor <- function(rds.sample=NULL,
                                 I.t=NULL,
                                 degree.in=NULL,
                                 degree.out=NULL,
                                 original.ordering=NULL,
                                 estimates=NULL){
  result <- list(rds.sample=rds.sample,
                 I.t=I.t,
                 degree.in=degree.in,
                 degree.out=degree.out,
                 original.ordering=original.ordering,
                 estimates=estimates)
  class(result) <- 'rds-object'
  return(result)
}
## Testing
# rdsObjectConstructor()

initializeRdsObject <- function(rds.sample, bin=1L, seeds=1L){
  ## Verification:
  if(any(table(rds.sample[,'interviewDt'])>1)) warning('Non unique interview times. Ignoring and proceeding...')
  
  ## Initialization:
  
  ord <- order(rds.sample[,'interviewDt'])
  rds.sample <-rds.sample[ord,] 
  if(bin>1){
    cuts <- cut(rds.sample$NS1, breaks = seq(from=0, to = (max(rds.sample$NS1)+bin), by = bin))
    rds.sample$NS1 <- as.numeric(cuts)
  }
    
  I.t <- makeSnowBall(rds.sample, seeds=seeds)
  result <- rdsObjectConstructor(rds.sample=rds.sample,
                                 I.t = I.t$I.t,
                                 degree.in = I.t$degree.in,
                                 degree.out = I.t$degree.out,
                                 original.ordering = ord)
  
  
  
  return(result)
}
## Testing:
# rds.object <- initializeRdsObject(rds.sample)
# ls.str(rds.object)







formatArrivalTimes <- function(arrival.times){
    
  ## Logic:
  # convert to numeric
  # set origin to sampling start
  # convert to informative scale
  arrival.times.numeric <- as.numeric(arrival.times)
  arrival.times.origin <- arrival.times.numeric - min(arrival.times.numeric, na.rm = TRUE)
  for(k in 1:10){
    if(max(arrival.times.origin %% 10^k)>0) break
    }
  arrival.times.clean <- arrival.times.origin /10^(k-1)
  return(arrival.times.clean)
}
## Testing:
# chords:::formatArrivalTimes(rds.object$rds.sample$interviewDt)



smoothDegrees <- function(degree.counts){
  lambda <- median(degree.counts)
  xs <- as.integer(names(degree.counts))
  ns <- sum(degree.counts)
  #   poisons <- ns * dpois(xs ,lambda)
  #   plot(degree.counts)
  dispersion <- mad(degree.counts)
  poisons <- ns * dnbinom(xs, size=dispersion, mu = lambda)
  #   points(poisons)
  .attribs <- attributes(degree.counts)
  degree.counts <- poisons
  attributes(degree.counts) <- .attribs
  return(degree.counts)
}



imputeEstimates <- function(Nk.estimates, n.k.counts, convergence){
  conv.inds <- as.logical(convergence==0)
  impute.inds <- as.logical(convergence==1L)
  
  if( sum(conv.inds, na.rm=TRUE)>2 && any(na.omit(impute.inds))){
    .data.frame <- na.omit(data.frame(y=Nk.estimates, 
                                      x=n.k.counts, 
                                      conv.ind=conv.inds,
                                      impute.ind=impute.inds))
    lm.2 <- rlm(y ~ x - 1, data=.data.frame, subset=conv.inds , maxit=100L)
    Nk.estimates[which(impute.inds)]  <- predict(lm.2, newdata = subset(.data.frame, impute.inds))  
  }
  return(Nk.estimates)                                                
}








## Estimate theta assuming beta_k=beta * k^theta:
getTheta <- function(rds.object, bin=1, robust=TRUE){
  nk.estimates <- rds.object$estimates
  
  log.bks <- nk.estimates$log.bk.estimates
  log.ks <- log(seq_along(log.bks)*bin)
  ok.inds <- is.finite(log.bks)& !is.na(log.bks)
  
  if(robust) {
    lm.1 <- rlm(log.bks[ok.inds]~log.ks[ok.inds])
  } else {
    lm.1 <- lm(log.bks~log.ks)
  }
  
  coefs <- as.list(coef(lm.1) )
  
  return(list(
    log.beta_0 = coefs$`(Intercept)`,
    theta = coefs$log.ks,
    model=lm.1))
}


## Recovering Nk with smoothed Nk:
thetaSmoothingNks <- function(rds.object,...){
  theta <- getTheta(rds.object, ...)
  smooth.bks <- exp(predict(theta$model))
  A.ks <- na.omit(rds.object$estimates$A.ks)
  B.ks <- na.omit(rds.object$estimates$B.ks)
  n.k.counts  <- na.omit(rds.object$estimates$n.k.counts)
  smooth.Nks <- (n.k.counts/smooth.bks + B.ks)/A.ks
  
  ## FIXME: 
  # if log.b.ks is infinite, smooth.Nks cannot be computed because
  # A.ks and B.ks will exist but smooth.bks will not. 
  
  return(smooth.Nks)
}
## Testing:
# thetaSmoothingNks(rds.object)






likelihood <- function(
  log.bk, Nk.estimates, I.t, 
  n.k.counts, degree.in, degree.out, 
  arrival.intervals, arrival.degree){
  ### Verification:
  
  ### Initialization:
  uniques <- which(!is.na(n.k.counts))
  n.k.t <- makeNKT(uniques, degree.in, degree.out)
  betas <- exp(log.bk)
  
  
  ## Computation
  result <- 0.
  for(i in seq_along(arrival.degree)){
    if(i==1) next()
    for(j in seq_along(uniques)){ 
      #       i <- 5; j <- 5
      k <- uniques[[j]]
      lambda <-  betas[k] * (Nk.estimates[k] - n.k.t[j,i-1]) * I.t[i-1]
      lambda <- max(lambda, .Machine$double.eps) 
      
      A <- ifelse(arrival.degree[i]==k, log(lambda), 0) 
      B <- lambda * arrival.intervals[i-1] 
      result <- result + A - B 
    }
  }
  result <- -result
  return(result)
}
## Testing:
# source('temp/Uganda_example.R')
# chords:::likelihood(rds.object$estimates$log.bk.estimates, 
#                     rds.object$estimates$Nk.estimates, 
#                     rds.object$I.t, 
#                     rds.object$estimates$n.k.counts, 
#                     rds.object$degree.in, 
#                     rds.object$degree.out, 
#                     rds.object$estimates$arrival.intervals, 
#                     rds.object$estimates$arrival.degree,
#                     const = 100)

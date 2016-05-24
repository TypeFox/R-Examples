.beta.i <-
function(z, params, sd.t, sd.tb, lower.limit, upper.limit){
  return(params[1]*exp(-(z/sd.tb)^2)*(pnorm((upper.limit*sd.tb - params[1]^2*z/sd.tb)/(params[1]*sd.t)) - pnorm((lower.limit*sd.tb - params[1]^2*z/sd.tb)/(params[1]*sd.t)))/sd.tb  )
}

.cfun1 <-
function(theta,ratings.1,ratings.2,i.limits){  
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  f.table <- mat.or.vec(n_classes, n_classes)
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
                
      intg <- .f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
      
      if(intg < 0){
        stop('f.table < 0')}else{
                f.table[lb1, lb2] <- - log(intg)
        }
      f.table[lb2, lb1] <- f.table[lb1, lb2]
    }
  } 
    
  nonzero <- table(ratings.1,ratings.2)>0
  f = f + sum(table(ratings.1,ratings.2)[nonzero]*f.table[nonzero]);   
  return(f);
}

.cfun2 <-
function(theta,ratings.1,ratings.2,i.limits){  
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  f.table <- mat.or.vec(n_classes, n_classes)
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
                
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }else{
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }
      
      if(intg < 0){
        stop('f.table < 0')}else{
                f.table[lb1, lb2] <- - log(intg)
        }
      f.table[lb2, lb1] <- f.table[lb1, lb2]
    }
  } 
    
  nonzero <- table(ratings.1,ratings.2)>0
  f = f + sum(table(ratings.1,ratings.2)[nonzero]*f.table[nonzero]);   
  return(f);
}

.cfun3 <-
function(theta,ratings.1,ratings.2){  
    
  n = nrow(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  for (i.n in 1:n){
    
    i1.lb <- ratings.1[i.n,1]
    i1.ub <- ratings.1[i.n,2]
    i2.lb <- ratings.2[i.n,1]
    i2.ub <- ratings.2[i.n,2]     
                
    intg <- .f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
    
    if(intg < 0){
      stop('f.table < 0')}else{
        f <- f - log(intg)
      }
  }    
  return(f);
}

.cfun4 <-
function(theta,ratings.1,ratings.2){  
    
  n = nrow(ratings.1);
  f = n*log(theta[1]);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  
  for (i.n in 1:n){
    
    i1.lb <- ratings.1[i.n,1]
    i1.ub <- ratings.1[i.n,2]
    i2.lb <- ratings.2[i.n,1]
    i2.ub <- ratings.2[i.n,2]     
                
    if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
      intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
      }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
        intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
        }else{
          intg <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
    
    if(intg < 0){
      stop('f.table < 0')}else{
        f <- f - log(intg)
      }
  }    
  return(f);
}

.chi2.int.test <-
function(d.lower,d.upper,bins = 10,mu,stdev){
  
  # test if interval-censored data follows normal distribution with mean mu and standard deviation stdev
  # using chi-squared test
  #
  # INPUT:
  # d.lower - lower limits of data censoring intervals
  # d.upper - upper limits of data censoring intervals
  # bins - number of categories in chi-squared test
  # mu - mean of baseline distribution
  # stdev - standard deviation of baseline distribution
  
  breaks <- qnorm(seq(0,1,1/bins),mean=mu,sd=stdev)
  nresp <- length(d.lower) 
  lims <- sort(union(d.lower,d.upper))
  dis_lim <- sort(union(lims, breaks))
  
  n_counts <- length(dis_lim)-1
  dis_prob <- rep(0, n_counts)
  mids <- rep(0, n_counts)
  for(i in 1:n_counts){
    dis_prob[i] <- pnorm(dis_lim[i+1], mean=mu, sd=stdev) - pnorm(dis_lim[i], mean=mu, sd=stdev)  
    mids[i] <- 0.5*(dis_lim[i] + dis_lim[i+1])
  }

  counts <- rep(0,n_counts)
  vals <- rep(0,nresp)

  for(i in 1:nresp){
    ind1 <- which(dis_lim==d.lower[i])
    ind2 <- which(dis_lim==d.upper[i])
    normf <- sum(dis_prob[ind1:(ind2-1)])
    for(j in ind1:(ind2-1)){
      counts[j] <- counts[j] + dis_prob[j]/normf
    }  
  }
  
  n_counts.eq <- length(breaks)-1
  counts.eq <- rep(0, n_counts.eq)
  dis_prob.eq <- rep(0, n_counts.eq)
  mids.eq <- rep(0, n_counts.eq)
  
  for(i in 1:n_counts.eq){
    ind1 <- which(dis_lim == breaks[i])
    ind2 <- which(dis_lim == breaks[i+1])
    mids.eq[i] <- 0.5*(breaks[i] + breaks[i+1])
    counts.eq[i] <- sum(counts[ind1:(ind2-1)])
    dis_prob.eq[i] <- sum(dis_prob[ind1:(ind2-1)])
  }  
  
  t <- chisq.test(c(counts.eq),p = c(dis_prob.eq))  
  return(t);
}

.f.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){  
  I <- function(u){
    return(exp(-(u/params[1])^2)*pnorm((sgn*u + z - params[3])/sd.t));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}

.fmod.int <-
function(z1, sgn1, z2, sgn2, params, sd.t, lower.limit, upper.limit){   
  I <- function(u){
    return(exp(-(u/params[1])^2)*(pnorm((sgn1*u + z1 - params[3])/sd.t) - pnorm((sgn2*u + z2 - params[3])/sd.t)));
  } 
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}

.gradcfun1 <-
function(theta,ratings.1,ratings.2,i.limits){
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
 
  Jf1.table <- mat.or.vec(n_classes, n_classes)
  Jf2.table <- mat.or.vec(n_classes, n_classes)
  Jf3.table <- mat.or.vec(n_classes, n_classes)
  
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
      
      corr.n <- (.f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)));

      Jf1.table[lb1, lb2] <- - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf1.table[lb2, lb1] <- Jf1.table[lb1, lb2]
      Jf2.table[lb1, lb2] <- - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf2.table[lb2, lb1] <- Jf2.table[lb1, lb2]
      Jf3.table[lb1, lb2] <- (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf3.table[lb2, lb1] <- Jf3.table[lb1, lb2]
    }
  }
  nonzero <- table(ratings.1,ratings.2)>0
  Jf[1] = Jf[1] + sum(table(ratings.1,ratings.2)[nonzero]*Jf1.table[nonzero]);
  Jf[2] = Jf[2] + sum(table(ratings.1,ratings.2)[nonzero]*Jf2.table[nonzero]);
  Jf[3] = Jf[3] + sum(table(ratings.1,ratings.2)[nonzero]*Jf3.table[nonzero]);
  
  return(Jf);
}

.gradcfun2 <-
function(theta,ratings.1,ratings.2,i.limits){
  
  n_classes = nrow(i.limits);
  n = length(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
 
  Jf1.table <- mat.or.vec(n_classes, n_classes)
  Jf2.table <- mat.or.vec(n_classes, n_classes)
  Jf3.table <- mat.or.vec(n_classes, n_classes)
  
  for (lb1 in 1:n_classes){
   for (lb2 in lb1:n_classes){
      i1.lb <- i.limits[lb1,1]
      i1.ub <- i.limits[lb1,2]
      i2.lb <- i.limits[lb2,1]
      i2.ub <- i.limits[lb2,2]      
      
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
        corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }else{
            corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
          

      Jf1.table[lb1, lb2] <- - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf1.table[lb2, lb1] <- Jf1.table[lb1, lb2]
      Jf2.table[lb1, lb2] <- - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf2.table[lb2, lb1] <- Jf2.table[lb1, lb2]
      Jf3.table[lb1, lb2] <- (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf3.table[lb2, lb1] <- Jf3.table[lb1, lb2]
    }
  }
  nonzero <- table(ratings.1,ratings.2)>0
  Jf[1] = Jf[1] + sum(table(ratings.1,ratings.2)[nonzero]*Jf1.table[nonzero]);
  Jf[2] = Jf[2] + sum(table(ratings.1,ratings.2)[nonzero]*Jf2.table[nonzero]);
  Jf[3] = Jf[3] + sum(table(ratings.1,ratings.2)[nonzero]*Jf3.table[nonzero]);  
  
  return(Jf);
}

.gradcfun3 <-
function(theta,ratings.1,ratings.2){
  
  n = nrow(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
    
  for (i.n in 1:n){
      i1.lb <- ratings.1[i.n,1]
      i1.ub <- ratings.1[i.n,2]
      i2.lb <- ratings.2[i.n,1]
      i2.ub <- ratings.2[i.n,2]
      
      corr.n <- (.f.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .f.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .f.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .f.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)));

      Jf[1] <- Jf[1] - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[2] <- Jf[2] - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[3] <- Jf[3] + (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
       
  }  
  return(Jf);
}

.gradcfun4 <-
function(theta,ratings.1,ratings.2){
  
  n = nrow(ratings.1);
  Jf = c(n/theta[1], 0, 0);  
  sdev.th = sqrt(0.5*theta[1]^2 + theta[2]^2);  
  sdev.thb = sqrt(2*(theta[1]^2 + theta[2]^2));  
    
  for (i.n in 1:n){
      i1.lb <- ratings.1[i.n,1]
      i1.ub <- ratings.1[i.n,2]
      i2.lb <- ratings.2[i.n,1]
      i2.ub <- ratings.2[i.n,2]
      
      if((i1.ub - i2.ub)<(i1.lb - i2.lb)){
        corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i1.lb,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb))
        }else if((i1.ub - i2.ub)>(i1.lb - i2.lb)){
          corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) + .fmod.int(i2.ub,1,i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }else{
            corr.n <- .fmod.int(i2.ub,1,i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .fmod.int(i1.ub,-1,i2.lb,1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb))
          }
          
      Jf[1] <- Jf[1] - (.sigma.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigma.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigma.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigma.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[2] <- Jf[2] - (.sigmab.int(i2.ub,1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .sigmab.int(i1.ub,-1,theta,sdev.th,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .sigmab.int(i1.lb,-1,theta,sdev.th,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .sigmab.int(i2.lb,1,theta,sdev.th,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;
      Jf[3] <- Jf[3] + (.beta.i(theta[3]-i2.ub,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.ub-i2.ub)) + .beta.i(i1.ub-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.ub-i2.ub),0.5*(i1.ub-i2.lb)) - .beta.i(i1.lb-theta[3],theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.ub),0.5*(i1.lb-i2.lb)) - .beta.i(theta[3]-i2.lb,theta,sdev.th,sdev.thb,0.5*(i1.lb-i2.lb),0.5*(i1.ub-i2.lb)))/corr.n;       
  }  
  return(Jf);
}

.intervalICC.est1 <-
function(ratings, classes, c.limits, theta0){

  costFunc <- function(theta){
    return(.cfun1(theta,ratings$t1,ratings$t2,c.limits))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun1(theta,ratings$t1,ratings$t2,c.limits))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=2"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}

.intervalICC.est2 <-
function(ratings, classes, c.limits, theta0){

  costFunc <- function(theta){
    return(.cfun2(theta,ratings$t1,ratings$t2,c.limits))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun2(theta,ratings$t1,ratings$t2,c.limits))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}

.intervalICC.est3 <-
function(ratings1, ratings2, theta0){

  costFunc <- function(theta){
    return(.cfun3(theta,ratings1,ratings2))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun3(theta,ratings1,ratings2))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings1)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}

.intervalICC.est4 <-
function(ratings1, ratings2, theta0){

  costFunc <- function(theta){
    return(.cfun4(theta,ratings1,ratings2))
  }
  grad.costFunc <- function(theta){
    return(.gradcfun4(theta,ratings1,ratings2))
  }
  
  th_optim <- tryCatch({constrOptim(theta=theta0, f=costFunc, grad=grad.costFunc, ui=rbind(c(1,0,0),c(0,1,0)), ci=c(.Machine$double.eps,.Machine$double.eps), method = "BFGS")}, error = function(e) {print("Numerical optimization failed; try optim.method=1"); return(0)})
    
  sigma2.b <- th_optim$par[2]^2
  sigma2.w <- th_optim$par[1]^2
  mu <- th_optim$par[3]
  icc <- sigma2.b/(sigma2.b + sigma2.w)
  loglikelihood <- -th_optim$value - nrow(ratings1)*log(pi)/2
  
  list(icc = icc,
       sigma2.b = sigma2.b,
       sigma2.w = sigma2.w,
       mu = mu,
       loglikelihood = loglikelihood)
}

.sigma.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){  
  I <- function(u){
    arg.in = sgn*u + z - params[3];    
    return(exp(-(u/params[1])^2)*(2*(u^2/params[1]^3)*pnorm(arg.in/sd.t) - params[1]*arg.in/(2*sd.t^3)*dnorm(arg.in/sd.t)));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}

.sigmab.int <-
function(z, sgn, params, sd.t, lower.limit, upper.limit){ 
  I <- function(u){
    arg.in = sgn*u + z - params[3];    
    return(exp(-(u/params[1])^2)*(-params[2]*arg.in*dnorm(arg.in/sd.t)/sd.t^3));
  }
  return(integrate(I,lower.limit,upper.limit,stop.on.error=1)$value);
}


stheoremme <- function(distribution0, distribution1, criterionly=FALSE){
  
  # Date: 13-02-2015
  
  # Create and check some items
  f0 <- distribution0/sum(distribution0)
  f1 <- distribution1/sum(distribution1)
  if (abs(length(f1)-length(f0))>0) {
    stop('!Distributions have not equal lenghts!')
  }
  if (length(f0[which(f0<0)])>0) {
    stop('!Distribution0 is not valid: probably contains negatives')
  }
  if (length(f1[which(f1<0)])>0) {
    stop('!Distribution1 is not valid: probably contains negatives')
  }
  #if (alpha>=1 | alpha<=0) {stop('Alpha value must be 0<|alpha|<1!')}
  
  #new distribution functions with zeros repaired
  f0[which(f0<1e-15)]=1e-15 
  f0 = f0/sum(f0)
  f1[which(f1<1e-15)]=1e-15
  f1 = f1/sum(f1)
  
  ### func entropy calulation
  dS <- function(f0,f1){
    
    #new distribution functions with zeros repaired
    f0[which(f0<1e-15)]=1e-15 
    f0 = f0/sum(f0)
    f1[which(f1<1e-15)]=1e-15
    f1 = f1/sum(f1)
    
    # zero step of analysis
    f0s <- f0^(1/1); f0s = f0s/sum(f0s)
    dH0to1 <- sum(-f0s*log(f0))-sum(-f1*log(f0))
    f1s <- f1^(1/1); f1s = f1s/sum(f1s)
    dH1to0 <- sum(-f1s*log(f1))-sum(-f0*log(f1))
    greenFlag <- TRUE
    D <- 1001 #D is temperature
    
    #check if s-criterion works at all
    if (dH0to1*dH1to0 > 0) { 
      greenFlag = FALSE
      dS <- 1001 #absurd high if criterion doesn't work
    } 
    if (dH0to1<0 & greenFlag) {
      Tm <- 1+c(0:200)/10 #{if Ham.mean diff is negative => D up(?)}
      dH <- NaN
      for (i in 1:200) {
        D <- Tm[i]
        f0s <- f0^(1/D); f0s = f0s/sum(f0s)
        dH[i] = abs(sum(-f0s*log(f0))-sum(-f1*log(f0)))
      }
      D=Tm[which(dH==min(dH))]
      f0s <- f0^(1/D); f0s = f0s/sum(f0s)
      dS <- sum(-f1*log(f1))-sum(-f0s*log(f0s))
    }
    if (dH0to1>0 & greenFlag) {
      Tm <- 1+c(0:200)/10
      dH <- NaN
      for (i in 1:200) {
        D <- Tm[i]
        f1s <- f1^(1/D); f1s = f1s/sum(f1s)
        dH[i] = abs(sum(-f1s*log(f1))-sum(-f0*log(f1)))
      }
      D=Tm[which(dH==min(dH))]
      f1s <- f1^(1/D); f1s = f1s/sum(f1s)
      dS <- sum(-f1s*log(f1s))-sum(-f0*log(f0))
    }
    if (abs(dH0to1)==0 | abs(dH1to0)==0) {
      dS <- 0
      D=1
    } 
    
    return(list(greenFlag=greenFlag,val=dS,D=D))
  }
  ###
  
  # P_value when if direct evolution is possible:
  p_val <- 1
  
  #search for the medium state if direct evolution is forbidden
  ds <- dS(f0,f1) 
  ref0_index <- 50 #initial params
  ref1_index <- 50
  ref0_ds <- ds$val
  ref1_ds <- 0
  
  if (ds$greenFlag == FALSE) {
    i <- c(1)
    while (i<100) { #generate new functions in between 
      fbtw1 <- ((i)*f0+(100-i)*f1)/100; fbtw1=fbtw1/sum(fbtw1) 
      ds_forward <- dS(f0,fbtw1) 
      if (ds_forward$greenFlag == FALSE) {
        i=i+1
      } else{     
        ref0_index <- i
        ref0_ds <- ds_forward$val
        i=101
      }  
    }
    i <- c(1)
    while (i<100) {#generate new functions in between 
      fbtw1 <- ((i)*f1+(100-i)*f0)/100; fbtw1=fbtw1/sum(fbtw1) 
      ds_back <- dS(fbtw1,f1) 
      if (ds_back$greenFlag  == FALSE) {
        i=i+1
      } else {
        ref1_index <- 100-i
        ref1_ds <- ds_back$val
        i=101
      }
    }
    
    if ((ref0_index-ref1_index) == 1) {
      p_val = (max(abs(50-ref0_index),abs(50-ref1_index))/50)^2
      #a bit speculative residual function    
    }
    
    if ((ref0_index-ref1_index) < 1) {
      #need to calculate fbtw1 again to match the indices to be at 
      #the same medium point
      fbtw1 <- ((ref0_index)*f0+(100-ref0_index)*f1)/100; fbtw1=fbtw1/sum(fbtw1) 
      ds_back <- dS(fbtw1,f1)
      ref1_ds <- ds_back$val
      p_val = (max(abs(50-ref0_index),abs(50-ref1_index))/50)^2
      #a bit speculative residual function  
    }
    
    if ((ref0_index-ref1_index) > 1) {
      ref0_ds = 1001
      p_val = 0
    }
  } 
  
  ###
  
  #outputs:
  dS_02 <- ref0_ds;  if (dS_02==1001) {dS_02=1e-15} #protection patch
  dS_21 <- ref1_ds;  if (dS_21==1001) {dS_21=1e-15} #protection patch
  dH_val <- -sum(f1*log(f1))+sum(f0*log(f0))
  dS_val <- dS_02+dS_21;  if (dS_val==1001) {dS_val=1e-15} #protection patch
  dH_ext <- dH_val-dS_val
  
  is_direct <- ds$greenFlag
  
  # Protection from abnormals when testing
  if (p_val<0) {p_val=1e-15}
  p_val = p_val
  ###
  
  ref_indices <- c(ref0_index,ref1_index)
  
  if (criterionly==FALSE) {
    formula <- paste(as.character(dH_val),c('{0to1} ='),as.character(dS_val),
                     c('{0to1} +'),as.character(dH_ext),c('{0to1}'))
  } else formula <- {'criterionly'}
  #identity <- 1
  #if (p_val < alpha) {identity = 0}
  
  
  return(list(is_direct=is_direct, r2_val = p_val, dH_val=dH_val, dS_val=dS_val,  
              dH_ext=dH_ext, dS_02=dS_02, dS_21=dS_21, formula=formula,
              ref_indices=ref_indices))
  
}


stheorem<- function(distribution0, distribution1, criterionly) UseMethod("stheorem")
stheorem.default <- function(distribution0, distribution1, criterionly=FALSE) {
  
  distribution0 <- as.vector(distribution0)
  distribution1 <- as.vector(distribution1)
  est <- stheoremme(distribution0, distribution1, criterionly)
  est$is_direct <- as.logical(est$is_direct)
  est$r2_val <- as.numeric(est$r2_val)
  est$dH_val <- as.numeric(est$dH_val)  
  est$dS_val <- as.numeric(est$dS_val)
  est$dH_ext <- as.numeric(est$dH_ext) 
  est$dS_02 <- as.numeric(est$dS_02)
  est$dS_21 <- as.numeric(est$dS_21)
  est$formula <- as.character(est$formula)
  est$ref_indices <- as.vector(est$ref_indices)
  
  est$call <- match.call()
  class(est) <- "stheorem"
  est
}

print.stheorem <- function(x,...){
  
  if (x$formula=='criterionly') {
    cat("\n")  
    cat(" S-theorem convergence criterion\n")
    cat("\nSystem evolution from state0 to state1 is thermodynamically")
    if (x$r2_val==1) {
      cat(" allowed ")
      cat("(R^2 = ")
      cat(as.character(x$r2_val)); cat(").\n")
    }
    if (x$r2_val<0.00043) {
      cat(" forbidden \n")
      cat("(R^2 = 0).\n")
    }
    if (x$r2_val<1 & x$r2_val>=0.00043) {
      cat(" possible through an indirect\n medium state2 ") 
      cat("(R^2 = ")
      cat(as.character(x$r2_val)); cat(").\n")
    }
  }
  
  ###
  if (x$formula!='criterionly') {
    if (x$r2_val<=1 & x$r2_val>=0.00043) {
      cat("\n S-theorem open system evolution model\n")
      cat("\n") 
      cat("Overall entropy shift ")
      cat("{H1-H0 = dS + dI}:\n") 
      cat(as.character(x$formula))
    }
    if (x$r2_val<1 & x$r2_val>=0.00043) {
      cat("\nwhere dS consists of two:\n")
      cat(as.character(x$dS_val),c('{0to1} ='),as.character(x$dS_02),
          c('{0to2} +'),as.character(x$dS_21),c('{2to1}\n'))
    }
    if (x$r2_val<0.00043) {
      cat("\nSystem evolution is thermodynamically forbidden. No valid model\n")
    }
  }
  cat("\n")   
}

crit.stheorem <- function(distribution0, distribution1) {
  return(stheorem(distribution0, distribution1, criterionly=TRUE))
}

cxds.stheorem <- function(distribution0, distribution1) {
  return(stheorem(distribution0, distribution1, criterionly=FALSE))
}

d1char.d1nat <- function(farr0, farr1, reject=c("")) {
  
  # Date: 23-02-2015
  
  for (i in 1:length(reject)) {
    farr0=farr0[which(farr0!=reject[i])]
    farr1=farr1[which(farr1!=reject[i])]    
  }
  
  cat("Separate objects with corresponding object numbers:\n")
  print(t(levels(as.factor(c(farr0,farr1)))))
  
  out <- d1nat(as.numeric(as.factor(c(farr0,farr1)))[1:length(farr0)], 
               as.numeric(as.factor(c(farr0,farr1)))[(length(farr0)+1):(length(farr0)+length(farr1))], 
               brks=max(as.numeric(as.factor(c(farr0,farr1)))))
  return(out)
}

d1native <- function(sample0, sample1, band=c(0,0), brks=0){
  
  # Date: 12-02-2015
  
  # set the window & adjust data to it:
  if (band[2]==band[1] & band[2]==0) {
    band = c(min(sample1, sample0), max(sample1, sample0))
  }
  data_0 = sample0[which(sample0>=band[1] & sample0<=band[2])]
  data_1  = sample1[which(sample1>=band[1] & sample1<=band[2])]
  if (length(data_1)>1 & length(data_0)>1) {
    full_range = (max(data_1, data_0)-min(data_1, data_0))
  } else stop("Band size doesn't match!")
  
  # set the breaks
  if (brks == 0) {
    num_breaks = floor(sqrt(length(data_0)))  #dance from zero!
    break_size = (max(data_0)-min(data_0))/num_breaks
    num_breaks = floor(full_range/break_size)
  } else {num_breaks = brks}
  break_size = full_range/num_breaks
  the_breaks = min(data_1, data_0)+c(0:num_breaks)*break_size
  
  #mid points
  the_mids <- the_breaks[1]-break_size/2+c(1:num_breaks)*break_size
  
  #distributions final - as quazi-histogram grouping
  f0 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f0[i] = sum(which(data_0>the_breaks[i] & data_0<=the_breaks[i+1])/
                  which(data_0>the_breaks[i] & data_0<=the_breaks[i+1]))
  }
  f0[1]=f0[1]+sum(which(data_0==the_breaks[1])/
                    which(data_0==the_breaks[1])) #take also zeros missed previously
  #
  f1 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f1[i] = sum(which(data_1>the_breaks[i] & data_1<=the_breaks[i+1])/
                  which(data_1>the_breaks[i] & data_1<=the_breaks[i+1]))
  }
  f1[1]=f1[1]+sum(which(data_1==the_breaks[1])/
                    which(data_1==the_breaks[1])) #take also zeros missed previously
  #
  f0=f0/sum(f0)
  f1=f1/sum(f1)
  
  # absolute chaotic entropy calc
  fneutral <- rep(1,length(f0))/length(f0)
  entropy_max <- -sum(log(fneutral)*fneutral)
  
  #statistic summary for new distributions
  f0_n_expected = sum(f0*the_mids)
  f0_n_stdev = sqrt(sum(f0*the_mids^2)-f0_n_expected^2)
  f0_n_entropy = -sum(f0*log(f0+1e-15))
  if (f0_n_entropy<0) {f0_n_entropy=0} #protection patch
  f0_n_1stmod = min(the_mids[which(f0==max(f0))]) #use 1st if where is > 1 modes
  ftemp<-f0; ftemp[which(the_mids==f0_n_1stmod)]=0 #temp without 1st mode 
  f0_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f0_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f0_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  #
  f1_n_expected = sum(f1*the_mids)
  f1_n_stdev = sqrt(sum(f1*the_mids^2)-f1_n_expected^2)
  f1_n_entropy = -sum(f1*log(f1+1e-15))
  if (f1_n_entropy<0) {f1_n_entropy=0} #protection patch
  f1_n_1stmod = min(the_mids[which(f1==max(f1))]) #use 1st if where is > 1 modes
  ftemp<-f1; ftemp[which(the_mids==f1_n_1stmod)]=0 #temp without 1st mode 
  f1_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f1_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f1_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  
  stat_summary <- cbind(expctd=c(f0=f0_n_expected,f1=f1_n_expected),
                        var=c(f0=f0_n_stdev^2,f1=f1_n_stdev^2), 
                        fsum=c(f0=sum(f0),f1=sum(f1)), 
                        xmin=c(f0=min(the_breaks),f1=min(the_breaks)), 
                        xmax=c(f0=max(the_breaks),f1=max(the_breaks)), 
                        n=c(f0=num_breaks,f1=num_breaks),
                        mod1=c(f0=f0_n_1stmod,f1=f1_n_1stmod),
                        mod2=c(f0=f0_n_2ndmod,f1=f1_n_2ndmod),
                        mod3=c(f0=f0_n_3rdmod,f1=f1_n_3rdmod),
                        H_val=c(f0=f0_n_entropy,f1=f1_n_entropy),
                        H_max=c(f0=entropy_max,f1=entropy_max))
  
  essence='Native distributions'
  
  #many happy returns!
  return(list(f0=f0, f1=f1, midpoints=the_mids, 
              stat_summary=stat_summary, essence=essence))
}


d1nat<- function(sample0, sample1, band=c(0,0), brks=0) UseMethod("d1nat")
d1nat.default <- function(sample0, sample1, band=c(0,0), brks=0) {
  
  sample0 <- as.vector(sample0)
  sample1 <- as.vector(sample1)
  band <- as.vector(band)
  brks <- as.numeric(brks)
  est <- d1native(sample0, sample1, band, brks)
  est$f0 <- as.vector(est$f0)
  est$f1 <- as.vector(est$f1)
  est$midpoints <- as.vector(est$midpoints)
  est$stat_summary <- as.matrix(est$stat_summary)
  est$essence <- as.character(est$essence)
  est$call <- match.call()
  class(est) <- "d1nat"
  est
}

print.d1nat <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Two probability mass functions ($f0,$f1) have ")
  cat("been generated at the common scale of values ($midpoints)\n")
  cat("Statistics summary:\n")
  print(x$stat_summary)
  plot(x)
}

plot.d1nat <- function(x,...){
  
  min_val <- x$stat_summary[1,4]
  max_val <- x$stat_summary[1,5]
  the_mids <- x$midpoints
  the_mids_1 <- the_mids+(max_val-min_val)/(x$stat_summary[1,6])/7
  plot(-1,-1, xlim = c(min_val,max_val), 
       ylim = c(0,(1.4*max(x$f0,x$f1))), ylab="Probability Mass", 
       xlab = "midpoints")
  lines(the_mids, x$f0, type='h', lwd = 3, col=1)
  lines(the_mids_1, x$f1, type='h', lwd = 3, col=8)
  legend(mean(the_mids),(1.4*max(x$f0,x$f1)), c("$f0 to sample0","$f1 to sample1"),
         lty=c(1,1), col=c(1,8))
  title(x$essence)
}

d1spectrum <- function(sample0, sample1, band=c(0,0), brks=0, meansub=TRUE) {
  
  # Date: 08-02-2015
  
  #mean subtraction
  if (meansub==TRUE){
    sample0 = sample0-mean(sample0)
    sample1 = sample1-mean(sample1)
  }
  
  # power spectrum and freq spaces calc
  pwR0 <- abs(fft(sample0))[1:(floor(length(sample0)/2))]
  freq0 <- 0.5*c(1:length(pwR0))/length(pwR0)
  pwR1 <- abs(fft(sample1))[1:(floor(length(sample1)/2))]
  freq1 <- 0.5*c(1:length(pwR1))/length(pwR1)
  
  # set the window & adjust data to it:
  if (band[2]==band[1] & band[2]==0) {
    band = c(0, 0.5) # freq band by default: full
  }
  
  #check the band
  if (min(band[1],band[2])<0) {stop("Wrong band!")}  
  if (max(band[1],band[2])>0.5) {stop("Wrong band!")}  
  
  #choose powers and freqs from within band only
  pwR0 = pwR0[which(freq0>=band[1] & freq0<=band[2])]
  freq0 = freq0[which(freq0>=band[1] & freq0<=band[2])]
  pwR1  = pwR1[which(freq1>=band[1] & freq1<=band[2])]
  freq1  = freq1[which(freq1>=band[1] & freq1<=band[2])]
  if (length(pwR0)>1 & length(pwR1)>1) {
    full_range = band[2]-band[1]
  } else stop("Band size doesn't match!")
  
  # set the breaks
  if (brks == 0) {
    num_breaks = min(length(freq0),length(freq1)) # choose as min of 2 posibilities
    break_size = full_range/num_breaks
    num_breaks = floor(full_range/break_size)
  } else {num_breaks = brks}
  break_size = full_range/num_breaks
  the_breaks = band[1]+c(0:num_breaks)*break_size
  
  #midpoints
  the_mids <- band[1]-break_size/2+c(1:num_breaks)*break_size
  
  #final quazi-histogramming distributions from within the breaks
  f0 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f0[i] = sum(pwR0[which(freq0>the_breaks[i] & freq0<=the_breaks[i+1])])
  }
  f0[1]=f0[1]+sum(pwR0[which(
    freq0==the_breaks[1])]) #take also zeros missed previously
  #
  f1 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f1[i] = sum(pwR1[which(freq1>the_breaks[i] & freq1<=the_breaks[i+1])])
  }
  f1[1]=f1[1]+sum(pwR1[which(
    freq1==the_breaks[1])]) #take also zeros missed previously
  #
  
  if (sum(f0)==0 | sum(f1)==0) {
    stop('One of the samples is a constant. Analysis is not relevant!')
  }
  
  f0=f0/sum(f0)
  f1=f1/sum(f1)
  
  # absolute chaotic entropy calc
  fneutral <- rep(1,length(f0))/length(f0)
  entropy_max <- -sum(log(fneutral)*fneutral)
  
  #statistic summary for new distributions
  f0_n_expected = sum(f0*the_mids)
  f0_n_stdev = sqrt(sum(f0*the_mids^2)-f0_n_expected^2)
  f0_n_entropy = -sum(f0*log(f0+1e-15))
  if (f0_n_entropy<0) {f0_n_entropy=0} #protection patch
  f0_n_1stmod = min(the_mids[which(f0==max(f0))]) #use 1st if where is > 1 modes
  ftemp<-f0; ftemp[which(the_mids==f0_n_1stmod)]=0 #temp without 1st mode 
  f0_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f0_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f0_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  #
  f1_n_expected = sum(f1*the_mids)
  f1_n_stdev = sqrt(sum(f1*the_mids^2)-f1_n_expected^2)
  f1_n_entropy = -sum(f1*log(f1+1e-15))
  if (f1_n_entropy<0) {f1_n_entropy=0} #protection patch
  f1_n_1stmod = min(the_mids[which(f1==max(f1))]) #use 1st if where is > 1 modes
  ftemp<-f1; ftemp[which(the_mids==f1_n_1stmod)]=0 #temp without 1st mode 
  f1_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f1_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f1_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  
  stat_summary <- cbind(expctd=c(f0=f0_n_expected,f1=f1_n_expected),
                        var=c(f0=f0_n_stdev^2,f1=f1_n_stdev^2), 
                        fsum=c(f0=sum(f0),f1=sum(f1)), 
                        xmin=c(f0=min(the_breaks),f1=min(the_breaks)), 
                        xmax=c(f0=max(the_breaks),f1=max(the_breaks)), 
                        n=c(f0=num_breaks,f1=num_breaks),
                        mod1=c(f0=f0_n_1stmod,f1=f1_n_1stmod),
                        mod2=c(f0=f0_n_2ndmod,f1=f1_n_2ndmod),
                        mod3=c(f0=f0_n_3rdmod,f1=f1_n_3rdmod),
                        H_val=c(f0=f0_n_entropy,f1=f1_n_entropy),
                        H_max=c(f0=entropy_max,f1=entropy_max))
  
  #many happy returns!
  return(list(f0=f0, f1=f1, freqs=the_mids, 
              stat_summary=stat_summary, essence='Spectral distributions'))
}

d1spec<- function(sample0, sample1, band=c(0,0), brks=0, meansub=TRUE) UseMethod("d1spec")
d1spec.default <- function(sample0, sample1, band=c(0,0), brks=0, meansub=TRUE) {
  
  sample0 <- as.vector(sample0)
  sample1 <- as.vector(sample1)
  band <- as.vector(band)
  brks <- as.numeric(brks)
  meansub <- as.logical(meansub)
  est <- d1spectrum(sample0, sample1, band, brks, meansub)
  est$f0 <- as.vector(est$f0)
  est$f1 <- as.vector(est$f1)
  est$freqs <- as.vector(est$freqs)
  est$stat_summary <- as.matrix(est$stat_summary)
  est$essence <- as.character(est$essence)
  est$call <- match.call()
  class(est) <- "d1spec"
  est
}

print.d1spec <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Two spectral probability mass functions ($f0,$f1) have ")
  cat("been generated at the common scale of values ($freqs)\n")
  cat("Statistics summary:\n")
  print(x$stat_summary)
  plot(x)
}

plot.d1spec <- function(x,...){
  min_val <- x$stat_summary[1,4]
  max_val <- x$stat_summary[1,5]
  the_mids <- x$freqs
  the_mids_1 <- the_mids+(max_val-min_val)/(x$stat_summary[1,6])/7
  plot(-1,-1, xlim = c(min_val,max_val), 
       ylim = c(0,(1.4*max(x$f0,x$f1))), ylab="Probability Mass", 
       xlab = "Frequency, 1/(sampling_interval)")
  lines(the_mids, x$f0, type='h', lwd = 3, col=1)
  lines(the_mids_1, x$f1, type='h', lwd = 3, col=8)
  legend(mean(the_mids),(1.4*max(x$f0,x$f1)), c("$f0 to sample0","$f1 t0 sample1"),
         lty=c(1,1), col=c(1,8))
  title(x$essence)
}

d2nat.d1nat <- function(d2arr0, d2arr1, band=c(0,0), brks=64, method='default'){
  
  # Date: 08-02-2015
  
  #3 methods: default image hist, cols and rows
  d2arr0.as.d1 <- c(d2arr0) #default
  d2arr1.as.d1 <- c(d2arr1) #default
  if (method=='cols'){
    d2arr0.as.d1 <- NaN
    for (i in 1:length(d2arr0[1,])){
      d2arr0.as.d1[i]=mean(d2arr0[,i]) 
    }
    d2arr1.as.d1 <- NaN
    for (i in 1:length(d2arr1[1,])){
      d2arr1.as.d1[i]=mean(d2arr1[,i]) 
    }
  }
  if (method=='rows'){
    d2arr0.as.d1 <- NaN
    for (i in 1:length(d2arr0[,1])){
      d2arr0.as.d1[i]=mean(d2arr0[i,]) 
    }
    d2arr1.as.d1 <- NaN
    for (i in 1:length(d2arr1[,1])){
      d2arr1.as.d1[i]=mean(d2arr1[i,]) 
    }
  }
  
  #now simply apply d1nat
  out <- d1nat(sample0=d2arr0.as.d1, sample1=d2arr1.as.d1, band=band, brks=brks)
  
  return(out)
}

d2spectrum <- function(d2arr0, d2arr1, band=c(0,0), brks=0, 
                       method='rad', meansub=TRUE) {
  
  # Date: 09-02-2015
  
  # check if arrays are squared
  L0 <- c(length(d2arr0[1,]),length(d2arr0[,1])) 
  L1<- c(length(d2arr1[1,]),length(d2arr1[,1])) 
  if ((max(L0)-min(L0)) > 0) {cat('Caution! Array0 isnt square. Will be cut!\n')}
  if ((max(L1)-min(L1)) > 0) {cat('Caution! Array1 isnt square. Will be cut!\n')}
  
  # cut to the same size/square
  L0 = min(L0)
  L1 = min(L1)
  d2arr0 <- d2arr0[1:L0,1:L0]
  d2arr1 <- d2arr1[1:L1,1:L1]
  
  #mean subtraction
  if (meansub==TRUE){
    d2arr0 = d2arr0-mean(d2arr0)
    d2arr1 = d2arr1-mean(d2arr1)
  } 
  
  #calc the size of future fft quarter space
  pw_size0 <- floor(L0/2-1)   
  pw_size1 <- floor(L1/2-1) 
  
  #spectr_calc
  pw0 <- abs(fft(d2arr0))
  pleft <- pw0[1:(pw_size0),(pw_size0+1):(1)]
  pright <- pw0[1:(pw_size0),L0:(L0-pw_size0+1)]
  pw0 = cbind(pleft,pright)
  pw0=pw0^2
  #
  pw1 <- abs(fft(d2arr1))
  pleft <- pw1[1:(pw_size1),(pw_size1+1):(1)]
  pright <- pw1[1:(pw_size1),L1:(L1-pw_size1+1)]
  pw1 = cbind(pleft,pright)
  pw1=pw1^2
  
  #angular pixel gouping
  if (method == 'ang' | method == 'ang90') {
    quanglos0 <- array(0,c(pw_size0, (2*pw_size0+1)))
    for (i in 1:pw_size0){
      for (j in 1:(2*pw_size0+1)) {
        quanglos0[i,j]=(180/pi)*atan((j-pw_size0-1)/(i-0.999999999999999)) 
      }
    }
    #
    quanglos1 <- array(0,c(pw_size1, (2*pw_size1+1)))
    for (i in 1:pw_size1){
      for (j in 1:(2*pw_size1+1)) {
        quanglos1[i,j]=(180/pi)*atan((j-pw_size1-1)/(i-0.999999999999999)) 
      }
    }
    
    #3 arrays to define with angular space
    angle <- c(1:180)
    freq0 <- c(0:180) 
    freq1 <- c(0:180)
    
    #angular powers calc
    pwR0<-rep(0,180); pwR1<-rep(0,180)
    for (k in 1:180) {
      pwR0[k]=sum(pw0[which(angle[k]==(91+as.integer(quanglos0)))])
      pwR1[k]=sum(pw1[which(angle[k]==(91+as.integer(quanglos1)))])
    }
    
    # set the window & adjust data to it:
    if (band[2]==band[1] & band[2]==0) {
      band = c(0, 180)  #freq band by default: full
    }
    if (min(band[1],band[2])<0) {stop("Wrong band!")}  
    if (max(band[1],band[2])>180) {stop("Wrong band!")}  
    pwR0 = pwR0[which(freq0>=band[1] & freq0<=band[2])]
    freq0 = freq0[which(freq0>=band[1] & freq0<=band[2])]
    pwR1  = pwR1[which(freq1>=band[1] & freq1<=band[2])]
    freq1  = freq1[which(freq1>=band[1] & freq1<=band[2])]
    if (length(pwR0)>1 & length(pwR1)>1) {
      full_range = band[2]-band[1]
    } else stop("Band size doesn't match!")
    
    # set the breaks:
    if (brks == 0) {
      num_breaks = full_range #it is 0:180 degrees
      break_size = full_range/num_breaks
      num_breaks = floor(full_range/break_size)
    } else {num_breaks = brks}
    break_size = full_range/num_breaks
    the_breaks = band[1]+c(0:num_breaks)*break_size
    ####
    
    # radial method grouping
  } else {
    method='rad'
    quadrantos0 <- array(0,c(pw_size0,(2*pw_size0+1)))
    for (i in 1:pw_size0) {
      for (j in 1:(2*pw_size0+1)) {
        quadrantos0[i,j]=(i^2+(j-pw_size0-1)^2) 
      }
    }
    #
    quadrantos1 <- array(1,c(pw_size1,(2*pw_size1+1)))
    for (i in 1:pw_size1) {
      for (j in 1:(2*pw_size1+1)) {
        quadrantos1[i,j]=(i^2+(j-pw_size1-1)^2) 
      }
    }
    
    #4 arrays to define with radial space
    radius0 <- c(0:pw_size0)
    radius1 <- c(0:pw_size1)
    freq0 <- 0.5*c(0:(pw_size0-1))/pw_size0
    freq1 <- 0.5*c(0:(pw_size1-1))/pw_size1
    
    #radial power grouped by radii
    pwR0<-rep(0,pw_size0)
    for (k in 1:pw_size0) {
      pwR0[k] = sum(pw0[which(quadrantos0<((radius0[k+1])^2) 
                              & quadrantos0>=((radius0[k])^2))])
    }
    #
    pwR1<-rep(0,pw_size1)
    for (k in 1:pw_size1) {
      pwR1[k] = sum(pw1[which(quadrantos1<((radius1[k+1])^2) 
                              & quadrantos1>=((radius1[k])^2))])
    }
    
    # set the window & adjust data to it:
    if (band[2]==band[1] & band[2]==0) {
      band = c(0, 0.5) #freq band by default: full
    }
    if (min(band[1],band[2])<0) {stop("Wrong band!")}  
    if (max(band[1],band[2])>0.5) {stop("Wrong band!")}  
    pwR0 = pwR0[which(freq0>=band[1] & freq0<=band[2])]
    freq0 = freq0[which(freq0>=band[1] & freq0<=band[2])]
    pwR1  = pwR1[which(freq1>=band[1] & freq1<=band[2])]
    freq1  = freq1[which(freq1>=band[1] & freq1<=band[2])]
    if (length(pwR0)>1 & length(pwR1)>1) {
      full_range = band[2]-band[1]
    } else stop("Band size doesn't match!")
    
    # set the breaks:
    if (brks == 0) {
      num_breaks = min(length(pwR0),length(pwR1)) # dance from what is less size!
      break_size = full_range/num_breaks
      num_breaks = floor(full_range/break_size)
    } else {num_breaks = brks}
    break_size = full_range/num_breaks
    the_breaks = band[1]+c(0:num_breaks)*break_size
  }
  
  #distribution midpoints
  the_mids <- band[1]-break_size/2+c(1:num_breaks)*break_size
  
  #quazi-histogram grouping to form final f0 & f1 distributions
  f0 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f0[i] = sum(pwR0[which(freq0>=the_breaks[i] & freq0<the_breaks[i+1])])
  }
  #f0[num_breaks]=f0[num_breaks]+sum(pwR0[which(
  #freq0==the_breaks[num_breaks+1])]) #take also zeros missed previously
  f1 <- rep(0,num_breaks)
  for (i in 1:(num_breaks)) {
    f1[i] = sum(pwR1[which(freq1>=the_breaks[i] & freq1<the_breaks[i+1])])
  }
  #special case : angle is shifted to 90degrees
  if (method=='ang90') {
    the_mids = the_mids-90 
    for (i in 1:num_breaks) {
      if (the_mids[i]<0) {the_mids[i]=180+the_mids[i]}
    }
  }
  
  #
  if (sum(f0)==0 | sum(f1)==0) {
    stop('One of the samples is a constant. Analysis is not relevant!')
  }
  
  f0=f0/sum(f0)
  f1=f1/sum(f1)
  
  # absolute chaotic entropy calc
  fneutral <- rep(1,length(f0))/length(f0)
  entropy_max <- -sum(log(fneutral)*fneutral)
  
  #statistic summary for new distributions
  f0_n_expected = sum(f0*the_mids)
  f0_n_stdev = sqrt(sum(f0*the_mids^2)-f0_n_expected^2)
  f0_n_entropy = 0.0001*round(-sum(f0*log(f0+1e-15))*10000)
  if (f0_n_entropy<0) {f0_n_entropy=0} #protection patch
  f0_n_1stmod = min(the_mids[which(f0==max(f0))]) #use 1st if where is > 1 modes
  ftemp<-f0; ftemp[which(the_mids==f0_n_1stmod)]=0 #temp without 1st mode 
  f0_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f0_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f0_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  #
  f1_n_expected = sum(f1*the_mids)
  f1_n_stdev = sqrt(sum(f1*the_mids^2)-f1_n_expected^2)
  f1_n_entropy = 0.0001*round(-sum(f1*log(f1+1e-15))*10000)
  if (f1_n_entropy<0) {f1_n_entropy=0} #protection patch
  f1_n_1stmod = min(the_mids[which(f1==max(f1))]) #use 1st if where is > 1 modes
  ftemp<-f1; ftemp[which(the_mids==f1_n_1stmod)]=0 #temp without 1st mode 
  f1_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f1_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f1_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  
  #
  stat_summary <- cbind(expctd=c(f0=f0_n_expected,f1=f1_n_expected),
                        var=c(f0=f0_n_stdev^2,f1=f1_n_stdev^2), 
                        fsum=c(f0=sum(f0),f1=sum(f1)), 
                        xmin=c(f0=min(the_breaks),f1=min(the_breaks)), 
                        xmax=c(f0=max(the_breaks),f1=max(the_breaks)), 
                        n=c(f0=num_breaks,f1=num_breaks),
                        mod1=c(f0=f0_n_1stmod,f1=f1_n_1stmod),
                        mod2=c(f0=f0_n_2ndmod,f1=f1_n_2ndmod),
                        mod3=c(f0=f0_n_3rdmod,f1=f1_n_3rdmod),
                        H_val=c(f0=f0_n_entropy,f1=f1_n_entropy),
                        H_max=c(f0=entropy_max,f1=entropy_max))
  
  #many happy returns!
  return(list(f0=f0, f1=f1, efs=the_mids, method=method,
              stat_summary=stat_summary, essence=c('Spectral distributions')))
}

d2spec<- function(d2arr0, d2arr1, band=c(0,0), brks=0, 
                  method='rad', meansub=TRUE) UseMethod("d2spec")
d2spec.default <- function(d2arr0, d2arr1, band=c(0,0), 
                           brks=0, method='rad', meansub=TRUE) {
  
  # check they are matrices
  if (class(d2arr0)!='matrix' | class(d2arr1)!='matrix') {
    stop('!one or both arrays are not of matrix class!')  
  }
  d2arr0 <- as.matrix(d2arr0)
  d2arr1 <- as.matrix(d2arr1)
  band <- as.vector(band)
  brks <- as.numeric(brks)
  method <- as.character(method)
  meansub <- as.logical(meansub)
  est <- d2spectrum(d2arr0, d2arr1, band, brks, method, meansub)
  est$f0 <- as.vector(est$f0)
  est$f1 <- as.vector(est$f1)
  est$efs <- as.vector(est$efs)
  est$method <- as.character(est$method)
  est$stat_summary <- as.matrix(est$stat_summary)
  est$essence <- as.character(est$essence)
  est$call <- match.call()
  class(est) <- "d2spec"
  est
}

print.d2spec <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Two spectral probability mass functions ($f0,$f1) have ")
  cat("been generated at the common scale of ")
  if (x$method=='ang90' | x$method=='ang') {
    cat("spectral angles ($efs)\n")
  } else {cat("radial freqs ($efs)\n")}
  cat("Statistics summary:\n")
  print(x$stat_summary)
  plot(x)
}

plot.d2spec <- function(x,...){
  min_val <- x$stat_summary[1,4]
  max_val <- x$stat_summary[1,5]
  the_mids <- x$efs
  the_mids_1 <- the_mids+(max_val-min_val)/(x$stat_summary[1,6])/7
  
  if (x$method=='ang90') {xlab <- "Spectr_angles-90, degrees"}  
  if (x$method=='ang') {xlab <- "Spectr_angles, degrees"}
  if (x$method=='rad') {xlab <- "Radial_frequency, 1/pixel"}
  plot(-1,-1, xlim = c(min_val,max_val), 
       ylim = c(0,(1.4*max(x$f0,x$f1))), ylab="Probability Mass", 
       xlab=xlab)
  lines(the_mids, x$f0, type='h', lwd = 3, col=1)
  lines(the_mids_1, x$f1, type='h', lwd = 3, col=8)
  legend(mean(the_mids),(1.4*max(x$f0,x$f1)), c("$f0 to d2arr0","$f1 to d2arr1"),
         lty=c(1,1), col=c(1,8))
  title(x$essence)
}

pvaligner <- function(distribution0, distribution1, 
                      outcomes0=c(0), outcomes1=c(0), tolerance=0) {
  
  # Date: 04-02-2015
  
  f0 <- distribution0
  f1 <- distribution1
  # check the consistence of the things 
  if (min(f0)<0 | min(f1)<0) {
    stop('!One of the distributions is invalid (<0?)')
  }
  if (tolerance == 0) {tolerance = 0.001} #default tolerance <-0.001
  if (tolerance >= 1 | tolerance <= 0) {
    stop('!Tolerance must be: 0<|tolerance|<1 !')
  }  
  
  # set f0,f1 and new scale
  if (length(outcomes0)<=1) {
    x0 <- seq(1:length(f0))
  } else {x0 <- outcomes0}
  if (length(outcomes1)<=1) {
    x1 <- seq(1:length(f1))
  } else {x1 <- outcomes1}
  
  #check other things
  if ((abs(length(x0)-length(f0))>0) | (abs(length(x1)-length(f1))>0)) {
    stop('!Lengths of distribution vector(s) dont match the length(s)
         of vector(s) of outcomes')
  }
  
  #create equidistant integer scale of values
  absolute_min <- min(x0,x1)   
  if (min(x0,x1) <= 0) { 
    x0=x0-absolute_min+1 #shift values all to positives
    x1=x1-absolute_min+1 #shift values all to positives
  }
  
  #check if scales of value arrayas are equidistant and compatible
  #using the tolerance value
  if ((x0[1]>x0[2]) | (x1[1]>x1[2])) {
    stop('Vector of outcomes is in decreasing order! Re-sort please')
  }
  dif_array1 <- NaN
  for (i in 1:length(x0)-1) {dif_array1[i]=x0[i+1]-x0[i]}
  ref_value1 <- (max(dif_array1)-min(dif_array1))/mean(dif_array1)
  if (ref_value1>tolerance) {
    stop('Scale of 1st vector of outcomes is not equidistant!')
  } 
  dif_array2 <- NaN
  for (i in 1:length(x1)-1) {dif_array2[i]=x1[i+1]-x1[i]}
  ref_value2 <- (max(dif_array2)-min(dif_array2))/mean(dif_array2)
  if (ref_value2>tolerance) {
    stop('Scale of 2nd vector of outcomes is not equidistant!')
  } 
  if (abs(mean(dif_array2)-mean(dif_array1))>tolerance*mean(dif_array1)) {
    stop('Vectors of outcomes are not compatible!')
  }
  
  #create the common scale of values
  xvector <- c(1:(1+as.integer(round((max(x0,x1)-min(x0,x1))/((x0[2]-x0[1]))))))
  
  #new distributions at common scale  
  if (min(x0)<min(x1)) {
    temp_array <- rep(0,length(xvector)-length(x0))
    f0 <- c(f0,temp_array)
    #f0[which(f0<1e-15)]=0
    temp_array <- rep(0,length(xvector)-length(x1))
    f1 <- c(temp_array,f1)
    
  } else {
    temp_array <- rep(0,length(xvector)-length(x1))
    f1<- c(f1,temp_array)
    temp_array <- rep(0,length(xvector)-length(x0))
    f0 <- c(temp_array,f0)
  }
  #
  f0=f0/sum(f0)
  f1=f1/sum(f1)
  
  #common x scale:
  x0=x0-min(x0,x1)+1
  the_mids=xvector+absolute_min-1
  
  # absolute chaotic entropy calc
  fneutral <- rep(1,length(f0))/length(f0)
  entropy_max <- -sum(log(fneutral)*fneutral)
  
  #statistic summary for new distributions
  f0_n_expected = sum(f0*the_mids)
  f0_n_var = (sum(f0*the_mids^2)-f0_n_expected^2)
  f0_n_entropy = -sum(f0*log(f0+1e-15))
  if (f0_n_entropy<0) {f0_n_entropy=0} #protection patch
  f0_n_1stmod = min(the_mids[which(f0==max(f0))]) #use 1st if where is > 1 modes
  ftemp<-f0; ftemp[which(the_mids==f0_n_1stmod)]=0 #temp without 1st mode 
  f0_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f0_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f0_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  #
  f1_n_expected = sum(f1*the_mids)
  f1_n_var = (sum(f1*the_mids^2)-f1_n_expected^2)
  f1_n_entropy = -sum(f1*log(f1+1e-15))
  if (f1_n_entropy<0) {f1_n_entropy=0} #protection patch
  f1_n_1stmod = min(the_mids[which(f1==max(f1))]) #use 1st if where is > 1 modes
  ftemp<-f1; ftemp[which(the_mids==f1_n_1stmod)]=0 #temp without 1st mode 
  f1_n_2ndmod = min(the_mids[which(ftemp==max(ftemp))])
  ftemp[which(the_mids==f1_n_2ndmod)]=0 #temp without 1st & 2nd modes 
  f1_n_3rdmod = min(the_mids[which(ftemp==max(ftemp))])
  
  stat_summary <- cbind(expctd=c(f0=f0_n_expected,f1=f1_n_expected),
                        var=c(f0=f0_n_var,f1=f1_n_var), 
                        fsum=c(f0=sum(f0),f1=sum(f1)), 
                        xmin=c(f0=min(the_mids),f1=min(the_mids)), 
                        xmax=c(f0=max(the_mids),f1=max(the_mids)), 
                        n=c(f0=length(the_mids),f1=length(the_mids)),
                        mod1=c(f0=f0_n_1stmod,f1=f1_n_1stmod),
                        mod2=c(f0=f0_n_2ndmod,f1=f1_n_2ndmod),
                        mod3=c(f0=f0_n_3rdmod,f1=f1_n_3rdmod),
                        H_val=c(f0=f0_n_entropy,f1=f1_n_entropy),
                        H_max=c(f0=entropy_max,f1=entropy_max))
  
  #many happy returns!
  return(list(new_scale=the_mids, f0=f0, f1=f1, stat_summary=stat_summary))
  }


pvalign<- function(distribution0, distribution1, 
                   outcomes0=c(0), outcomes1=c(0), tolerance=0) UseMethod("pvalign")
pvalign.default <- function(distribution0, distribution1, 
                            outcomes0=c(0), outcomes1=c(0), tolerance=0) {
  
  distribution0 <- as.vector(distribution0)
  distribution1 <- as.vector(distribution1)
  outcomes0 <- as.vector(outcomes0)
  outcomes1 <- as.vector(outcomes1)
  tolerance <- as.numeric(tolerance)
  est <- pvaligner(distribution0, distribution1, 
                   outcomes0, outcomes1, tolerance)
  est$new_scale <- as.vector(est$new_scale)
  est$f0 <- as.vector(est$f0)
  est$f1 <- as.vector(est$f1)
  est$stat_summary <- as.matrix(est$stat_summary)
  est$call <- match.call()
  class(est) <- "pvalign"
  est
}

print.pvalign <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Two spectral probability mass functions ($f0,$f1) have ")
  cat("been generated at the new common scale ($new_scale)\n")
  cat("Statistics summary:\n")
  print(x$stat_summary)
  plot(x)
}

plot.pvalign <- function(x,...){
  
  min_val <- x$stat_summary[1,4]
  max_val <- x$stat_summary[1,5]
  the_mids <- x$new_scale
  the_mids_1 <- the_mids+(max_val-min_val)/(x$stat_summary[1,6])/7
  
  plot(-1,-1, xlim = c(min_val,max_val), 
       ylim = c(0,(1.4*max(x$f0,x$f1))), ylab="Probability Mass", 
       xlab = "New scale of outcomes")
  lines(the_mids, x$f0, type='h', lwd = 3, col=1)
  lines(the_mids_1, x$f1, type='h', lwd = 3, col=8)
  legend(mean(the_mids),(1.4*max(x$f0,x$f1)), 
         c("f0 to outcomes0","f1 to outcomes1"),
         lty=c(1,1), col=c(1,8))
  title("Re-scaled distributions")
}

utild1bin <- function(arr0, arr1, method='mean', trsh=0, inverted=FALSE, d2=FALSE){
  
  # Date: 12-02-2015
  
  thr0<-mean(c(arr0,arr1)) #default thresholds is mean for two
  thr1<-thr0
  
  if (method == 'mean') {thr1=mean(c(arr0,arr1));thr0=thr1} #default
  if (method == 'mean0') {thr1=mean(c(arr0));thr0=thr1}
  if (method == 'mean1') {thr1=mean(c(arr1));thr0=thr1}
  if (method == 'med')  {thr1=median(c(arr0,arr1));thr0=thr1} 
  if (method == 'med0')  {thr1=median(c(arr0));thr0=thr1} 
  if (method == 'med1')  {thr1=median(c(arr1));thr0=thr1} 
  if (method == 'spec') {
    if (length(trsh)==2) {thr0=trsh[1]; thr1=trsh[2]}
    if (length(trsh)==1) {thr0=trsh; thr1=trsh}
    if (length(trsh)>2) {stop('!more than 2 threshold levels introduced')}
  }
  if (method == 'imean') {thr0=mean(c(arr0)); thr1=mean(c(arr1))}
  if (method == 'imed') {thr0=median(c(arr0)); thr1=median(c(arr1))}
  if (method == 'i1stQ') {
    thr0=summary(c(arr0))[2]
    thr1=summary(c(arr1))[2]
  }
  if (method == 'i3rdQ') {
    thr0=summary(c(arr0))[5]
    thr1=summary(c(arr1))[5]
  }
  if (method == 'isubsd') {
    thr0=mean(c(arr0))-sd(c(arr0))
    thr1=mean(c(arr1))-sd(c(arr1))
  }
  if (method == 'iaddsd') {
    thr0=mean(c(arr0))+sd(c(arr0))
    thr1=mean(c(arr1))+sd(c(arr1))
  }
  
  if (d2==FALSE) {
    b0 <- NaN
    b1 <- NaN
    if (class(arr0)=='matrix' | class(arr1)=='matrix') {
      stop ('One of the arrays has > 1 dimensions')
    }
  } else {
    b0 <- array(0,dim(arr0))
    b1 <- array(0,dim(arr1))
  }
  
  if (inverted) {
    b0[which(arr0<=thr0)]=TRUE
    b0[which(arr0>thr0)]=FALSE
    b1[which(arr1<=thr1)]=TRUE
    b1[which(arr1>thr1)]=FALSE
  } else {
    b0[which(arr0<=thr0)]=FALSE
    b0[which(arr0>thr0)]=TRUE
    b1[which(arr1<=thr1)]=FALSE
    b1[which(arr1>thr1)]=TRUE
  }
  
  #counts:
  ntrue0 = sum(b0[which(b0==TRUE)])
  ntrue1 = sum(b1[which(b1==TRUE)])
  nfalse0 = length(b0)-ntrue0
  nfalse1 = length(b1)-ntrue1
  
  #
  counts<-cbind(n_false=c(bin0=nfalse0,bin1=nfalse1),
                n_true=c(bin0=ntrue0, bin1=ntrue1))
  
  return(list(bin0=b0,bin1=b1,counts=counts))
}


########
utild2bin <- function(d2arr0, d2arr1, method='mean', trsh=0, inverted=FALSE){
  
  if (class(d2arr0)!='matrix' | class(d2arr1)!='matrix') {
    stop ('One of the arrays has < 2 dimensions')
  }
  
  return(utild1bin(d2arr0, d2arr1, 
                   method=method, trsh=trsh, inverted=inverted, d2=TRUE))
}

utild1clean <- function(arr0, arr1, method='bothends', nsigma=3, d2=FALSE) {
  
  # Date 15-02-2015 
  
  M0=mean(arr0); M1=mean(arr1)
  S0=sd(arr0); S1=sd(arr1)
  
  left0 <- M0-nsigma*S0; right0 <- M0+nsigma*S0
  left1 <- M1-nsigma*S1; right1 <- M1+nsigma*S1
  
  if (d2==FALSE) {
    c0 <- arr0
    c1 <- arr1
    if (class(arr0)=='matrix' | class(arr1)=='matrix') {
      stop ('One of the arrays has > 1 dimensions')
    }
  } else {
    c0 <- arr0
    c1 <- arr1
  }
  
  Ncount0 <- 0; Ncount1 <- 0
  
  if (method=='left') {
    Ncount0=length(c0[which(c0<left0)])
    c0[which(arr0<left0)]=left0
    Ncount1=length(c1[which(c1<left1)])
    c1[which(arr1<left1)]=left1
  } 
  if (method=='right') {
    Ncount0=length(c0[which(c0>right0)])
    c0[which(arr0>right0)]=right0
    Ncount1=length(c1[which(c1>right1)])
    c1[which(arr1>right1)]=right1
  }
  if (method!='left' & method!='right') {
    Ncount0=length(c0[which(c0<left0)])+length(c0[which(c0>right0)])
    c0[which(arr0<left0)]=left0
    c0[which(arr0>right0)]=right0
    Ncount1=length(c1[which(c1<left1)])+length(c1[which(c1>right1)])
    c1[which(arr1<left1)]=left1
    c1[which(arr1>right1)]=right1
  }
  
  cat(as.character(Ncount0)); cat(" outliers in array0 and ");
  cat(as.character(Ncount1)); cat(" outliers in array1");
  cat(" have been replaced by critical values\n"); cat("\n")
  
  return(list(clean0=c0, clean1=c1))
}

utild2clean <- function(d2arr0, d2arr1, method='bothends', nsigma=3){
  
  if (class(d2arr0)!='matrix' | class(d2arr1)!='matrix') {
    stop ('One of the arrays has < 2 dimensions')
  }
  
  return(utild1clean(d2arr0, d2arr1, method=method, nsigma=nsigma, d2=TRUE))
}

utild1filt <- function(arr0, arr1, outsize=2, strong=1) {
  
  # Date 23-02-2015
  
  if (outsize<0) {stop('! bandpass frame is invalid (must be >= 0)')}
  if (outsize>0.4*min(length(arr0),length(arr1))) {
    stop('! outsize is too high')
  }
  if (strong<=0) {strong=1} #stupid protection
  if (class(arr0)=="matrix" | class(arr1)=="matrix") {
    stop('! one/both arrays have >1 dims!')
  }
  
  #gauss kernel calc func
  d1gausscoeffs <- function(s, n) {
    t <- seq(-n,n,1)
    out <- exp(-(t^2/(2*s^2))); out=out/sum(out)
    return (out)
  }
  
  #convolution for lowpass
  sigma <- ceiling(outsize)*strong
  radius <- ceiling(outsize)
  if (outsize==0) {sigma=1}
  kern <- d1gausscoeffs(sigma, radius)
  outlow0 <- arr0
  outlow1 <- arr1
  if (outsize==0) {radius=1}
  outlow0[(radius+1):(length(arr0)-radius)] = 
    convolve(arr0,kern,type='f')
  outlow1[(radius+1):(length(arr1)-radius)] = 
    convolve(arr1,kern,type='f')
  outlow0[1:(radius)] = outlow0[radius+1] #border values left
  outlow0[(length(arr0)-radius+1):length(arr0)] = outlow0[(length(arr0)-radius)] # also border 
  outlow1[1:(radius)] = outlow1[radius+1] # also border
  outlow1[(length(arr1)-radius+1):length(arr1)] = outlow1[(length(arr1)-radius)] # also border 
  
  
  out0<-outlow0
  out1<-outlow1
  
  return(list(filt0=out0, filt1=out1))
}


#####

utild2filt <- function(d2arr0, d2arr1, outsize=2, strong=1) {
  
  # Date 23-02-2015
  
  if (class(d2arr0)!="matrix" | class(d2arr1)!="matrix") {
    stop('! one/both arrays have <2 dims!')
  }
  if (outsize<0) {stop('! bandpass frame is invalid (must be >= 0)')}
  if (outsize>0.4*min(length(d2arr0[,1]),length(d2arr1[,1]),
                      length(d2arr0[1,]),length(d2arr1[1,]))) {
    stop('! outsize is too high')
  }
  if (strong<=0) {strong=1} #stupid protection
  
  #gauss kernel calc func
  # 2D
  d2gausscoeffs <- function(s, n) {
    t<-array(0, c((2*n+1),(2*n+1)))
    for (i in 1:(2*n+1)) {
      for (j in 1:(2*n+1)) {
        t[i,j] <- sqrt((i-n-1)^2+(j-n-1)^2)
      }
    }
    out <- exp(-(t^2/(2*s^2))); out=out/sum(out)
    return (out)
  }
  
  # convolve 2d function
  conv2d <- function(arr2d, kernel2) {
    gap <- (length(kernel2[1,])-1)/2
    Wn <- length(arr2d[1,])-2*gap
    Hn <- length(arr2d[,1])-2*gap
    
    out <- arr2d
    for (i in (1+gap):(Hn+gap)) {
      for (j in (1+gap):(Wn+gap)) {
        out[i,j]=0
        for (k in -gap:gap) {
          for (l in -gap:gap) {
            out[i,j]=out[i,j]+arr2d[i+k,j+l]*kernel2[k+gap+1,l+gap+1] 
          }
        }
      }
    }
    return(out)
  }
  
  #convolution for lowpass
  sigma <- ceiling(outsize)*strong
  radius <- ceiling(outsize)
  if (outsize==0) {sigma=1}
  kern <- d2gausscoeffs(sigma, radius)
  if (outsize==0) {radius=1}
  outlow0 = conv2d(d2arr0,kern)
  outlow1 = conv2d(d2arr1,kern)
  
  out0<-outlow0
  out1<-outlow1
  
  return(list(filt0=out0, filt1=out1))
  
}

utild1group <- function(arr0, arr1, radius=1, method='split1') {
  
  #Date 24-02-15
  
  if (radius<=0) {stop('!radius is invalid (must be > 0)!')}
  if (radius>0.4*min(length(arr0),length(arr1))) {
    stop('! radius is too high!')
  }
  
  groupi=floor(2*radius+1)
  
  if (class(arr0)=="matrix" | class(arr1)=="matrix") {
    stop('! one/both arrays have >1 dims!')
  }
  
  
  if (method!='splitN') {
    
    new_size0 <- floor(length(arr0)/groupi)
    new_size1 <- floor(length(arr1)/groupi)
    s0 <- c(1:new_size0)*0
    s1 <- c(1:new_size1)*0
    
    for (i in 1:new_size0) {
      s0[i] = mean(arr0[((i-1)*groupi+1):(i*groupi)])
    }
    for (i in 1:new_size1) {
      s1[i] = mean(arr1[((i-1)*groupi+1):(i*groupi)])
    }
    
  } 
  
  if (method=='splitN') {
    
    new_size0 <- length(arr0)-groupi+1
    new_size1 <- length(arr1)-groupi+1
    s0 <- c(1:new_size0)*0
    s1 <- c(1:new_size1)*0
    
    for (i in 1:new_size0) {
      s0temp <- c(1:groupi)*0
      for (k in 1:groupi) {
        #j = ((k-1)*(new_size0/radius)+i) # index of a new array
        #shift = k+i-1
        s0temp[k] = arr0[i+k-1]
      }
      s0[i] = mean(s0temp)
    }
    
    for (i in 1:new_size1) {
      s1temp <- c(1:groupi)*0
      for (k in 1:groupi) {
        #j = ((k-1)*(new_size0/radius)+i) # index of a new array
        #shift = k+i-1
        s1temp[k] = arr1[i+k-1]
      }
      s1[i] = mean(s1temp)
    }
    
  }
  
  
  out <- list(group0=s0, group1=s1)
  return(out)
  
}


######
utild2group <- function(d2arr0, d2arr1, radius=1, method='split1') {
  
  #Date 24-02-15
  
  if (class(d2arr0)!="matrix" | class(d2arr1)!="matrix") {
    stop('! one/both arrays have <2 dims!')
  }
  if (radius<=0) {stop('!radius is invalid (must be > 0)!')}
  if (radius>0.4*min(length(d2arr0[1,]),length(d2arr1[1,]),
                     length(d2arr0[,1]),length(d2arr1[,1]))) {
    stop('! radius is too high!')
  }
  
  groupi=floor(2*radius+1)
  
  if (method!='splitN') {
    
    new_width0 <- floor(length(d2arr0[,1])/groupi)
    new_width1 <- floor(length(d2arr1[,1])/groupi)
    new_heigth0 <- floor(length(d2arr0[1,])/groupi)
    new_heigth1 <- floor(length(d2arr1[1,])/groupi)
    
    s0 <- array(0, c(new_width0, new_heigth0))
    s1 <- array(0, c(new_width1, new_heigth1))
    
    for (i in 1:new_width0) {
      for (j in 1:new_heigth0) {
        s0[i,j] = mean(d2arr0[((i-1)*groupi+1):(i*groupi),
                              ((j-1)*groupi+1):(j*groupi)])
      }
    }
    for (i in 1:new_width1) {
      for (j in 1:new_heigth1) {
        s1[i,j] = mean(d2arr1[((i-1)*groupi+1):(i*groupi),
                              ((j-1)*groupi+1):(j*groupi)])
      }
    }
    
  } 
  
  if (method=='splitN') {
    
    new_width0 <- length(d2arr0[,1])-groupi+1
    new_width1 <- length(d2arr1[,1])-groupi+1
    new_heigth0 <- length(d2arr0[1,])-groupi+1
    new_heigth1 <- length(d2arr1[1,])-groupi+1
    
    s0 <- array(0, c(new_width0, new_heigth0))
    s1 <- array(0, c(new_width1, new_heigth1))
    
    for (i in 1:new_width0) {
      for (j in 1:new_heigth0) {
        
        s0temp <- array(0, c(groupi,groupi))
        
        for (k in 1:groupi) {
          for (m in 1:groupi) {
            s0temp[k,m] = d2arr0[i+k-1,j+m-1]
          }
        }
        
        s0[i,j] = mean(s0temp)
        
      }
    }
    
    for (i in 1:new_width1) {
      for (j in 1:new_heigth1) {
        
        s1temp <- array(0, c(groupi,groupi))
        
        for (k in 1:groupi) {
          for (m in 1:groupi) {
            s1temp[k,m] = d2arr1[i+k-1,j+m-1]
          }
        }
        s1[i,j] = mean(s1temp)
      }
    }
    
  }
  
  
  out <- list(group0=s0, group1=s1)
  return(out)
  
}

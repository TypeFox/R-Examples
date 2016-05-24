#' @keywords internal
deconvolve_inputs <- function(inputs, sep, header, out, control, njobs) {
  ## Helper function to deconvolve all inputs prior to running GIMME
  #require(parallel) #should be moved to import in NAMESPACE for R package
  #require(foreach)
  #require(doSNOW)
  #require(doParallel)
  
  if (missing(njobs)) {
    njobs <- parallel::detectCores(logical=TRUE) #on hyperthreading CPUs, use all threads (typically 2/core)
    message("In deconvolve_inputs, njobs not specified. Defaulting to number of cores: ", njobs, ".")
  }
  
  if (njobs > 1) {
    setDefaultClusterOptions(master="localhost")
    clusterobj <- makeSOCKcluster(njobs)
    registerDoSNOW(clusterobj)
    #clusterobj <- makeCluster(njobs)
    #registerDoParallel(clusterobj)
    
    on.exit(try(stopCluster(clusterobj)))    
  }
  
  stopifnot(is.numeric(control$TR) && control$TR > 0) #user must specify TR of data for deconvolution to be well specified
  if (!exists('method', where=control)) {
    message("Defaulting to Bush and Cisler 2013 deconvolution")
    control$method  <- "bush"
  }
  
  if (control$method == "bush") {
    if (!exists('nev_lr', where=control)) { 
      message("nev_lr (learning rate) not specified for Bush deconvolution. Default to .01")
      control$nev_lr  <- .01  #neural events learning rate
    }
    if (!exists('epsilon', where=control)) { 
      message("epsilon not specified for Bush deconvolution. Default to .005")
      control$epsilon <- .005 #convergence criterion
    }
    if (!exists('kernel', where=control)) { 
      message("HRF kernel not specified for Bush deconvolution. Default to spm double gamma")
      control$kernel  <- spm_hrf(control$TR)$hrf #default to canonical SPM difference of gammas with specified TR
    }
  }
  
  if (control$method == "wu") {
    if (!exists('threshold', where=control)) { 
      message("Activation threshold not specified for Wu deconvolution. Default to 1.0 SD")
      control$threshold <- 1.0 #SD
    }
    if (!exists('max_lag', where=control)) { 
      message("Maximum event-to-neural onset lag not specified for Wu deconvolution. Default to 10 seconds.")
      control$max_lag = ceiling(10/control$TR) #10 seconds maximum lag from neural event to HRF onset
    }
  }
  
  ##now proceed on deconvolution
  stopifnot(file.exists(out))
  deconv.dir <- file.path(out, "deconvolved_inputs")
  dir.create(deconv.dir, showWarnings=FALSE)
  
  #deconvolve inputs in parallel (over subjects/inputs)
  #for (f in inputs) {
  #using doSNOW, need to explicitly export functions and subfunctions for deconvolution to each worker
  exportFuncs <- c("deconvolve_nlreg", "sigmoid", "dsigmoid", "generate_feature", "wgr_deconv_canonhrf_par",
      "wgr_adjust_onset", "wgr_trigger_onset", "CanonicalBasisSet", "Fit_Canonical_HRF", "get_parameters2",
      "spm_hrf", "spm_get_bf", "spm_gamma_bf", "spm_orth")
  f = NULL
  files <- foreach(f=inputs, .combine=c, .export=exportFuncs, .inorder=TRUE) %dopar% {
    stopifnot(file.exists(f))
    d <- as.matrix(read.table(f, header=header, sep=sep))
    if (control$method == "wu") {
      d_deconv <- wgr_deconv_canonhrf_par(d, thr=control$threshold, event_lag_max=control$max_lag, TR=control$TR)$data_deconv #already supports multi-variable (column) deconvolution
    } else if (control$method == "bush") {
      #parApply? Already parallel above, so would need to rework this or do the %:% nested parallel
      d_deconv <- apply(d, 2, deconvolve_nlreg, kernel=control$kernel, nev_lr=control$nev_lr, epsilon=control$epsilon)
    }
    
    outfile <- file.path(deconv.dir, paste0(tools::file_path_sans_ext(basename(f)), "_deconvolve.txt"))
    write.table(d_deconv, file=outfile, sep=sep, quote=FALSE, row.names=FALSE)
    
    return(outfile)
  }
  
  return(files)
}



## R versions of deconvolution algorithms for testing in GIMME
## Ported from MATLAB by Michael Hallquist, September 2015

## R port of Bush and Cisler 2013, Magnetic Resonance Imaging
## Adapted from the original provided by Keith Bush

## Author:      Keith Bush, PhD
## Institution: University of Arkansas at Little Rock
## Date:        Aug. 9, 2013

deconvolve_nlreg <- function(BOLDobs, kernel, nev_lr=.01, epsilon=.005) {
  ## Description:
  ## This function deconvolves the BOLD signal using Bush 2011 method
  ##
  ## Inputs:
  ##    BOLDobs - observed BOLD timeseries
  ##    kernel  - assumed kernel of the BOLD signal
  ##    nev_lr  - learning rate for the assignment of neural events
  ##    epsilon - relative error change (termination condition)
  ##
  ## Outputs:
  ##    encoding - reconstructed neural events
  
  ## Determine time series length
  N = length(BOLDobs)
  
  ##Calc simulation steps related to simulation time
  K = length(kernel)
  A = K - 1 + N
  
  ##Termination Params
  preverror = 1e9 #previous error
  currerror = 0   #current error
  
  ##Construct random activation vector (fluctuate slightly around zero between -2e-9 and 2e-9)
  activation = rep(2e-9, A)*runif(A) - 1e-9
  
  ##Presolve activations to fit target_adjust as encoding
  max_hrf_id_adjust = which.max(kernel) - 1 #element of kernel 1 before max
  BOLDobs_adjust = BOLDobs[max_hrf_id_adjust:N]
  pre_encoding = BOLDobs_adjust - min(BOLDobs_adjust)
  pre_encoding = pre_encoding/max(pre_encoding) #unit normalize
  encoding = pre_encoding
  activation[K:(K-1 + length(BOLDobs_adjust))] = log(pre_encoding/(1-pre_encoding))
  
  while (abs(preverror-currerror) > epsilon) {
    
    ##Compute encoding vector
    encoding = sigmoid(activation)
    
    ##Construct feature space
    feature = generate_feature(encoding,K)
    
    ##Generate virtual bold response by multiplying feature (N x K) by kernel (K x 1) to get N x 1 estimated response
    ytilde = feature[K:nrow(feature),] %*% kernel
    
    ##Convert to percent signal change
    meanCurrent = mean(ytilde)
    brf = ytilde - meanCurrent
    brf = brf/meanCurrent
    
    ##Compute dEdbrf
    dEdbrf = brf - BOLDobs
    
    ##Assume normalization does not impact deriv much.
    dEdy = dEdbrf
    
    ##Precompute derivative components
    dEde = diag(K) %*% kernel
    back_error = c(rep(0, K-1), dEdy, rep(0, K-1))
    
    ##Backpropagate Errors
    delta = c()
    for (i in 1:A) {
      active = activation[i]
      deda = dsigmoid(active);
      dEda = dEde * deda
      this_error = back_error[i:(i-1+K)]
      delta = c(delta, sum(dEda * this_error))
    }
    
    ##Update estimate
    activation = activation - nev_lr * delta
    
    ##Iterate Learning
    preverror = currerror
    currerror = sum(dEdbrf^2)
  }
  
  ## remove the initial timepoints corresponding to HRF (so that returned signal matches in time and length)
  encoding <- encoding[K:length(encoding)]
  
  return(encoding)    
}


## Support functions
sigmoid <- function(x) {
  y <- 1/(1+exp(-x))
  return(y)
}

dsigmoid <- function(x) {                        
  y=(1-sigmoid(x))*sigmoid(x)
  return(y)
}

generate_feature <- function(encoding, K) {
  fmatrix = matrix(0, length(encoding), K)
  fmatrix[,1] = encoding
  
  for (i in 2:K) {
    fmatrix[,i] = c(rep(0, i-1), encoding[1:(length(encoding) - (i-1))])
  }
  return(fmatrix)
}


#####
## Wu code

## R port of Wu et al. 2013, Medical Image Analysis
## Adapted from the original provided by Daniele Marinazzo

wgr_deconv_canonhrf_par <- function(data, thr=1.0, event_lag_max, TR) {
  ### function [data_deconv event HRF adjust_global PARA] = wgr_deconv_canonhrf_par(data,thr,event_lag_max,TR)
  
  ### this function implements the method described in
  ### Wu et al,
  ### A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data,
  ### Med Image Anal. 2013 Jan 29. pii: S1361-8415(13)00004-2. doi: 10.1016/j.media.2013.01.003
  
  ### input
  ### data, dimensions time points x number of voxels, normalized
  ### threshold, assuming data are normalized
  ### event_lag_max: maximum time from neural event to BOLD event in bins, not time
  ### (e.g. if we assume 10seconds, and TR=2s,  set the value to 10/2=5)
  ### TR is the TR parameter, in seconds.
  
  ### Some parts of the code (subfunction: Fit_Canonical_HRF, CanonicalBasisSet, get_parameters2) were modified from the hemodynamic response estimation toolbox(http://www.stat.columbia.edu/~martin/HRF_Est_Toolbox.zip).
  ###
  ### the code uses the parallel for loop parfor. In case of older matlab versions, parfor can be changed to standard for.
  ###
  ### The default is using canonical hrf and two derivatives, as described in the paper.
  ### The function can be modified to use instead point processes, FIR model, etc.
  
  ##force data to explicit matrix form (time x variables) for proper looping
  ##this will allow data to be passed in as a 1-D vector for single variable problems
  if (!inherits(data, "matrix")) { data <- matrix(data, ncol=1) }
  
  N = nrow(data); nvar = ncol(data)
  even_new = wgr_trigger_onset(data,thr)
  
  p_m=3 #define what HRF to fit
  ## Options: p_m=1 - only canonical HRF
  ##          p_m=2 - canonical + temporal derivative
  ##          p_m=3 - canonical + time and dispersion derivative
  
  T = round(30/TR) ## assume HRF effect lasts 30s.
  data_deconv = matrix(0, nrow=N, ncol=nvar)
  HRF  = matrix(0, nrow=T, ncol=nvar)
  PARA = matrix(0, nrow=3, ncol=nvar)
  event = list()
  event_lag <- rep(0, nvar)
  ##warning off
  
  ##can parallelize over nvar
  for (i in 1:nvar) {    
    out = wgr_adjust_onset(data[,i], even_new[[i]], event_lag_max, TR, p_m, T, N)
    data_deconv[,i] <- out$data_deconv
    HRF[,i] <- out$hrf
    event[[i]] <- out$events
    event_lag[i] <- out$event_lag
    PARA[,i] <- out$param        
  }
  
  ##warning on
  return(list(data_deconv=data_deconv, event=event, HRF=HRF, event_lag=event_lag, PARA=PARA))
  
}

wgr_adjust_onset <- function(dat,even_new,event_lag_max,TR,p_m,T,N) {  
  ## global adjust.
  kk=1 #this is just a 1-based version of the loop iterator event_lag...
  hrf   = matrix(NA_real_, nrow=T, ncol=event_lag_max+1)
  param = matrix(NA_real_, nrow=p_m, ncol=event_lag_max+1)
  Cov_E = rep(NA_real_, event_lag_max+1)
  for (event_lag in 0:event_lag_max) {
    RR = even_new - event_lag; RR = RR[RR >= 0]
    design = matrix(0, nrow=N, ncol=1)
    design[RR,1] = 1 #add pseudo-events to design matrix
    fit = Fit_Canonical_HRF(dat,TR,design,T,p_m);
    hrf[,kk] <- fit$hrf
    param[,kk] <- fit$param
    Cov_E[kk] <- cov(fit$e) #covariance of residual
    kk = kk+1;
  }
  
  C   = min(Cov_E)
  ind = which.min(Cov_E)
  ad_global = ind - 1 #begin with 0.
  even_new = even_new - ad_global
  even_new = even_new[even_new>=0]
  hrf = hrf[,ind] #keep only best HRF (minimize error of pseudo-event timing)
  param = param[,ind] #keep only best params
  
  ## linear deconvolution.
  H = fft(c(hrf, rep(0, N-T)))  ##    H=fft([hrf; zeros(N-T,1)]);
  M = fft(dat)
  data_deconv = Re(fft(Conj(H)*M/(H*Conj(H)+C), inverse=TRUE)/length(H)) ## only keep real part -- there is a tiny imaginary residue in R
  
  return(list(data_deconv=data_deconv, hrf=hrf, events=even_new, event_lag=ad_global, param=param))
}

wgr_trigger_onset <- function(mat, thr) {
  ##function [oneset] = wgr_trigger_onset(mat,thr)
  N = nrow(mat); nvar = ncol(mat)
  mat = apply(mat, 2, scale) #z-score columns
  oneset <- list()
  ## Computes pseudo event.
  for (i in 1:nvar) {
    oneset_temp = c()
    for (t in 2:(N-1)) {
      if (mat[t,i] > thr && mat[t-1,i] < mat[t,i] && mat[t,i] > mat[t+1,i]) { ## detects threshold
        oneset_temp = c(oneset_temp, t)
      }
    }
    oneset[[i]] = oneset_temp
  }
  
  return(oneset)
}




#####
## Helper functions for Wu deconvolution algorithm
## Original code from Lindquist and Wager HRF Toolbox
CanonicalBasisSet <- function(TR) {
  len = round(30/TR) # 30 secs worth of images
  xBF <- list()
  xBF$dt = TR
  xBF$length= len
  xBF$name = 'hrf (with time and dispersion derivatives)'
  xBF = spm_get_bf(xBF)
  
  v1 = xBF$bf[1:len,1]
  v2 = xBF$bf[1:len,2]
  v3 = xBF$bf[1:len,3]
  
  ## orthogonalize
  h = v1
  dh =  v2 - (v2 %*% v1/norm(v1, "2")^2)*v1
  dh2 =  v3 - (v3 %*% v1/norm(v1, "2")^2)*v1 - (v3 %*% dh/norm(dh, "2")^2)*dh
  
  ## normalize amplitude
  h = h/max(h)
  dh = dh/max(dh)
  dh2 = dh2/max(dh2)
  
  return(list(h=h, dh=dh, dh2=dh2))
}

Fit_Canonical_HRF <- function(tc, TR, Run, T, p) {
  ##function [hrf, e, param] = Fit_Canonical_HRF(tc,TR,Run,T,p)
  ##
  ## Fits GLM using canonical hrf (with option of using time and dispersion derivatives)';
  ##
  ## INPUTS:
  ##
  ## tc    - time course
  ## TR    - time resolution
  ## Runs  - expermental design
  ## T     - length of estimated HRF
  ## p     - Model type
  ##
  ## Options: p=1 - only canonical HRF
  ##          p=2 - canonical + temporal derivative
  ##          p=3 - canonical + time and dispersion derivative
  ##
  ## OUTPUTS:
  ##
  ## hrf   - estimated hemodynamic response function
  ## fit   - estimated time course
  ## e     - residual time course
  ## param - estimated amplitude, height and width
  
  len = length(Run)
  
  X = matrix(0, nrow=len, ncol=p)
  
  bf = CanonicalBasisSet(TR)
  h  = bf$h; dh = bf$dh; dh2 = bf$dh2
  
  v = convolve(Run,rev(h), type="open") #this is the R equivalent of conv(Run, h)
  X[,1] = v[1:len]
  
  if (p > 1) {        
    v = convolve(Run,rev(dh), type="open")
    X[,2] = v[1:len]
  }
  
  if (p > 2) {
    v = convolve(Run,rev(dh2), type="open")
    X[,3] = v[1:len]
  }
  
  X   = cbind(rep(1, len), X) #add intercept
  b   = MASS::ginv(X) %*% tc
  e   = tc - X %*% b
  fit = X %*% b
  
  b = b[2:length(b)]
  
  if (p == 2) {
    bc = sign(b[1])*sqrt(b[1]^2 + b[2]^2)
    H = cbind(h, dh)
  } else if (p == 1) {
    bc = b[1]
    H = matrix(h, ncol=1)
  } else if (p>2) {   
    bc = sign(b[1])*sqrt(b[1]^2 + b[2]^2 + b[3]^2)
    H = cbind(h, dh, dh2)
  }
  
  hrf = H %*% b
  
  param = get_parameters2(hrf,T)
  
  return(list(hrf=hrf, e=e, param=param))
}


get_parameters2 <- function(hdrf, t) {
  
  ##function [param] = get_parameters2(hdrf,t)
  
  ## Find model parameters
  ##
  ## Height - h
  ## Time to peak - p (in time units of TR seconds)
  ## Width (at half peak) - w
  
  ## Calculate Heights and Time to peak:
  
  ## n = t(end)*0.6;
  n = round(t*0.8)
  
  p = which.max(abs(hdrf[1:n]))
  h = hdrf[p]
  
  ##if (p > t(end)*0.6), warning('Late time to peak'), end;
  
  if (h > 0) {
    v = as.numeric(hdrf >= h/2)
  } else {
    v = as.numeric(hdrf <= h/2)
  }
  
  b = which.min(diff(v))
  v[(b+1):length(v)] = 0
  w = sum(v)
  
  cnt = p-1
  g = hdrf[2:length(hdrf)] - hdrf[1:(length(hdrf)-1)]
  
  while(cnt > 0 && abs(g[cnt]) < 0.001) {      
    h = hdrf[cnt]
    p = cnt
    cnt = cnt-1
  }
  
  param = rep(0,3)
  param[1] = h
  param[2] = p
  param[3] = w
  
  return(param)
}
###########################################################################
# Modelling faecal egg count data (one-sample case) using Stan
###########################################################################

fec_stan<-function(fec,rawCounts=FALSE,CF=50,
                    zeroInflation=TRUE,muPrior,kappaPrior,phiPrior,
                    nsamples=12000,nburnin=2000,thinning=1,nchain=1,ncore=1,adaptdelta=0.9,verbose=FALSE){
  # checks from FECR_PoGa.R -------------------------------------------------
  #   if (sys.parent() == 0) env <- asNamespace("eggCounts") else env <- parent.frame()
  #   assign(".verboselevel", verbose*.verboselevel, envir = env)
  
  # number of faecal samples
  n <- length(fec)
  # check correction factors
  if (any(CF < 0)|(ceiling(CF)!=floor(CF)))    stop("correction factor(s) should be a positive integer", call.=FALSE)
  if(length(CF)>1 && length(CF)!=n) stop("Lengths of the vectors for FEC and correction factors do not match\n")
  
  # raw counts or EpGs?
  dilution <- CF;
  if(rawCounts){
    dilution <- 1
  }
  
  # check iteration parameters
  checkpars(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptdelta, zeroInflation, verbose)
  
  # divide data by correction factor
  if(any(fec %% dilution !=0)) stop(paste(c("Correction factor preCF does not match the given pre-treatment faecal egg counts. \n
  A possible correction factor (the largest common divisor) is"),mGCD(fec)))
  fec <- fec/dilution 
  
  # check model and set default values
  priors <- fec_setPrior(muPrior=muPrior, kappaPrior=kappaPrior, phiPrior=phiPrior)
  
  # set stan code and model name
  if(zeroInflation){code<-zinb_stan(priors);model<-"Zero-inflated Bayesian model"}
  if(!zeroInflation){code<-nb_stan(priors);model<-"Bayesian model without zero-inflation"}
                 
  # create data for stan use
  if (length(CF)==1) CF<-rep(CF,n)
  epg_data <-list(J=n, ystarraw = fec, CF=CF)
  if (length(setdiff(priors,fec_setPrior()))==0){
    if(zeroInflation){stanModel<-stanmodels$zinb}
    if(!zeroInflation){stanModel<-stanmodels$nb}
  } else {
  stanModel <- stan_model(model_name=paste(model),model_code=code)}
  # whether or not to suppress progress information and errors
  if (verbose){
    samples <- sampling(stanModel,data=epg_data,iter=nsamples, warmup=nburnin, chains=nchain,thin=thinning,control = list(adapt_delta = adaptdelta),cores=ncore)
  } else {
  samples <- suppressMessages(suppressWarnings(sampling(stanModel, data=epg_data,iter=nsamples, warmup=nburnin, chains=nchain,thin=thinning,control = list(adapt_delta = adaptdelta),cores=ncore,refresh=-1)))
  }
  
  # calculate samples
  if(zeroInflation){
  meanEPG<-rowMeans(extract(samples,"mui")[[1]])*(1-extract(samples,"phi")$phi)
  } else {meanEPG<-rowMeans(extract(samples,"mui")[[1]])}
  
  cat("Model: ", model,"\n","Number of Samples: ",nsamples, "\n","Warm-up samples: ",nburnin,"\n","Thinning: ",thinning,"\n")
  printSummary(cbind(meanEPG))
  return(invisible(samples))
}
###########################################################################
# Modelling the reduction in faecal egg count data (two-sample case) using Stan
###########################################################################

# check if the arguments in fecr_stan is sensible -------------------------

checkpars <- function(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptdelta, zeroInflation, verbose)
{
  if ((nburnin < 0)|(ceiling(nburnin)!=floor(nburnin)))
    stop("'nburnin' should be a positive integer", call.=FALSE)
  if ((nsamples < 1)|(ceiling(nsamples)!=floor(nsamples)))
    stop("'nsamples' should be a positive integer", call.=FALSE)
  if ((nsamples < nburnin))
    stop("The total number of samples must be greater than the number of burn-in samples", call.=FALSE)
  if ((thinning < 1)|(ceiling(thinning)!=floor(thinning)))
    stop("'thinning' should be a positive integer", call.=FALSE)
  if ((nchain < 1)|(ceiling(nchain)!=floor(nchain)))
    stop("The number of chains must be a positive integer", call.=FALSE)
  if ((ncore < 1)|(ceiling(ncore)!=floor(ncore)))
    stop("The number of cores must be a positive integer", call.=FALSE)
  if (!is.logical(rawCounts))
    stop("The rawCounts argument must be a logical", call.=FALSE)
  if (!is.logical(zeroInflation))
    stop("The zeroInflation argument must be a logical", call.=FALSE)
  if (!is.logical(verbose))
    stop("The verbose argument must be a logical", call.=FALSE)
  if ((adaptdelta<0)||(adaptdelta>1))
    stop("adaptdelta must be between 0 and 1", call.=FALSE)
  if (nburnin < 100)      cat("NOTE: 'nburnin' seems low.\n")
  if (nsamples < 1000)     cat("NOTE: 'nsamples' seems low.\n")
  invisible()
}

fecr_stan<-function(preFEC,postFEC,rawCounts=FALSE,preCF=50,postCF=preCF,
                    paired=TRUE,zeroInflation=TRUE,muPrior,kappaPrior,deltaPrior,phiPrior,
                    nsamples=12000,nburnin=2000,thinning=1,nchain=1,ncore=1,adaptdelta=0.9,verbose=FALSE){
# checks from FECR_PoGa.R -------------------------------------------------
#   if (sys.parent() == 0) env <- asNamespace("eggCounts") else env <- parent.frame()
#   assign(".verboselevel", verbose*.verboselevel, envir = env)
  
  # number of faecal samples
  preN <- length(preFEC)
  postN <- length(postFEC)

  # check correction factors
  if ((preCF < 0)|(ceiling(preCF)!=floor(preCF)))    stop("correction factor(s) should be a positive integer", call.=FALSE)
  if ((postCF < 0)|(ceiling(postCF)!=floor(postCF)))    stop("correction factor(s) should be a positive integer", call.=FALSE)
  if(length(preCF)>1 && length(preCF)!=preN) stop("Lengths of the vectors preCF and preFEC do not match\n")
  if(length(postCF)>1 && length(postCF)!=postN) stop("Lengths of the vectors postCF and postFEC do not match\n")
  
  # raw counts or EpGs?
  preDilution <- preCF; postDilution <- postCF
  if(rawCounts){
    preDilution <- postDilution <- 1
  }
  
  # check function arguments
  checkpars(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptdelta, zeroInflation, verbose)
  if (!is.logical(paired))
    stop("The paired argument must be a logical", call.=FALSE)
  
  # divide data by correction factor
  if(any(preFEC %% preDilution !=0)) stop(paste(c("Correction factor preCF does not match the given pre-treatment faecal egg counts. \n
  A possible correction factor (the largest common divisor) is"),mGCD(preFEC)))
  preFEC <- preFEC/preDilution
  if(any(postFEC %% postDilution !=0)) stop(paste(c("Correction factor postCF does not match the given post-treatment faecal egg counts. \n
  A possible correction factor (the largest common divisor) is"),mGCD(postFEC)))
  postFEC <- postFEC/postDilution
  
  if ( mean( preFEC) < mean( postFEC))
    cat("NOTE: mean of pre-treatment is smaller of post-treatment. Results may be unreliable.\n")
  if ( median( preFEC) < median( postFEC))
    cat("NOTE: median of pre-treatment is smaller of post-treatment. Results may be unreliable.\n")
  
  
  # set default values
  priors <- fecr_setPrior(muPrior=muPrior, kappaPrior=kappaPrior, deltaPrior=deltaPrior,phiPrior=phiPrior)
  
  # set update functions
  if(paired & zeroInflation){code<-ZI_paired_stan(priors);model<-"Zero-inflated Bayesian model for paired design"}
  if(paired & !zeroInflation){code<-paired_stan(priors);model<-"Bayesian model without zero-inflation for paired design"}
  if(!paired & zeroInflation){code<-ZI_unpaired_stan(priors);model<-"Zero-inflated Bayesian model for unpaired design"}
  if(!paired & !zeroInflation){code<-unpaired_stan(priors);model<-"Bayesian model without zero-inflation for unpaired design"}
  
  if(paired && preN != postN){
    stop("post sample size different to pre sample size\n")
  }
  # check counts
  if (preN==postN){if(sum( (preFEC-postFEC)^2) <= .Machine$double.eps)
    cat("Note: the pre-treatment and post-treatment counts are identical\n")}
  
  if (length(preCF)==1) preCF<-rep(preCF,preN)
  if (length(postCF)==1) postCF<-rep(postCF,postN)
  
  # create data list for stan use
  epg_data <-if(paired){
    list(J=preN, ystarbraw = preFEC, ystararaw = postFEC, fpre = preCF, fpost = postCF)
  } else {
      list(Jb=preN, Ja=postN, ystarbraw = preFEC, ystararaw = postFEC, fpre = preCF, fpost = postCF)}
  
  # whether or not to suppress progress information and errors
  if (length(setdiff(priors,fecr_setPrior()))==0){
    if(paired & zeroInflation){stanModel<-stanmodels$zipaired}
    if(paired & !zeroInflation){stanModel<-stanmodels$paired}
    if(!paired & zeroInflation){stanModel<-stanmodels$ziunpaired}
    if(!paired & !zeroInflation){stanModel<-stanmodels$unpaired}
  } else {
  stanModel <- stan_model(model_name=paste(model),model_code=code)}
  if (verbose){
    samples <- sampling(stanModel,data=epg_data,iter=nsamples, warmup=nburnin, chains=nchain,
                        thin=thinning,control = list(adapt_delta = adaptdelta),cores=ncore)
   } else {
  samples <- suppressMessages(suppressWarnings(sampling(stanModel,data=epg_data,iter=nsamples, warmup=nburnin, chains=nchain,
                              thin=thinning,control = list(adapt_delta = adaptdelta),cores=ncore,refresh=-1)))}
  
  # generate samples according to different models
  if(paired & !zeroInflation){
           meanEPG.untreated<-rowMeans(extract(samples,"mub")[[1]])
           meanEPG.treated<-rowMeans(extract(samples,"mub")[[1]])*extract(samples,"delta")$delta
           fecr<-1-extract(samples,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
         }
  if(!paired & !zeroInflation){
           meanEPG.untreated<-rowMeans(extract(samples,"mub")[[1]])
           meanEPG.treated<-rowMeans(extract(samples,"mua")[[1]])*extract(samples,"delta")$delta
           fecr<-1-extract(samples,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
         }
  if(paired & zeroInflation){
           meanEPG.untreated<-rowMeans(extract(samples,"mub")[[1]])*(1-extract(samples,"phi")$phi)
           meanEPG.treated<-rowMeans(extract(samples,"mub")[[1]])*extract(samples,"delta")$delta*(1-extract(samples,"phi")$phi)
           fecr<-1-extract(samples,"delta")$delta
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
         }
   if(!paired & zeroInflation){
           meanEPG.untreated<-rowMeans(extract(samples,"mub")[[1]])*(1-extract(samples,"phi")$phi)
           meanEPG.treated<-rowMeans(extract(samples,"mua")[[1]])*extract(samples,"delta")$delta*(1-extract(samples,"phi")$phi)
           fecr<-1-extract(samples,"delta")$delta
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
   }
  cat("Model: ", model,"\n","Number of Samples: ",nsamples, "\n","Warm-up samples: ",nburnin,"\n","Thinning: ",thinning,"\n")
  printSummary(result)
  return(invisible(samples))
}
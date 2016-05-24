## THIS COLLECTION OF FUNCTIONS IMPLEMENTS EMD VIA THE "EMD" PACKAGE
##AND ALSO RUNS THE ENSEMBLE EMD METHOD
##AND COMPLETE ENSEMBLE EMD METHOD

CEEMD <- function(sig, tt, noise.amp, trials, verbose = TRUE, spectral.method = "arctan", diff.lag = 1, tol = 5, max.sift = 200, stop.rule = "type5", boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max.imf = 100, interm = NULL, noise.type = "gaussian", noise.array = NULL)
{
        #Performs the Complete Ensemble Empirical Mode Decomposition as described in Torres et al (2011) A Complete Empirical Mode Decomposition with Adaptive Noise 
        #It runs EMD on a given signal for N=TRIALS
        #Each IMF set is saved to disk in TRIALS.DIR, which is created if it does not exist already.
        #Finally the EEMD function averages IMFs from all the trials together to produce an ensemble average and saves it in TRIALS.DIR 
        #INPUTS
        #   SIG is the time series
        #   TT is the sample times 
        #   NOISE.AMP is the amplitude of the noise distribution (upper/lower cutoff if uniform, standard deviation if gaussian, ignored otherwise)
        #   TRIALS is the number of times to run EMD
        #   VERBOSE - if TRUE, print progress
        #   SPECTRAL.METHOD defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation ("arctan")  
        #   or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain"); chain is dangerous at high frequencies
        #   DIFF.LAG specifies if you want to do naive differentiation (DIFF.LAG = 1), central difference method (DIFF.LAG = 2) or higher difference methods (DIFF.LAG > 2)
        #   MAX.SIFT stop sifting after this many times
        #   STOP.RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
        #   of envelope mean must be less than the user-specified tolerance level in the sense
        #   that the local average of upper and lower envelope is zero. The stopping rules
        #   type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
        #   of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
        #   and Wu (2008), respectively."
        #   BOUNDARY - how the beginning and end of the signal are handled
        #   SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
        #   SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
        #   SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
        #   MAX.IMF - How many IMFs are allowed, IMFs above this number will not be recorded
        #   INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.
        #   NOISE.TYPE - zero mean gaussian with standard deviation NOISE.AMP (if "gaussian" or unspecified), or "uniform" (absolute maximum/minimum amplitude), or "custom" (user defined noise matrix, expects NOISE.MATRIX)
        #   NOISE.ARRAY - A TRIALS x LENGTH(TT) matrix of user-defined noise to use in place of gaussian or uniform noise.  NOISE.TYPE must be set to "custom" for this input to be used.
     
    if(!(noise.type  %in% c("uniform", "gaussian", "custom"))) {
        stop(paste("Did not recognise noise.type option", noise.type,  "Please choose either ''uniform'' or ''gaussian''"))
    }

    if(noise.type == "custom") {
        if(!is.null(noise.array)) {
            if((dim(noise.array)[1] != trials) | dim(noise.array)[2] != length(tt)) {
                stop("You requested a custom noise array but either the number of rows did not equal the number of CEEMD trials or the number of columns did not equal the signal length, or both.")
            }
        } else {
            stop("If noise.type = \"custom\", then you must set noise.array equal to an array with the same number of rows as CEEMD trials and the same number of columns as signal samples.")
        }
     }

     #Make noise matrix and extract IMF 1
    if(noise.type == "uniform") {
        noise <- t(array(noise.amp * runif(length(sig) * trials), dim = c(length(sig), trials)))
        noise <- noise - mean(noise)
    } else if (noise.type == "gaussian") {
        noise <- t(array(noise.amp * rnorm(length(sig) * trials), dim = c(length(sig), trials)))
    } else if (noise.type == "custom") {
       noise <- noise.array
    }
     
     if(verbose) {
         print("Extracting IMF 1 from each noise/signal realization...")
     }

     imfs <- rep(0, length(sig))
     for(k in 1:trials) {
         imfs <- imfs + Sig2IMF(sig + noise[k, ], tt, spectral.method = spectral.method, diff.lag = diff.lag,
             tol = tol, max.sift = max.sift, stop.rule = stop.rule, boundary = boundary,
             sm = sm, smlevels = smlevels, spar = spar, max.imf = 1, interm = interm)$imf[, 1]
         if(verbose) {
             print(paste0("Trial ", k, " complete."))
         }
     }
     imfs <- imfs/trials

     if(verbose) {
         print("IMF 1 extracted.")
     }

     #Now decompose the noise into IMFs
     noise.imfs <- NULL
     if(verbose) {
          print("Decomposing noise series...")
     }

     for(k in 1:trials) {
         noise.imfs[[k]] <- Sig2IMF(noise[k, ], tt, spectral.method = spectral.method, diff.lag = diff.lag,
             tol = tol, max.sift = max.sift, stop.rule = stop.rule, boundary = boundary,
             sm = sm, smlevels = smlevels, spar = spar, max.imf = max.imf, interm = interm)$imf
         if(verbose) {
             print(paste0("Noise trial ", k, " complete."))
         }
     }
 
     #Continue extracting IMFs until either the maximum limit of IMFs is reached or the signal no longer has oscillatory elements
     r <- sig - imfs
     n.i <- 1
     noise.imf <- rep(0, length(sig))
     escape <- FALSE
     while(n.i < max.imf & EMD::extrema(r)$nextreme > 2) {
         imf.avg <- rep(0, length(sig))
         for(k in 1:trials) {
             #Extract only the first IMF of the series
             if(dim(noise.imfs[[k]])[2] < n.i ) {
                 warning(paste0("Attempted to extract more IMFs from the signal than are present in the noise series for trial ", k, "."))
                 noise.imf <- rep(0, length(sig))
             } else {
                 noise.imf <- noise.imfs[[k]][, n.i]
             }
             if(EMD::extrema(r + noise.imf)$nextreme > 2) {
                 imf.avg <- imf.avg + Sig2IMF(r + noise.imf, tt, 
                     spectral.method = spectral.method, diff.lag = diff.lag,
                     tol = tol, max.sift = max.sift, stop.rule = stop.rule, boundary = boundary,
                     sm = sm, smlevels = smlevels, spar = spar, max.imf = 1, interm = interm)$imf[, 1]
             } else {
                 imf.avg <- imf.avg + r + noise.imf
             }
             if(verbose) {
                print(paste("IMF", n.i + 1, "TRIAL", k))
             }
        }
        
        imf.avg <- imf.avg/trials #Average the IMF trials together
        imfs <- cbind(imfs, imf.avg)
        r <- r - imf.avg
        n.i <- n.i + 1
   }

   imfs.arr <- array(imfs, dim = c(length(sig), n.i))
   ceemd.result <- NULL
   ceemd.result$original.signal <- sig
   ceemd.result$residue <- r
   ceemd.result$tt <- tt
   ceemd.result$max.sift <- max.sift
   ceemd.result$tol <- tol
   ceemd.result$stop.rule <- stop.rule
   ceemd.result$boundary <- boundary
   ceemd.result$sm <- sm
   ceemd.result$smlevels <- smlevels
   ceemd.result$spar <- spar
   ceemd.result$max.imf <- max.imf
   ceemd.result$interm <- interm

   ceemd.result$hinstfreq <- array(0, dim = c(length(ceemd.result$original.signal), n.i))
   ceemd.result$hamp <- ceemd.result$hinstfreq
   ceemd.result$imf <- imfs.arr
   ceemd.result$nimf <- n.i
   for(i in 1:n.i)
   {
       imf = imfs.arr[,i]
       aimf = HilbertTransform(imf)
       ceemd.result$hinstfreq[, i] = InstantaneousFrequency(aimf, tt, method = spectral.method, lag = diff.lag)
       ceemd.result$hamp[, i] = HilbertEnvelope(aimf)
   }
   invisible(ceemd.result)

   
} 
CombineTrials <- function(in.dirs, out.dir, copy = TRUE)
{
   #Moves trial files from different directories, numbers them sequentially, and puts them in the specified directory.
   #This is important because the function EEMDCompile expects them to be numbered consecutively from 1, and will crash if this is not the case.
   #INPUTS
   #    IN.DIRS is a vector of paths to directories containing trial files.  Trial files will be found recursively (so trial files in subdirectories will be discovered and moved).
   #    OUT.DIR is a directory to put the combined trial set into.  If OUT.DIR does not exist, it will be created
   #    COPY asks if you want to copy files into OUT.DIR (TRUE) or move them from IN.DIRS to OUT.DIR (false)

   trial.file.pattern = "TRIAL_\\d+\\.?(RData|RDATA)$" 
   e1=simpleError(paste("Directory", out.dir, "is not empty!"))
   if(!file.exists(out.dir))
   {
       dir.create(out.dir, recursive = TRUE)
       cat("Created directory:", out.dir, "\n")
   }

   if(length(list.files(out.dir))>0)
   {
       stop(e1)
   }

   c = 1
   for(d in in.dirs)
   {
       if(!file.exists(d)) 
       {
           warning(cat("Trials directory", d, "does not exist and will be skipped."))
       }
       else
       {
            trial.files = list.files(d, pattern = trial.file.pattern, recursive = TRUE, full.names = TRUE)
            if(length(trial.files) == 0)
            {
                warning(cat("No EMD trial files found in directory", d))
            }
            else
            {
               for (trial.file in trial.files)
               {
                   if(copy)
                   {
                       res = file.copy(trial.file, paste(out.dir, "/", "TRIAL_",sprintf("%05i",c),".RData", sep = ""))
                       if(!res)
                       {
                           warning(cat("Failed to copy", trial.file))
                       }
                   }
                   else
                   {
                       file.rename(trial.file, cat(out.dir, "/", "TRIAL_",sprintf("%05i",c),".RData", sep = ""))
                   }

                   c = c + 1
               }
           }
       }
   }
}

EEMD <-function(sig, tt, noise.amp, trials, nimf, trials.dir = NULL, verbose = TRUE, spectral.method = "arctan", diff.lag = 1, tol = 5, max.sift = 200, stop.rule = "type5", boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max.imf = 100, interm = NULL, noise.type = "gaussian", noise.array = NULL) 
{
	#Performs the Ensemble Empirical Mode Decomposition as described in Huang and Wu (2008) A Review on the Hilbert Huang Transform Method and its Applications to Geophysical Studies
	#It runs EMD on a given signal for N=TRIALS
 	#Each IMF set is saved to disk in TRIALS.DIR, which is created if it does not exist already.
	#Finally the EEMD function averages IMFs from all the trials together to produce an ensemble average and saves it in TRIALS.DIR 
	#INPUTS
        #   SIG is the time series
        #   TT is the sample times 
        #   NOISE.AMP is the amplitude of the noise distribution (upper/lower cutoff if uniform, standard deviation if gaussian, ignored otherwise)
        #   TRIALS is the number of times to run EMD
        #   NIMF is the number of IMFs to be saved, IMFs past this number will be discarded
        #   TRIALS.DIR is the directory in which to put the EMD runs, if NULL, make a directory called "trials".  If TRIALS.DIR does not already exist, make it.
        #   VERBOSE - if TRUE, print progress
        #   SPECTRAL.METHOD defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation ("arctan")  
        #   or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain"); chain is dangerous at high frequencies
        #   DIFF.LAG specifies if you want to do naive differentiation (DIFF.LAG = 1), central difference method (DIFF.LAG = 2) or higher difference methods (DIFF.LAG > 2)
        #   MAX.SIFT stop sifting after this many times
        #   STOP.RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
        #   of envelope mean must be less than the user-specified tolerance level in the sense
        #   that the local average of upper and lower envelope is zero. The stopping rules
        #   type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
        #   of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
        #   and Wu (2008), respectively."
        #   BOUNDARY - how the beginning and end of the signal are handled
        #   SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
        #   SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
        #   SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
        #   MAX.IMF - How many IMFs are allowed, IMFs above this number will not be recorded
        #   INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.
        #   NOISE.TYPE - zero mean gaussian with standard deviation NOISE.AMP (if "gaussian" or unspecified), or "uniform" (absolute maximum/minimum amplitude), or "custom" (user defined noise matrix, expects NOISE.MATRIX)
        #   NOISE.ARRAY - A TRIALS x LENGTH(TT) matrix of user-defined noise to use in place of gaussian or uniform noise.  NOISE.TYPE must be set to "custom" for this input to be used.

        
	#OUTPUTS are saved to TRIALS.DIR in variable EMD.RESULT
	
	if(is.null(trials.dir))
	{
		trials.dir="trials"
	}

	if(!file.exists(trials.dir))
	{
		dir.create(trials.dir, recursive=TRUE)
		print(paste("Created trial directory:",trials.dir))
	}

        if(!(noise.type  %in% c("uniform", "gaussian", "custom"))) {
            stop(paste("Did not recognise noise.type option", noise.type,  "Please choose either ''uniform'' or ''gaussian''"))
        }

        if(noise.type == "custom") {
            if(!is.null(noise.array)) {
                if((dim(noise.array)[1] != trials) | dim(noise.array)[2] != length(tt)) {
                    stop("You requested a custom noise array but either the number of rows did not equal the number of EEMD trials or the number of columns did not equal the signal length, or both.")
                }
            } else {
                stop("If noise.type = \"custom\", then you must set noise.array equal to an array with the same number of rows as EEMD trials and the same number of columns as signal samples.")
            }
        }
     
	averaged.imfs=array(0,nimf*length(sig),dim=c(length(sig),nimf))
	averaged.noise=array(0,length(sig),dim=c(length(sig),1))
	averaged.residue=array(0,length(sig),dim=c(length(sig),1))
	for (j in 1:trials)
	{
                if(noise.type == "uniform") {	
		    noise=runif(length(sig),min=noise.amp*-1, max=noise.amp)
                } else if (noise.type == "gaussian") {
                    noise = rnorm(length(sig), mean = 0, sd = noise.amp)
                } else if (noise.type == "custom") {
                    noise <- noise.array[j, ]
                }
		tmpsig=sig+noise
		emd.result=Sig2IMF(tmpsig,tt, spectral.method = spectral.method, diff.lag = diff.lag, 
                   tol = tol, max.sift = max.sift, stop.rule = stop.rule, boundary = boundary, 
                   sm = sm, smlevels = smlevels, spar = spar, max.imf = max.imf, interm = interm)
		emd.result$noise=noise
		emd.result$original.signal=tmpsig-noise
		save(emd.result, file=paste(trials.dir, "/", "TRIAL_",sprintf("%05i",j),".RData",sep=""))
		if (emd.result$nimf<nimf)
		{
			trial.nimf=emd.result$nimf
		}
		else
		{
			trial.nimf=nimf
		}
		averaged.imfs[,1:trial.nimf]=averaged.imfs[,1:trial.nimf]+emd.result$imf[,1:trial.nimf]
		averaged.noise=averaged.noise+noise
		averaged.residue=averaged.residue+emd.result$residue
                if(verbose)
                {
		    print(paste("TRIAL",as.character(j),"OF",as.character(trials),"COMPLETE"))
                }
	}
}

EEMDCompile<-function(trials.dir, trials, nimf)
{
	#Averages trials together to produce a set of ensemble IMFs.
	#Produces the Hilbert spectrogram of each trial and puts it into a structure with all the other trials.
	#This is used later to generate an ensemble Hilbert spectrogram of the entire EEMD run.
	#INPUTS
	#	TRIALS.DIR is the location where the trial files produced by EEMD are stored
	#	TRIALS is the number of trials to average together
	#	NIMF is the number of IMFs to build
	#OUTPUTS
	#	EEMD.RESULT containes the ensemble IMFs

        trial.file.pattern = "TRIAL_\\d+\\.?(RData|RDATA)$"
	emd.result=NULL
	if(length(list.files(trials.dir, pattern = trial.file.pattern))==0)
	{
		stop(paste("No EMD trial files found in directory", trials.dir))
	}
			
	counter=1
	sind=1
        for(file.name in list.files(trials.dir, pattern=trial.file.pattern))
        {
		load(paste(trials.dir,"/",file.name,sep=""))
		
		
                if(counter==1)
                {
			siglen=length(emd.result$original.signal)
		        averaged.imfs=array(0,nimf*siglen,dim=c(siglen,nimf))
        		averaged.noise=array(0,siglen,dim=c(siglen,1))
        		averaged.residue=array(0,siglen,dim=c(siglen,1))
			hinstfreq=array(0,dim=c(length(emd.result$original.signal),nimf,trials))
			hamp=array(0,dim=c(length(emd.result$original.signal),nimf,trials))
                }
                if(emd.result$nimf>=nimf)
                {
                        imf.ind=nimf
                }
                else
                {
                        imf.ind=emd.result$nimf
                }
		averaged.imfs[,1:imf.ind]=averaged.imfs[,1:imf.ind]+emd.result$imf[,1:imf.ind]
		averaged.noise=averaged.noise+emd.result$noise
		averaged.residue=averaged.residue+emd.result$residue

                hinstfreq[,1:imf.ind,counter]=emd.result$hinstfreq[,1:imf.ind]
                hamp[,1:imf.ind,counter]=emd.result$hamp[,1:imf.ind]
                sind=sind+imf.ind

                counter=counter+1
                if(counter>trials)
                {
                        break
                }
        }
	counter=counter-1
        real.imfs = which(apply(array(as.logical(averaged.imfs), dim = dim(averaged.imfs)), 2, any) == TRUE) #Find out which IMFs have data in them
        if(length(real.imfs < nimf))
        {
            warning("The number of requested IMFs is greater than the maximum number of IMFs produced during individual EEMD trials.  Only IMFs with data will be recorded.")
        }

	averaged.imfs=averaged.imfs[, real.imfs] / counter
        hinstfreq = hinstfreq[, real.imfs, ]
        hamp = hamp[, real.imfs, ]
	averaged.noise=averaged.noise/counter
	averaged.residue=averaged.residue/counter
	
        if(counter<trials)
        {
                warning("Number of trials requested was greater than the number of trials found in the trials directory")
        }

        
	EEMD.result=c()
        EEMD.result$nimf=length(real.imfs)
        EEMD.result$tt=emd.result$tt
        EEMD.result$original.signal=emd.result$original.signal
	EEMD.result$averaged.imfs=averaged.imfs
	EEMD.result$averaged.noise=averaged.noise
	EEMD.result$averaged.residue=averaged.residue
        EEMD.result$hinstfreq=hinstfreq
        EEMD.result$hamp=hamp
	EEMD.result$trials=trials 
        invisible(EEMD.result)
}

EEMDResift <- function(EEMD.result, resift.rule, spectral.method = "arctan", diff.lag = 1, tol = 5, max.sift = 200, stop.rule = "type5", boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max.imf = 100, interm = NULL)
{
	#Resifts averaged IMFs generated by EEMD to generate valid IMFs for Hilbert Transform
	#INPUTS
	#	EEMD.RESULT contains the averaged IMF set generated from EEMD.
        #       RESIFT.RULE determines how the resifting occurs
        #       If resift.rule is numeric, get the nth IMF (so if resift.rule is 2, get the 2nd resifted IMF)
        #       If resift.rule is "last", return the last IMF.
        #       If resift.rule is "max.var", return the IMF with the most variance
        #       If resift.rule is "all", get all the IMFs returned by rerunning EMD on the averaged IMFs made by EEMD.
	#	This will likely be quite large.
        #       SPECTRAL.METHOD defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation ("arctan")  
        #       or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain"); chain is dangerous at high frequencies
        #       DIFF.LAG specifies if you want to do naive differentiation (DIFF.LAG = 1), central difference method (DIFF.LAG = 2) or higher difference methods (DIFF.LAG > 2)
        #       MAX.SIFT stop sifting after this many times
        #       STOP.RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
        #       of envelope mean must be less than the user-specified tolerance level in the sense
        #       that the local average of upper and lower envelope is zero. The stopping rules
        #       type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
        #       of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
        #       and Wu (2008), respectively."
        #       BOUNDARY - how the beginning and end of the signal are handled
        #       SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
        #       SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
        #       SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
        #       MAX.IMF - How many IMFs are allowed, IMFs above this number will not be recorded
        #       INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.

	#OUTPUTS
	#	EEMD.RESULT$IMF a set of IMFs that are generated from the EEMD imfs.

	resift.result=EEMD.result
	resift.result$imf=c()
        resift.result$hinstfreq = c()
        resift.result$hamp = c()
	resift.result$averaged.imfs=NULL

	if(!is.numeric(resift.rule) & !resift.rule %in% c("last", "max.var", "all"))
        {
               e=simpleError(paste("Did not recognize resift.rule:", resift.rule))
               stop(e)
        }		
	for(k in seq_len(dim(EEMD.result$averaged.imfs)[2]))
	{
		if(sum(EEMD.result$averaged.imfs[,k]==0)!=length(EEMD.result$averaged.imfs[,k]))
		{
			emd.result=Sig2IMF(EEMD.result$averaged.imfs[,k], EEMD.result$tt, spectral.method = spectral.method, diff.lag = diff.lag,
                            tol = tol, max.sift = max.sift, stop.rule = stop.rule, boundary = boundary,
                            sm = sm, smlevels = smlevels, spar = spar, max.imf = max.imf, interm = interm)
			if(is.numeric(resift.rule))
			{
				if(emd.result$nimf>=resift.rule)
				{
					resift.result$imf=cbind(resift.result$imf, emd.result$imf[,resift.rule])
                                        resift.result$hamp=cbind(resift.result$hamp, emd.result$hamp[,resift.rule])
                                        resift.result$hinstfreq=cbind(resift.result$hinstfreq, emd.result$hinstfreq[,resift.rule])
				}
				else
				{
					resift.result$imf=cbind(resift.result$imf, NA)
                                        resift.result$hamp=cbind(resift.result$hamp, NA)
                                        resift.result$hinstfreq=cbind(resift.result$hinstfreq, NA)

				}
			}
			else
			{
				if(resift.rule=="last")
				{
					resift.result$imf=cbind(resift.result$imf, emd.result$imf[,emd.result$nimf])
                                        resift.result$hamp=cbind(resift.result$hamp, emd.result$hamp[,emd.result$nimf])
                                        resift.result$hinstfreq=cbind(resift.result$hinstfreq, emd.result$hinstfreq[,emd.result$nimf])
				}
				
				if(resift.rule=="max.var")
				{
					var.list=c()
					for(j in seq_len(emd.result$nimf))
					{
						var.list=cbind(var.list, var(emd.result$imf[,j]))
					}
					
					resift.result$imf=cbind(resift.result$imf, emd.result$imf[,var.list==max(var.list)])
                                        resift.result$hamp=cbind(resift.result$hamp, emd.result$hamp[,var.list==max(var.list)])
                                        resift.result$hinstfreq=cbind(resift.result$hinstfreq, emd.result$hinstfreq[,var.list==max(var.list)])
				}

				if(resift.rule=="all")
				{
					resift.result$imf=cbind(resift.result$imf, emd.result$imf)
                                        resift.result$hamp=cbind(resift.result$hamp, emd.result$hamp)
                                        resift.result$hinstfreq=cbind(resift.result$hinstfreq, emd.result$hinstfreq)
				}
			}	
		}
	}

        resift.result$EEMD.trials = resift.result$trials 
        resift.result$trials = NULL
	resift.result$resift.max.sift = max.sift
        resift.result$resift.tol = tol
        resift.result$resift.stop.rule = stop.rule
        resift.result$resift.boundary = boundary
        resift.result$resift.sm = sm
        resift.result$resift.smlevels = smlevels
        resift.result$resift.spar = spar
        resift.result$resift.max.imf = max.imf
        resift.result$resift.interm = interm
	resift.result$resift.rule=resift.rule
	resift.result$nimf=dim(resift.result$imf)[2]
	invisible(resift.result)
}

Sig2IMF <- function(sig, tt, spectral.method = "arctan", diff.lag = 1, stop.rule = "type5", tol = 5, boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max.sift = 200, max.imf = 100, interm = NULL)

{
    #Extract IMFs
    #This function is intended to take data, recover IMFs, and save them
    #It calls and runs code developed by Donghoh Kim and Hee-Seok Oh as part of the "EMD" package available on CRAN.
    #Refer to their documentation for details on the EMD algorithm as implemented here.
    #INPUTS
    #	SIG is the time series
    #   TT is the sample times 
    #   SPECTRAL.METHOD defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation ("arctan")  
    #   or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain"); chain is dangerous at high frequencies
    #   DIFF.LAG specifies if you want to do naive differentiation (DIFF.LAG = 1), central difference method (DIFF.LAG = 2) or higher difference methods (DIFF.LAG > 2)
    #   TOL  determines what value is used to stop the sifting, this will depend on the chosen stop rule
    #   MAX.SIFT stop sifting after this many times
    #   STOP.RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
    #   of envelope mean must be less than the user-specified tolerance level in the sense
    #   that the local average of upper and lower envelope is zero. The stopping rules
    #   type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
    #   of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
    #   and Wu (2008), respectively."
    #   BOUNDARY - how the beginning and end of the signal are handled
    #   SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
    #   SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
    #   SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
    #   MAX.IMF - How many IMFs are allowed, IMFs above this number will not be recorded
    #   INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.
    #OUTPUT is a list containing the original signal, IMFs, and information on EMD parameters.
    #Danny Bowman
    #UNC Chapel Hill

    emd.result=EMD::emd(sig, tt, max.sift=max.sift, stoprule=stop.rule, tol=tol, 
        boundary=boundary,sm=sm,spar=spar, 
        check=FALSE, plot.imf=FALSE,max.imf=max.imf)
    emd.result$original.signal=sig
    emd.result$tt=tt
    emd.result$max.sift = max.sift
    emd.result$tol = tol
    emd.result$stop.rule = stop.rule
    emd.result$boundary = boundary
    emd.result$sm = sm
    emd.result$smlevels = smlevels
    emd.result$spar = spar
    emd.result$max.imf = max.imf
    emd.result$interm = interm
  
    emd.result$hinstfreq = array(0, dim = c(length(emd.result$original.signal), emd.result$nimf))
    emd.result$hamp = emd.result$hinstfreq
    
    if(emd.result$nimf > 0) {    
        for(i in seq(emd.result$nimf))
        {
            imf = emd.result$imf[,i]
            aimf = HilbertTransform(imf)
            emd.result$hinstfreq[, i] = InstantaneousFrequency(aimf, tt, method = spectral.method, lag = diff.lag)
            emd.result$hamp[, i] = HilbertEnvelope(aimf)
        }
    } else {
        warning("Sig2IMF says:  No IMFs were extracted, possibly because the signal lacks extrema.")
    }
    invisible(emd.result)
}

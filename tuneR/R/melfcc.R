# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

melfcc <- function(samples, sr=samples@samp.rate, wintime=0.025, hoptime=0.010,
            numcep=12, lifterexp=0.6, htklifter=FALSE,
            sumpower=TRUE, preemph=0.97, dither=FALSE,
            minfreq=0, maxfreq=sr/2, nbands=40, bwidth=1.0, 
            dcttype = c("t2", "t1", "t3", "t4"),
            fbtype = c("mel", "htkmel", "fcmel", "bark"), 
            usecmp=FALSE, modelorder=NULL, spec_out=FALSE,
            frames_in_rows=TRUE){

    dcttype <- match.arg(dcttype)
    fbtype <- match.arg(fbtype)
    
    if(!(is(samples, "Wave") || is(samples, "WaveMC"))) 
        stop("'samples' needs to be of class 'Wave' or 'WaveMC'")
    validObject(samples)

    if(nchannel(samples) > 1) 
        stop("Processing for more than one channel not yet implemented...")

    if(!is.null(modelorder) && !(modelorder==as.integer(modelorder) && modelorder > 0))
        stop("'modelorder' has to be a non-negative integer or NULL")

    if(!is.null(modelorder) && modelorder > 0 && numcep > modelorder+1)
        stop("No. of cepstra can't be larger than 'modelorder+1'")

    samples1 <- if(is(samples, "Wave")) samples@left else samples@.Data[,1]
    if(preemph != 0){
        ssamples <- as.vector(filter(samples1, filter=c(1, -preemph), 
            method="convolution", sides=1, circular=FALSE))
        ssamples[1] <- samples1[1]
    } else {
        ssamples <- samples1
    }

    # Compute FFT power spectrum
    pspectrum <- powspec(x=ssamples, sr=sr, wintime=wintime, steptime=hoptime, dither=dither)

    # Conversion to Mel/Bark scale
    aspectrum <- audspec(pspectrum=pspectrum, sr=sr, nfilts=nbands, 
                         fbtype=fbtype, minfreq=minfreq, maxfreq=maxfreq,
                         sumpower=sumpower, bwidth=bwidth)$aspectrum
    # PLP-like weighting and compression
    if(usecmp){
        aspectrum <- postaud(x=aspectrum, fmax=maxfreq, fbtype=fbtype)$y
    }

    #lpcas <- numeric(0)
    lpcas <- NULL
    if(!is.null(modelorder) && modelorder > 0){
        if(dcttype != "t1"){
            warning("PLP cepstra are implicitly dcttype 1")
        }

        # LPC/PLP
        lpcas <- dolpc(x=aspectrum, modelorder=modelorder)

        # Cepstra out of LPC/PLP
        cepstra <- lpc2cep(a=lpcas, nout=numcep)
    } else {
        # Cepstra via DCT
        cepstra <- spec2cep(spec=aspectrum, ncep=numcep, type=dcttype)$cep
    }

    # Liftering
    cepstra <- lifter(x=cepstra, lift=lifterexp, htk=htklifter)

    if(spec_out){
      if(frames_in_rows){
          res <- list(cepstra=t(cepstra), aspectrum=t(aspectrum),
                  pspectrum=t(pspectrum), lpcas = if(is.null(lpcas)) lpcas else t(lpcas))
      } else {
        res <- list(cepstra=cepstra, aspectrum=aspectrum, pspectrum=pspectrum,
                lpcas=lpcas)
      }
    } else {
      if(frames_in_rows){
        res <- t(cepstra)
      } else {
        res <- cepstra
      }
    }
    return(res)
}

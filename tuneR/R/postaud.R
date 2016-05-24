# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

postaud <- function(x, fmax, fbtype=c("bark", "mel", "htkmel", "fcmel"),
broaden=FALSE){

    if(!(is.numeric(x) && is.matrix(x)))
      stop("'x' has to be a numeric matrix")

    nbands <- nrow(x)
    nframes <- ncol(x)

    nfpts <- nbands + 2*broaden
    
    fbtype <- match.arg(fbtype)
    bandcfhz <- switch(fbtype,
                    bark = bark2hz(z=seq(0, hz2bark(f=fmax), length.out=nfpts)),
                    mel = mel2hz(z=seq(0, hz2mel(f=fmax), length.out=nfpts)),
                    htkmel = mel2hz(z=seq(0, hz2mel(f=fmax, htk=TRUE), length.out=nfpts), htk=TRUE),
                    fcmel = mel2hz(z=seq(0, hz2mel(f=fmax, htk=TRUE), length.out=nfpts), htk=TRUE),
    )
    bandcfhz <- bandcfhz[(1+broaden):(nfpts-broaden)]

    # Calculate equal loudness curve
    fsq <- bandcfhz^2
    ftmp <- fsq + 1.6e5
    eql <- ( (fsq/ftmp)^2) * ( (fsq + 1.44e6) / (fsq + 9.61e6) )

    # Apply equal loudness curve
    z <- matrix(rep(eql, nframes), ncol=nframes) * x
    
    # Compression
    z <- z^(.33)

    # Replicate first and last band (because they are unreliable as calculated)
    if(broaden){
        y <- z[c(1, 1:nbands, nbands),]
    } else {
        y <- z[c(2, 2:(nbands-1), nbands-1),]
    }
    return(list(y=y, eql=eql))
}


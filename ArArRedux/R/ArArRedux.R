
#' Create a new \code{\link{redux}} object
#'
#' Initialises a new \code{\link{redux}} object by packing a
#' \code{\link{logratios}} dataset together with all the parameters
#' needed for age calculation
#' 
#' @param X an object of class \code{\link{logratios}}
#' @param Jpos a vector of integers denoting the positions of the
#' fluence monitors in the irradiation stack
#' @param detectors a list of strings denoting the detectors for each
#' argon isotope
#' @return an object of class \code{\link{redux}}
#' @export
newredux <- function(X,Jpos,detectors=
            list(Ar36="H1",Ar37="L2",Ar38="L1",Ar39="AX",Ar40="H1")){
    out <- X
    out$Jpos <- Jpos
    out$detectors <- detectors
    out$param$l0 = 5.5492e-4 # 40K decay constant [Ma-1]
    out$param$sl0 = 0.0047e-4 # Renne et al. (2010)
    out$param$l7 = 365.25*log(2)/34.95 # 37Ar decay constant [yr-1]
    out$param$sl7 = out$param$l7*0.04/34.95 # Renne and Norman (2001)
    out$param$l9 = log(2)/269 # 39Ar decay constant [a-1]
    out$param$sl9 = out$param$l9*1.5/269 # Stoenner et al. (1965)
    out$param$l6 = log(2)/301200 # 36Cl decay constant [a-1]
    out$param$sl6 = out$param$l6*1000/301200 # uncertainty
    out$param$pcl = 252.7 # P(36Cl/38Cl) for OSTR reactor
    out$param$spcl = 1.8 # Renne et al. (2008)
    out$param$ts = 28.201 # age of the Fish Canyon Tuff
    out$param$sts = 0.023 # 1-sigma age uncertainty
    out$param$air = 298.56 # atmospheric 40/36-ratio
    out$param$sair = 0.155 # Lee et al. (2006)
    class(out) <- "redux"
    return(out)
}

#' Set or get Ar-Ar_Redux parameters
#'
#' This function is used to query and modify the half lives, standard
#' ages etc. associated with an object of class \code{\link{redux}}
#'
#' \code{\link{param}} grants access to the following parameters:
#'  
#' \code{l0}: 40K decay constant (default value = 5.5492e-4 Ma-1,
#' Renne et al. [2010])\cr
#' \code{sl0}: standard error of the 40K decay constant (default value
#' = 0.0047e-4 Ma-1)\cr
#' \code{l7}: 37Ar decay constant (default value = 7.2438 yr-1, Renne
#' and Norman [2001])\cr
#' \code{sl7}: standard error of the 37Ar decay constant (default
#' value = 0.0083 yr-1)\cr
#' \code{l9}: 39Ar decay constant (0.002577 yr-1 Stoenner et
#' al. [1965])\cr
#' \code{sl9}: standard error of the 39Ar decay constant (0.000014
#' yr-1)\cr
#' \code{l6}: 36Cl decay constant (default value = 2301.3e-9 yr-1)\cr
#' \code{sl6}: standard error of the 36Cl decay constant (default
#' value = 7.6e-9 yr-1\cr
#' \code{pcl}: (36Cl/38Cl)-production rate (default value = 252.7 for
#' OSTR reactor, Renne et al. [2008])\cr
#' \code{spcl}: standard error of the (36Cl/38Cl)-production rate
#' (default value = 1.8)\cr
#' \code{ts}: age of the fluence monitor (default = 28.201 Myr for the
#' Fish Canyon Tuff, Kuiper et al. [2008])\cr
#' \code{sts}: standard error of the fluence monitor age (default
#' value = 0.023 Myr)\cr
#' \code{air}: atmospheric 40Ar/36Ar ratio (default value = 298.56,
#' Lee et al. [2006])\cr
#' \code{sair}: standard error of the atmospheric 40Ar/36Ar ratio
#' (default value = 0.155)
#' 
#' @param X an object of class \code{\link{redux}}
#' @param ... any combination of the parameters given below
#' @return returns the modified \code{\link{redux}} object OR the
#' current parameter values if no optional arguments are supplied.
#' @examples
#' data(Melbourne)
#' param(Melbourne$X)$air
#' Y <- param(Melbourne$X,air=295.5)
#' param(Y)$air
#' @export
param <- function(X,...){
    arguments <- list(...)
    if (length(arguments)>0){
        X$param[names(arguments)] <- arguments
        return(X)
    } else {
        return(X$param)
    }
}

#' Plot a time resolved mass spectrometry signal
#'
#' Plots the raw signal of a given isotope against time.
#' 
#' @param x an object of class \code{\link{timeresolved}} or
#' \code{\link{PHdata}}
#' @param label a string with the name of the run
#' @param mass a string indicating the isotope of interest
#' @param ... optional parameters
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' mMC <- loaddata(samplefile,masses)
#' plot(mMC,"MD2-1a","Ar40")
#' @rdname plot
#' @export
plot.timeresolved <- function(x,label,mass,...){
    j <- which(x$labels==label)-1
    if (length(j)!=1){
        print('invalid input into plot function')
        return(NA)
    }
    k <- which(x$masses==mass)
    i <- j*nmasses(x)+k
    graphics::plot(x$thetime[,i],x$d[,i],type='p')
}
#' @examples
#' mPH <- loaddata(samplefile,masses,PH=TRUE)
#' plot(mPH,"MD2-1a","Ar40")
#' @rdname plot
#' @export
plot.PHdata <- function(x,label,mass,...){
    plot.timeresolved(x$signals[[mass]],label,mass,...)
}

# returns the indices of timevec2 which are closest to timevec1
nearest <- function(timevec1,timevec2){
    ii <- NULL
    for (i in 1:length(timevec1)) {
        ii <- c(ii, which.min(abs(timevec2 - timevec1[i])))
    }
    return(ii)
}

cast <- function(x,...){ UseMethod("cast",x) }
cast.default <- function(x,...){stop()}
cast <- function(from,to){
    out <- from
    class(out) <- to
    if (to == "logratios"){
        nruns <- nruns(out)
        num <- from$masses
        den <- rep(from$denmass,nmasses(from))
        out$num <- rep(num,nruns)
        out$den <- rep(den,nruns)
        out$masses <- NULL
        out$denmass <- NULL
        out$blankindices <- NULL
        out$nlr <- rep(nmasses(from),nruns)
    }
    return(out)
}

readthedate <- function(x){
    thedate <- strptime(x,"%d-%b-%Y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d/%m/%y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d-%b-%Y")
    return(as.numeric(thedate))
}

# masses = vector of strings
# cirr = column with irradiation can label
# cpos = column with position in irradiation stack
# clabel = column containing the sample labels
# cdate = column containing the date
# ci = matrix with column indices of the masses of interest
newtimeresolved <- function(thetable,masses,cirr,cpos,clabel,cdate,ci){
    out <- list()
    class(out) <- "timeresolved"
    out$masses <- masses
    out$irr <- as.character(thetable[,cirr])
    out$pos <- as.numeric(thetable[,cpos])
    out$labels <- as.character(thetable[,clabel])
    out$thedate <- readthedate(thetable[,cdate])
    nmasses <- nmasses(out)
    nruns <- nruns(out)
    ncycles <- dim(ci)[1]
    out$d <- matrix(0,nrow=ncycles,ncol=nmasses*nruns)
    out$thetime <- t(thetable[,ci[,1]+1])
    for (i in 1:ncycles){
        for (j in 1:nmasses){
            k <- seq(from=j,to=nmasses*nruns,by=nmasses)
            out$d[i,k] <- thetable[,ci[i,j]]
        }
    }
    return(out)
}

newPHdata <- function(thetable,masses,cirr,cpos,clabel,cdate,ci){
    out <- list()
    class(out) <- "PHdata"
    out$masses <- masses
    for (i in 1:nmasses(out)){
        out$signals[[masses[i]]] <-
            newtimeresolved(thetable,masses[i],
                cirr,cpos,clabel,cdate,as.matrix(ci[,i]))
    }
    return(out)
}

nmasses <- function(x){
    return(length(x$masses))
}

nruns <- function(x,...){ UseMethod("nruns",x) }
nruns.default <- function(x,...){
    length(x$labels)
}
nruns.PHdata <- function(x,...){
    return(nruns(x$signals[[1]]))
}

ncycles <- function(x,...){ UseMethod("ncycles",x) }
ncycles.default <- function(x,...){stop()}
ncycles.timeresolved <- function(x,...){
    return(dim(x$thetime)[1])
}

getsignal <- function(X,prefix,num=NULL){
    i <- getindices(X,prefix,num)
    return(cbind(X$intercepts[i], sqrt(X$covmat[i,i])))
}

getindices <- function(...){ UseMethod("getindices") }
getindices.default <- function(nmasses,nruns=NULL,
                               imasses=NULL,iruns=NULL,...){
    if (is.null(nruns)) nruns <- max(iruns)
    if (is.null(imasses)) imasses <- 1:nmasses
    if (is.null(iruns)) iruns <- 1:nruns
    i <- NULL
    for (irun in iruns){
        i <- c(i,(irun-1)*nmasses + imasses)
    }
    return(i)
}
getindices.logratios <- function(x,iruns,...){
    cs <- c(0,cumsum(x$nlr))
    i <- NULL
    for (irun in iruns){
        i <- c(i,(cs[irun]+1):cs[irun+1])
    }
    return(i)
}
getindices.redux <- function(X,prefix=NULL,num=NULL,den=NULL,
                             pos=NULL,exact=TRUE,invert=FALSE,...){
    i <- 1:length(X$intercepts)
    if (is.null(prefix)) {
        i1 <- i
    } else {
        if (exact){
            matches <- X$labels %in% prefix
            if (invert){
                j <- which(!matches)
            } else {
                j <- which(matches)
            }
        } else {
            j <- array(grep(prefix, X$labels, invert))
        }
        i1 <- getindices.logratios(X,j)
    }
    if (is.null(num)) {
        i2 <- i
    } else {
        if (exact){
            matches <- X$num %in% num
            if (invert){
                i2 <- which(!matches)
            } else {
                i2 <- which(matches)
            }
        } else {
            i2 <- array(grep(num, X$num, invert))
        }
    }
    if (is.null(den)) {
        i3 <- i
    } else {
        if (exact){
            matches <- X$den %in% den
            if (invert){
                i3 <- which(!matches)
            } else {
                i3 <- which(matches)
            }
        } else {
            i3 <- array(grep(den, X$den, invert))
        }
    }
    if (is.null(pos)) {
        i4 <- i
    } else {
        matches <- X$pos %in% pos
        if (invert){
            i4 <- which(!matches)
        } else {
            i4 <- which(matches)
        }
    }
    return(which((i %in% i1) & (i %in% i2) &
                 (i %in% i3) & (i %in% i4)))
}

getrunindices <- function(X,prefixes,invert=FALSE){
    ns <- nruns(X)
    i <- c() # vector of intercepts
    for (prefix in prefixes){
        i <- c(i,grep(prefix, X$labels))
    }
    j <- sort(i)
    if (invert){
        return((1:ns)[-j])
    } else {
        return(j)
    }
}

getruns <- function(x,...){ UseMethod("getruns",x) }
getruns.default <- function(x,...){stop()}
getruns.timeresolved <- function(x,i,...){
    ii <- getindices(nmasses=nmasses(x),iruns=i)
    return(x$d[,ii])
}

#' Select a subset of some data
#'
#' Extracts those intercepts, covariances etc. that match a given list
#' of indices or labels.
#' 
#' @param x an object of class \code{\link{timeresolved}},
#' \code{\link{logratios}}, \code{\link{redux}} or
#' \code{\link{results}}
#' @param i a vector with indices of the selected runs
#' @param labels a string or a vector of strings with sample names
#' @param ... other arguments
#' @return an object of the same class as \code{x}
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' MD <- subset(ages,labels=c("MD2-1","MD2-2","MD2-3","MD2-4","MD2-5"))
#' plotcorr(MD)
#' @rdname subset
#' @export
subset.timeresolved <- function(x,i=NULL,labels=NULL,...){
    if (is.null(i)) i <- which(x$labels %in% labels)
    out <- x
    out$d <- getruns(x,i)
    out$thetime <- x$thetime[,i]
    out$thedate <- x$thedate[i]
    out$irr <- x$irr[i]
    out$pos <- x$pos[i]
    out$labels <- x$labels[i]
    if (methods::is(x,"blankcorrected")){ out$blankindices <- x$blankindices[i] }
    return(out)
}
#' @rdname subset
#' @export
subset.logratios <- function(x,i=NULL,labels=NULL,...){
    if (is.null(i)) i <- which(x$labels %in% labels)
    out <- x
    out$irr <- x$irr[i]
    out$pos <- x$pos[i]
    out$labels <- x$labels[i]
    out$thedate <- x$thedate[i]
    out$nlr <- x$nlr[i]
    j <- getindices.logratios(x,i)
    out$num <- x$num[j]
    out$den <- x$den[j]
    out$intercepts <- x$intercepts[j]
    out$covmat <- x$covmat[j,j]
    return(out)
}
#' @rdname subset
#' @export
subset.redux <- function(x,i=NULL,labels=NULL,...){
    if (is.null(i)) i <- which(x$labels %in% labels)
    return(subset.logratios(x,i))
}
#' @rdname subset
#' @export
subset.results <- function(x,i=NULL,labels=NULL,...){
    out <- x
    if (is.null(i)) i <- which(x$labels %in% labels)
    out$labels <- x$labels[i]
    out$thedate <- x$thedate[i]
    out$ages <- x$ages[i]
    out$covmat <- x$covmat[i,]
    out$covmat <- out$covmat[,i]
    return(out)
}

#' Apply a blank correction
#'
#' Applies a blank correction to some time-resolved mass spectrometer data
#' 
#' @param x an object of class \code{\link{timeresolved}} or
#' \code{\link{PHdata}}
#' @param ... other arguments
#' @param blanklabel as string denoting the prefix of the blanks
#' @param prefix a string to be prepended to the non-blanks
#' @return an object of class \code{\link{blankcorrected}}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' blanklabel <- "EXB#"
#' l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
#' plotcorr(l)
#' @export
blankcorr <- function(x,...){ UseMethod("blankcorr",x) }
#' @rdname blankcorr
#' @export
blankcorr.default <- function(x,...){stop()}
#' @rdname blankcorr
#' @export
blankcorr.timeresolved <- function(x,blanklabel=NULL,prefix='',...){
    if (is.null(blanklabel)){
        out <- x
        out$blankindices <- 1:nruns(x)
    } else {
        # find indices of the blanks and non-blanks
        iblanks <- array(grep(blanklabel,x$labels))
        iothers <- array(grep(blanklabel,x$labels,invert=TRUE))
        blanks <- subset(x,iblanks)
        others <- subset(x,iothers)
        out <- others
        out$labels <- unlist(lapply(prefix,paste0,others$labels))
        inearestblanks <- nearest(others$thedate,blanks$thedate)
        out$d <- others$d - getruns(blanks,inearestblanks)
        out$blankindices <- as.vector(inearestblanks)
    }
    class(out) <- append(class(out),"blankcorrected")
    return(out)
}
#' @rdname blankcorr
#' @export
blankcorr.PHdata <- function(x,blanklabel=NULL,prefix='',...){
    out <- x
    for (mass in out$masses){
        out$signals[[mass]] <-
            blankcorr.timeresolved(out$signals[[mass]],blanklabel,prefix)
    }
    return(out)
}

#' Select a subset of isotopes from a dataset
#'
#' Extracts the intercepts, covariance matrix, etc. of a selection of
#' isotopes from a larger dataset
#' 
#' @param x an object of class \code{\link{logratios}},
#' \code{\link{timeresolved}}, \code{\link{PHdata}} or
#' \code{\link{redux}}.
#' @param ... other arguments
#' @return an object of the same class as x
#' @examples
#' kfile <- system.file("K-glass.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' mk <- loaddata(kfile,masses)
#' lk <- fitlogratios(blankcorr(mk,"EXB#","K:"),"Ar40")
#' k <- getmasses(lk,"Ar39","Ar40") # subset on the relevant isotopes
#' plotcorr(k)
#' @export
getmasses <- function(x,...){ UseMethod("getmasses",x) }
#' @rdname getmasses
#' @export
getmasses.default <- function(x,...){ stop() }
#' @param mass a vector of strings denoting the masses of interest
#' @param invert boolean parameter indicating whether the selection
#' should be inverted (default = FALSE)
#' @rdname getmasses
#' @export
getmasses.timeresolved <- function(x,mass,invert=FALSE,...){
    out <- x
    if (invert){ imasses <- which(x$masses != mass) }
    else       { imasses <- which(x$masses == mass) }
    out$masses <- out$masses[imasses]
    ii <- getindices(nmasses(x),nruns(x),imasses)
    out$d <- x$d[,ii]
    return(out)
}
#' @param num vector of strings indicating the numerator isotopes
#' @param den vector of string indicating the denominator isotopes
#' @rdname getmasses
#' @export
getmasses.logratios <- function(x,num,den,invert=FALSE,...){
    out <- x
    if (invert){
        i <- which(!((x$num %in% num) & (x$den %in% den)))
    } else {
        i <- which((x$num %in% num) & (x$den %in% den))
    }
    out$num <- x$num[i]
    out$den <- x$den[i]
    out$intercepts <- x$intercepts[i]
    out$covmat <- x$covmat[i,i]
    out$nlr <- graphics::hist(i,breaks=c(0,cumsum(x$nlr)),
                    plot=FALSE)$counts
    return(out)
}

setmasses <- function(x,mass,value){ UseMethod("setmasses",x) }
setmasses.default <- function(x,mass,value){stop()}
setmasses.timeresolved <- function(x,mass,value){
    imasses <- which(x$masses == mass)
    ii <- getindices(nmasses(x),nruns(x),imasses)
    x$d[,ii] <- value
    return(x)
}
setmasses.fit <- function(x,mass,value){
    imasses <- which(x$masses == mass)
    ii <- getindices(nmasses(x),nruns(x),imasses)
    x$intercepts[ii] <- value$intercepts
    for (i in 1:length(ii)){
        x$covmat[ii[i],ii] <- value$covmat[i,]
    }
    return(x)
}

replacenegatives <- function(x){
    out <- x
    nmasses <- nmasses(x)
    nruns <- nruns(x)
    isnegative <- apply(x$d<0,2,"sum")>0
    ntoreplace <- sum(isnegative)
    out$d[,isnegative] <- # effectively set to zero
    seq(from=1e-18,to=1e-20,length.out=ncycles(x)*ntoreplace)
    return(out)
}

Jtakeratios <- function(nruns,inum,iden){
    nlr <- length(inum)
    nmasses <- nlr+1
    J <- matrix(0,nrow=nlr,ncol=nmasses) # elementary matrix
    for (i in 1:length(inum)){
        J[i,inum[i]] <- 1
        J[,iden] <- -1
    }
    out <- matrix(0,nrow=nlr*nruns,ncol=nmasses*nruns)
    for (i in 1:nruns){
        irow <- ((i-1)*nlr+1):(i*nlr)
        icol <- ((i-1)*nmasses+1):(i*nmasses)
        out[irow,icol] <- J
    }
    return(out)
}

# x = object of class "logratios"
# ireplace = vector of indices with runs that need replacing
Jmean <- function(nlr,ireplace){
    nruns <- length(nlr)
    idontreplace <- which(!((1:nruns) %in% ireplace))
    nlrmean <- nlr[ireplace[1]]
    nlrout <- c(nlrmean,nlr[-ireplace])
    n <- length(ireplace)
    J <- matrix(0,nrow=sum(nlrout),ncol=sum(nlr))
    for (i in ireplace){ # first calculate the mean
        j <- 1:nlrmean
        k <- getindices.logratios(list(nlr=nlr),i)
        J[j,k] <- diag(nlrmean)/n
    }
    if (length(idontreplace)==0) return(J)
    for (i in 1:length(idontreplace)){
        j <- getindices.logratios(list(nlr=nlrout),i+1)
        k <- getindices.logratios(list(nlr=nlr),idontreplace[i])
        J[j,k] <- diag(nlrout[i+1])
    }
    return(J)
}

takeratios <- function(x,...){ UseMethod("takeratios",x) }
takeratios.default <- function(x,...){stop()}
takeratios.timeresolved <- function(x,denmass,...){
    den <- getmasses(x,denmass)
    num <- getmasses(x,denmass,invert=TRUE)
    out <- num
    for (mass in num$masses){ # loop through the isotopes
        ratio <- getmasses(num,mass)$d/den$d
        out <- setmasses(out,mass,ratio)
    }
    out$denmass <- denmass
    class(out) <- append(class(out),"ratio")
    return(out)
}
takeratios.fit <- function(x,denmass,...){
    out <- x
    iden <- which(x$mass == denmass)
    inum <- which(x$mass != denmass)
    J <- Jtakeratios(nruns(x),inum,iden)
    out$masses <- x$masses[inum]
    out$intercepts <- J %*% x$intercepts
    out$covmat <- J %*% x$covmat %*% t(J)
    out$denmass <- denmass
    return(out)
}

takelogs <- function(x){
    out <- x
    out <- replacenegatives(x)
    out$d <- log(out$d)
    class(out) <- append(class(out),"logged")
    return(out)
}

#' Extrapolation to 'time zero'
#'
#' This function extrapolates time resolved mass spectrometer data to
#' t=0. When fed with multicollector data, it forms the ratios of the
#' raw signals, forms their logs and performs linear regression to t=0
#' When fed with single collector data, the function first takes their
#' logs and extrapolates them to t=0 before taking ratios, unless
#' \code{denmass}=NULL, in which case the logs of the raw signals are
#' extrapolated.
#' 
#' @param x an object of class \code{timeresolved} or \code{PHdata}
#' @param ... further arguments (see below)
#' @return an object of class \code{logratios}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' blanklabel <- "EXB#"
#' l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
#' plotcorr(l)
#' @export
fitlogratios <- function(x,...){ UseMethod("fitlogratios",x) }
#' @rdname fitlogratios
#' @export
fitlogratios.default <- function(x,...){ stop() }
#' @param denmass a string denoting the denominator isotope
#' @rdname fitlogratios
#' @export
fitlogratios.timeresolved <- function(x,denmass,...){
    r <- takeratios(x,denmass)
    l <- takelogs(r)
    f <- timezero(l)
    return(cast(f,"logratios"))
}
#' @rdname fitlogratios
#' @export
fitlogratios.PHdata <- function(x,denmass=NULL,...){
    z <- newfit(x)
    for (i in 1:nmasses(x)){
        l <- takelogs(x$signals[[i]])
        z <- setmasses(z,z$masses[i],timezero(l))
    }
    if (is.null(denmass)) {
        f <- z
        f$denmass <- NA
    } else {
        f <- takeratios(z,denmass)
    }
    return(cast(f,"logratios"))
}

newfit <- function(x,...){ UseMethod("newfit",x) }
newfit.default <- function(x,...){ stop() }
newfit.timeresolved <- function(x,nmasses=NULL,nruns=NULL,...){
    out <- x
    class(out) <- "fit"
    if (is.null(nmasses)) nmasses <- nmasses(x)
    if (is.null(nruns)) nruns <- nruns(x)
    out$d <- NULL
    out$thetime <- NULL
    out$intercepts <- rep(0,nmasses*nruns)
    out$covmat <- matrix(0,nrow=nmasses*nruns,ncol=nmasses*nruns)
    return(out)
}
newfit.PHdata <- function(x,...){
    x1 <- x$signals[[1]]
    out <- newfit(x1,nmasses(x),nruns(x1))
    out$masses <- x$masses
    out$signals <- NULL
    return(out)
}

timezero <- function(x){
    nmasses <- nmasses(x)
    out <- newfit(x,nmasses)
    bi <- rle(as.vector(x$blankindices))$values # blank indices
    irunsout <- 0
    for (i in bi){ # loop through the groups
        irunsx <- which(i==x$blankindices)
        irunsout <- utils::tail(irunsout,n=1) + (1:length(irunsx))
        g <- subset(x,irunsx) # extract group
        out <- setfit(out,fit(g),nmasses,irunsout)
    }
    return(out)
}

# x and f are both objects of class "fit"
# nmasses = number of masses or isotope ratios per sample
# iruns = indicates which samples need replacing
setfit <- function(x,f,nmasses,iruns){
    out <- x
    ii <- getindices(nmasses=nmasses,iruns=iruns)
    out$intercepts[ii] <- f$intercepts
    out$covmat[ii,ii] <- f$covmat
    return(out)
}

fit <- function(x){
    out <- newfit(x)
    if (nruns(x)>1) { # average time for all samples in group
        thetime <- apply(x$thetime,1,"mean")
    } else {
        thetime <- x$thetime
    }
    f <- stats::lm(x$d ~ thetime)
    if (nruns(x) == 1 & nmasses(x) == 1) { # only one intercept
        out$intercepts <- stats::coef(f)["(Intercept)"]
        covcolname <- "(Intercept)"
    } else { # an entire row of intercepts
        out$intercepts <- stats::coef(f)["(Intercept)",]
        covcolname <- ":(Intercept)"
    }
    myvcov <- stats::vcov(f)
    j <- which(colnames(myvcov)==covcolname)
    out$covmat <- myvcov[j,j]
    return(out)
}

#' Load mass spectrometer data
#'
#' Loads a .csv file with raw mass spectrometer data
#' 
#' @param fname the file name, must end with .csv
#' @param masses a vector of strings denoting the order of the
#' isotopes listed in the table
#' @param MS the type of mass spectrometer
#' @param PH a boolean indicating whether the data are to be treated
#' as multicollector (PH=FALSE) or 'peak hopping' (PH=TRUE) data. The
#' default is PH=FALSE.
#' @return if PH=FALSE: an object of class \code{timeresolved}\cr
#' if PH=TRUE: an object of class \code{PHdata}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' plot(m,"MD2-1a","Ar40")
#' @export
loaddata <- function(fname,masses,MS="ARGUS-VI",PH=FALSE){
    thetable <- utils::read.csv(file=fname,header=FALSE,skip=3)
    nrows <- dim(thetable)[1] # number of MS runs
    ncols <- dim(thetable)[2] # number of columns
    nmass <- length(masses) # number of isotopes
    ncycles <- (ncols-3)/(2*nmass)
    cirr <- 1
    cpos <- 2
    clabel <- 3
    cdate <- 4 # column with the date
    ci <- matrix(NA,nrow=ncycles,ncol=nmass) # column indices
    for (i in 1:nmass){
        ci[,i] <- seq(from=3+2*i,to=ncols,by=2*nmass)
    }
    if (PH) {
        out <- newPHdata(thetable,masses,cirr,cpos,clabel,cdate,ci)
    } else {
        out <- newtimeresolved(thetable,masses,cirr,cpos,clabel,cdate,ci)
    }
    return(out)
}


theday <- function(thedate){
    dateinseconds <- as.POSIXlt(thedate,origin="1970-01-01 00:00:00")
    dateindays <- (1900+dateinseconds$year-1970)*365 + dateinseconds$yday
    dayinseconds <- dateindays*24*3600
    return(dayinseconds)
}

#' Calculate the arithmetic mean
#'
#' Calculate the arithmetic mean of some logratio data
#' 
#' @param x an object of class \code{redux} or \code{logratios}
#' @param i (optional) vector of sample indices
#' @param newlabel (optional) string with the new label to be assigned
#' to the averaged values
#' @return an object of the same class as \code{x}
#' @examples
#' data(Melbourne)
#' K <- average(Melbourne$X,grep("K:",Melbourne$X$labels),newlabel="K-glass")
#' plotcorr(K)
#' @export
average <- function(x,i=NULL,newlabel=NULL){
    if (length(i)==0) return(x)
    if (is.null(i)) i <- 1:nruns(x)
    if (is.null(newlabel)) newlabel <- as.character(x$labels[i[1]])
    j <- which(!((1:nruns(x)) %in% i))
    out <- subset(x,c(i[1],j))
    out$labels[1] <- newlabel
    J <- Jmean(x$nlr,i)
    out$intercepts <- J %*% x$intercepts
    out$covmat <- J %*% x$covmat %*% t(J)
    return(out)
}

#' Average all the data collected on the same day.
#'
#' This function is useful for grouping a number of replicate air
#' shots or calibration experiments
#' 
#' @param x an object of class \code{timeresolved}, \code{logratios},
#' \code{PHdata} or \code{redux}
#' @param newlabel a string with the new label that should be given to
#' the average
#' @return an object of the same class as x
#' @examples
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' dlabels <- c("H1","AX","L1","L2")
#' md <- loaddata(dfile,dlabels,PH=TRUE)
#' ld <- fitlogratios(blankcorr(md))
#' d <- averagebyday(ld,"DCAL")
#' plotcorr(d)
#' @export
averagebyday <- function(x,newlabel){
    out <- x
    thedays <- rle(theday(x$thedate))$values
    for (i in 1:length(thedays)){
        j <- which(theday(out$thedate)==thedays[i])
        thelabel <- paste(newlabel,i,sep="-")
        out <- average(out,j,thelabel)
    }
    return(out)
}

averagebypos <- function(X,pos,newlabel){
    out <- X
    for (i in 1:length(pos)){
        j <- which(out$pos==pos[i])
        thelabel <- paste(newlabel,i,sep="-")
        out <- average(out,j,thelabel)
    }
    return(out)
}

#' Merge a list of logratio data
#'
#' Recursively concatenates a list of logratio data into one big dataset
#' 
#' @param lrlist a list containing items of class
#' \code{\link{logratios}} or \code{\link{redux}}
#' @return an object of the same class as \code{x} containing the
#' merged dataset
#' @examples
#' samplefile <-  system.file("Samples.csv",package="ArArRedux")
#' kfile <- system.file("K-glass.csv",package="ArArRedux")
#' cafile <- system.file("Ca-salt.csv",package="ArArRedux")
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' blanklabel <- "EXB#"
#' Jpos <- c(3,15)
#' dlabels <- c("H1","AX","L1","L2")
#'  
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' mk <- loaddata(kfile,masses) # K-interference data
#' mca <- loaddata(cafile,masses) # Ca interference data
#' md <- loaddata(dfile,dlabels,PH=TRUE) # detector intercalibrations
#'  
#' # form and fit logratios
#' l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
#' lk <- fitlogratios(blankcorr(mk,blanklabel,"K:"),"Ar40")
#' k <- getmasses(lk,"Ar39","Ar40") # subset on the relevant isotopes
#' lca <- fitlogratios(blankcorr(mca,blanklabel,"Ca:"),"Ar37")
#' ca <- getmasses(lca,c("Ar36","Ar39"),c("Ar37","Ar37")) # subset
#' ld <- fitlogratios(blankcorr(md))
#' d <- averagebyday(ld,"DCAL")
#' 
#' # merge all data (except air shots) into one big logratio structure
#' X <- newredux(concat(list(l,k,ca,d)),Jpos)
#' data(Melbourne)
#' if (isTRUE(all.equal(Melbourne$X,X))) {
#'    print("We just reconstructed the built-in dataset Melbourne$X")}
#' @export
concat <- function(lrlist){
    if (length(lrlist)==2) {
        x <- lrlist[[1]]
        y <- lrlist[[2]]
        out <- x
        out$irr <- c(x$irr,y$irr)
        out$pos <- c(x$pos,y$pos)
        out$labels <- c(x$labels,y$labels)
        out$thedate <- c(x$thedate,y$thedate)
        out$num <- c(x$num,y$num)
        out$den <- c(x$den,y$den)
        out$nlr <- c(x$nlr,y$nlr)
        out$intercepts <- c(x$intercepts,y$intercepts)
        out$covmat <- mergematrices(x$covmat,y$covmat)
    } else {
        x <- lrlist[[1]]
        rest <- lrlist[2:length(lrlist)]
        y <- concat(rest)
        out <- concat(list(x,y))
    }
    return(out)
}

#' define the interference corrections
#'
#' create a new object of class \code{logratios} containing
#' the interferences from neutron reactions on Ca and K
#'
#' @param intercepts a vector with logratios
#' @param covmat the covariance matrix of the logratios
#' @param num a vector of strings marking the numerator isotopes of
#' \code{intercepts}
#' @param den a vector of strings marking the denominator isotopes of
#' \code{intercepts}
#' @param irr an object of class \code{irradiations}
#' @param label a string with a name which can be used to identify the
#' interference data in subsequent calculations
#' @return an object of class \code{logratios}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' irrfile <- system.file("irradiations.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' X <- read(samplefile,masses,blabel="EXB#",Jpos=c(3,15))
#' irr <- loadirradiations(irrfile)
#'# assume log(36Ar/37Ar) = log(39Ar/37Ar) = 1 in co-irradiate Ca-salt
#'# with variances of 0.0001 and zero covariances
#' ca <- interference(intercepts=c(1,1),
#'                    covmat=matrix(c(0.001,0,0,0.001),nrow=2),
#'                    num=c("Ar39","Ar36"),den=c("Ar37","Ar37"),
#'                    irr=X$irr[1],label="Ca-salt")
#'# assume log(39Ar/40Ar) = 4.637788 in co-irradiate K-glass
#'# with variance 7.9817e-4
#' k <- interference(intercepts=4.637788,covmat=7.9817e-4,
#'                   num="Ar39",den="Ar40",irr=X$irr[1],
#'                   label="K-glass")
#' ages <- process(X,irr,ca=ca,k=k)
#' summary(ages)
#' @export
interference <- function(intercepts,covmat,num,den,irr,label){
    out <- list()
    class(out) <- "logratios"
    out$irr <- irr
    out$pos <- 0
    out$labels <- label
    out$thedate <- as.numeric(as.Date("1970-01-01 00:00:00"))
    out$num <- num
    out$den <- den
    out$nlr <- length(intercepts)
    out$intercepts <- intercepts
    out$covmat <- covmat
    return(out)
}

mergematrices <- function(xcovmat,ycovmat){
    if (is.null(xcovmat)) return(ycovmat)
    if (is.null(ycovmat)) return(xcovmat)
    nx <- max(1,nrow(xcovmat))
    ny <- max(1,nrow(ycovmat))
    covmat <- matrix(0,nrow=nx+ny,ncol=nx+ny)
    covmat[1:nx,1:nx] <- xcovmat
    covmat[(nx+1):(nx+ny),(nx+1):(nx+ny)] <- ycovmat
    return(covmat)
}

# X = object of class "redux"
Jcal <- function(X,clabel,detectors){
    ci <- array(grep(clabel, X$labels)) # calibration indices
    si <- array(grep(clabel, X$labels, invert=TRUE))
    nc <- length(ci) # number of calibrations
    ns <- length(si)
    nci <- sum(X$nlr[ci]) # number of calibration intercepts
    nsi <- sum(X$nlr[-ci]) # number of sample intercepts
    J <- matrix(0,nrow=nsi,ncol=nsi+nci)
    for (i in si){ # loop through all the samples
        # indices of the nearest calibration data
        ic <- ci[nearest(X$thedate[i],X$thedate[ci])]
        # intercept indices of current sample
        sj <- getindices(X,prefix=X$labels[i])
        for (js in sj){ # loop through all the masses
            num <- X$num[js]
            den <- X$den[js]
            # get the detectors for the numerator and denominator masses
            ndet <- detectors[[num]]
            ddet <- detectors[[den]]
            if (ndet != ddet) {
                jn <- getindices(X,prefix=X$labels[ic],num=ndet)
                jd <- getindices(X,prefix=X$labels[ic],num=ddet)
                J[js,jn] <- -1
                J[js,jd] <- 1
            }
            J[js,js] <- 1
        }
    }
    return(J)
}

#' Detector calibration
#'
#' Apply the detector calibration for multicollector data
#' 
#' @param X an object of class \code{redux}
#' @param clabel the label of the detector calibration data
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' plotcorr(C)
#' @export
calibration <- function(X,clabel){
    j <- grep(clabel,X$labels,invert=TRUE)
    if (length(j)==length(X$labels)) return(X)
    J <- Jcal(X,clabel,X$detectors)
    out <- subset(X,j)
    out$intercepts <- J %*% X$intercepts
    out$covmat <- J %*% X$covmat %*% t(J)
    return(out)
}

Jair <- function(X,detectors){
    ns <- nruns(X)
    di <- getrunindices(X,detectors)
    ai <- getrunindices(X,"air-ratio")
    dai <- c(di,ai)
    si <- (1:ns)[-dai] # sample indices
    ndai <- sum(X$nlr[dai])
    nsi <- sum(X$nlr[si])
    J <- matrix(0,nsi,nsi+ndai)
    for (i in si){ # loop through the samples
        # intercept indices of current sample
        sj <- getindices(X,prefix=X$labels[i])
        for (js in sj){ # loop through all the masses
            num <- X$num[js]
            den <- X$den[js]
            a <- as.integer(substr(num,start=3,stop=6))
            b <- as.integer(substr(den,start=3,stop=6))
            # indices of the nearest calibration data
            dj <- array(grep(detectors[[den]],X$labels))
            id <- di[nearest(X$thedate[i],X$thedate[dj])]
            jd <- getindices(X,prefix=X$labels[id])
            ja <- getindices(X,prefix="air-ratio")
            J[js,js] <- 1
            J[js,jd] <- (log(a)-log(b))/(log(40)-log(36))
            J[js,ja] <- (log(a)-log(b))/(log(40)-log(36))
        }
    }
    return(J)
}

air <- function(X){
    out <- list(
        num = "Ar40",
        den = "Ar36",
        intercepts = log(X$param$air), 
        covmat = (X$param$sair/X$param$air)^2, # variance of the air ratio
        irr = NULL,
        pos = NULL,
        labels = "air-ratio",
        thedate = as.numeric(as.Date("1970-01-01 00:00:00")),
        nlr = 1
    )
    return(out)
}

#' Apply the mass fractionation correction
#'
#' Applies the fractionation obtained from air shot data by
#' \code{\link{fractionation}} to the denominator detector in order to
#' correct it for the mass difference between the numerator and
#' denominator isotopes.
#' 
#' @param X an object of class \code{redux}
#' @param fract a list with fractionation data for the detectors used
#' to measure the denominator isotopes
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' A <- massfractionation(C,Melbourne$fract)
#' plotcorr(A)
#' @export
massfractionation <- function(X,fract){
    fdet <- X$detectors[c("Ar37","Ar40")]
    # add the air shot data
    errorlessair <- air(X)
    errorlessair$covmat <- 0
    if (methods::is(fract,"logratios")){
        Y <- concat(list(X,fract,errorlessair))
    } else {
        Y <- concat(c(list(X),fract,list(errorlessair)))
    }
    # apply the fractionation correction
    J <- Jair(Y,fdet)
    j <- getrunindices(Y,c(unlist(fdet),"air-ratio"),invert=TRUE)
    out <- subset(Y,j)
    out$intercepts <- J %*% Y$intercepts
    out$covmat <- J %*% Y$covmat %*% t(J)
    return(out)
}

#' Load the irradiation schedule
#'
#' Loads a .csv file with the schedule of a multi-stage neutron
#' irradiation
#' 
#' @param fname file name (in .csv format)
#' @return a list of irradiations, where each irradiation is a named
#' list containing:
#' 
#' \code{tin}: vector with the start times of irradiations \cr
#' \code{tout}: vector with the end times of irradiations \cr
#' \code{P}: vector with the power of the irradiations
#' @examples
#' irrfile <- system.file("irradiations.csv",package="ArArRedux")
#' irr <- loadirradiations(irrfile)
#' str(irr)
#' @export
loadirradiations <- function(fname){
    f <- file(fname)
    open(f)
    out <- list()
    while (TRUE) { # read the file line by line
        line <- readLines(f, n=1, warn=FALSE)
        if (length(line) <= 0){ # EOF
            close(f)
            return(out)
        }
        l <- unlist(strsplit(line, split=','))
        if (l[1]=='In') {
            out[[irrname]]$tin <- c(out[[irrname]]$tin,readthedate(l[2]))
            out[[irrname]]$P <- c(out[[irrname]]$P,as.numeric(l[3]))
        } else if (l[1]=='Out'){
            out[[irrname]]$tout <- c(out[[irrname]]$tout,readthedate(l[2]))
        } else {
            irrname <- l[1]
            out[[irrname]] <- list(P=c(),tin=c(),tout=c())
        }
    }
}

getD <- function(P,T,dt,lambda){
    num <- sum(P*dt)
    den <- sum(P*(exp(-lambda*T)-exp(-lambda*(T+dt)))/lambda)
    return(log(num)-log(den))
}

getdDdL <- function(P,T,dt,lambda){
    num <- sum(P*((lambda*T-lambda*dt+1)*exp(-lambda*(T+dt)) -
                  (lambda*T          +1)*exp(-lambda*T)))
    den <- sum(P*(exp(-lambda*(T+dt))-exp(-lambda*T)))
    return(num/den)
}

Jdecay <- function(X,isotope){
    ni <- length(X$intercepts)
    ns <- nruns(X)
    J <- matrix(0,nrow=ni,ncol=ni+ns)
    ii <- 0
    for (i in 1:ns){
        for (j in 1:X$nlr[i]){
            ii <- ii + 1
            J[ii,ii] <- 1
            if (!is.na(X$num[ii]) & X$num[ii]==isotope) J[ii,ni+i] <- 1
            if (!is.na(X$den[ii]) & X$den[ii]==isotope) J[ii,ni+i] <- -1
        }
    }
    return(J)
}

getTdt <- function(irr,thedate){
    dt <- (irr$tout-irr$tin)/(365*24*3600)
    T <- (thedate-irr$tout)/(365*24*3600)
    return(list(T=T,dt=dt))
}

# computes decay correction
getDmatrix <- function(X,irradiations,isotope){
    ns <- nruns(X)
    D <- rep(0,ns)
    dDdL <- rep(0,ns)
    if (isotope=="Ar37"){
        lambda <- X$param$l7
        vlambda <- X$param$sl7^2
    }
    if (isotope=="Ar39"){
        lambda <- X$param$l9
        vlambda <- X$param$sl9^2
    }
    for (i in 1:ns){ # loop through the samples
        irr <- irradiations[[X$irr[i]]]
        Tdt <- getTdt(irr,X$thedate[i])
        D[i] <- getD(irr$P,Tdt$T,Tdt$dt,lambda)
        dDdL[i] <- getdDdL(irr$P,Tdt$T,Tdt$dt,lambda)
    }
    covmatD <- (dDdL * vlambda) %*% t(dDdL)
    return(list(intercepts=D,covmat=covmatD))
}

#' Correct for radioactive decay occurred since irradiation
#'
#' Correct for radioactive decay of neutron-induced 37Ar and 39Ar
#' occurred since irradiation
#' 
#' @param X an objects of class \code{redux}
#' @param irr the irradiation schedule
#' @param isotope a string denoting the isotope that needs correcting
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' A <- massfractionation(C,Melbourne$fract)
#' D9 <- decaycorrection(A,Melbourne$irr,"Ar39")
#' plotcorr(D9)
#' @export
decaycorrection <- function(X,irr,isotope){
    out <- X
    D <- getDmatrix(X,irr,isotope)
    J <- Jdecay(X,isotope)
    intercepts <- c(X$intercepts,D$intercepts)
    covmat <- mergematrices(X$covmat,D$covmat)
    out$intercepts <- J %*% intercepts
    out$covmat <- J %*% covmat %*% t(J)
    return(out)
}

getE <- function(P,T,dt,lambda){
    if (any(T<0)) return(0)
    num <- sum(P*(exp(-lambda*(T+dt))-exp(-lambda*T)))
    den <- lambda*sum(P*dt)
    return(log(1+num/den))
}

getdEdL <- function(P,T,dt,lambda){
    if (any(T<0)) return(0)
    num <- sum(P*((lambda*T          +1)*exp(-lambda* T    ) -
                  (lambda*T+lambda*dt+1)*exp(-lambda*(T+dt))))
    den <- lambda*sum(P*(lambda*dt + exp(-lambda*(T+dt)) - exp(-lambda*T)))
    return(num/den)
}

# computes Cl correction
getEmatrix <- function(X,irradiations){
    out <- X
    ns <- nruns(X)
    lpcl <- log(X$param$pcl)
    slpcl <- X$param$spcl/X$param$pcl
    out$intercepts <- rep(lpcl,ns)
    out$num <- rep("Ar36",ns) 
    out$den <- rep("Ar38",ns)
    out$nlr <- rep(1,ns)
    dEdL <- rep(0,ns)
    lambda <- X$param$l6
    for (i in 1:ns){ # loop through the samples
        irr <- irradiations[[X$irr[i]]]
        Tdt <- getTdt(irr,X$thedate[i])
        out$intercepts[i] <- out$intercepts[i] +
                             getE(irr$P,Tdt$T,Tdt$dt,lambda)
        dEdL[i] <- getdEdL(irr$P,Tdt$T,Tdt$dt,lambda)
        out$labels[i] <- paste("Cl:",X$labels[i],sep='')
    }
    J <- cbind(rep(1,ns),dEdL)
    covmatE <- matrix(c(slpcl^2,0,0,X$param$sl6^2),nrow=2)
    out$covmat <- J %*% covmatE %*% t(J)
    return(out)
}

#' Cl-interference correction
#'
#' Apply the interference correction for the Cl-decay products
#' 
#' @param X an object of class \code{redux}
#' @param irr the irradiation schedule
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' Cl <- clcorrection(Melbourne$X,Melbourne$irr)
#' plotcorr(Cl)
#' @export
clcorrection <- function(X,irr){
    E <- getEmatrix(X,irr)
    return(concat(list(X,E)))
}

expired <- function(irr,thedate,l7){
    Tdt <- getTdt(irr,thedate)
    return(any(Tdt$T>5*log(2)/l7))
}

getJabcdef <- function(Z,Slabels,nl){
    J <- matrix(0,nrow=nl,ncol=length(Z$intercepts))
    for (i in 1:length(Slabels)){
        j <- (i-1)*5
        label <- Slabels[i]
        iair <- getindices(Z,"air-ratio","Ar40","Ar36")
        i60 <- getindices(Z,label,"Ar36","Ar40")
        i70 <- getindices(Z,label,"Ar37","Ar40")
        i80 <- getindices(Z,label,"Ar38","Ar40")
        i90 <- getindices(Z,label,"Ar39","Ar40")
        i67ca <- getindices(Z,"Ca-salt","Ar36","Ar37")
        i97ca <- getindices(Z,"Ca-salt","Ar39","Ar37")
        i68cl <- getindices(Z,paste("Cl:",label,sep=""),"Ar36","Ar38")
        J[j+1,c(iair,i60)] <- 1
        J[j+2,c(iair,i67ca,i70)] <- 1
        J[j+3,c(iair,i68cl,i80)] <- 1
        J[j+4,i90] <- 1
        J[j+5,c(i97ca,i70)] <- 1
    }
    i90k <- getindices(Z,"K-glass","Ar39","Ar40")
    J[nl,i90k] <- -1
    return(J)
}

getabcdef <- function(Cl){
    Z <- concat(list(Cl,air(Cl))) # matrix with everything
    si <- getrunindices(Z,c("Ca-salt","K-glass","Cl:","air-ratio"),
                        invert=TRUE)
    ns <- length(si)
    theS <- subset(Z,si) # contains only samples
    theK <- subset(Z,labels="K-glass") # contains only K-glass
    out <- Z
    out$irr <- c(theS$irr,theK$irr)
    out$pos <- c(theS$pos,theK$pos)
    out$labels <- c(theS$labels,theK$labels)
    out$num <- c(rep(c("a","b","c","d","e"),ns),"f")
    out$den <- rep(NA,5*ns+1)
    out$nlr <- c(rep(5,ns),1)
    out$thedate <- c(theS$thedate,theK$thedate)
    Jv <- getJabcdef(Z,theS$labels,length(out$num))
    out$intercepts <- exp(Jv %*% Z$intercepts)
    Jw <- apply(Jv,2,"*",out$intercepts)
    out$covmat <- Jw %*% Z$covmat %*% t(Jw)
    return(out)
}

#' Calculate the 40Ar*/39ArK-ratios
#'
#' Calculate the 40Ar*/39ArK-ratios of interference corrected logratio
#' intercept data
#' 
#' @param X an object of class \code{redux} containing some
#' interference corrected logratio intercept data
#' @param irr the irradiation schedule
#' @return an object of class \code{link{redux}} containing the
#' 40Ar*/39ArK-ratios as \code{intercepts} and its covariance matrix
#' as \code{covmat}
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' plotcorr(R)
#' @export
get4039 <- function(X,irr){
    Y <- getabcdef(X)
    ns <- nruns(Y)-1
    ni <- length(Y$intercepts)
    out <- subset(Y,1:ns)
    out$intercepts <- rep(0,ns)
    out$covmat <- matrix(0,ns,ns)
    out$num <- rep(NA,ns)
    out$den <- rep(NA,ns)
    out$nlr <- rep(1,ns)
    J <- matrix(0,nrow=ns,ncol=ni)
    hasKglass <- "K-glass" %in% X$labels
    hasCasalt <- "Ca-salt" %in% X$labels
    if (hasKglass) { ff <- utils::tail(Y$intercepts,n=1) }
    else { ff <- 0 }
    for (i in 1:ns){
        j <- (i-1)*5
        label <- Y$labels[i]
        aa <- Y$intercepts[getindices(Y,label,num='a')]
        bb <- Y$intercepts[getindices(Y,label,num='b')]
        cc <- Y$intercepts[getindices(Y,label,num='c')]
        dd <- Y$intercepts[getindices(Y,label,num='d')]
        ee <- Y$intercepts[getindices(Y,label,num='e')]
        if (!hasCasalt | expired(irr[[Y$irr[i]]],Y$thedate[i],Y$param$l7)){
            out$intercepts[i] <- (1-aa+cc)/dd-ff
            J[i,j+1] <- -1/dd           # dR/da
            J[i,j+2] <-  0              # dR/db
            J[i,j+3] <-  1/dd           # dR/dc
            J[i,j+4] <- -(1-aa+cc)/dd^2 # dR/dd
            J[i,j+5] <-  0              # dR/de
        } else {
            out$intercepts[i] <- (1-aa+bb+cc)/(dd-ee)-ff
            J[i,j+1] <- -1/(dd-ee)
            J[i,j+2] <-  1/(dd-ee)
            J[i,j+3] <-  1/(dd-ee)
            J[i,j+4] <- -(1-aa+bb+cc)/(dd-ee)^2
            J[i,j+5] <-  (1-aa+bb+cc)/(dd-ee)^2
        }
    }
    if (hasKglass) J[,ni] <- -1
    out$covmat <- J %*% Y$covmat %*% t(J)
    return(out)
}

interpolateRJ <- function(R){ 
    S <- averagebypos(R,R$Jpos,newlabel="J")
    ji <- getindices(S,pos=R$Jpos)
    si <- getindices(S,pos=R$Jpos,invert=TRUE)
    ns <- length(si)
    J <- matrix(0,nrow=2*ns,ncol=length(S$intercepts))
    for (i in 1:length(si)){ # loop through samples
        J[i,si[i]] <- 1
        px <- S$pos[si[i]]
        dx <- px-R$Jpos
        pm <- R$Jpos[which(dx==max(dx[dx<0]))]
        pM <- R$Jpos[which(dx==min(dx[dx>0]))]
        mi <- getindices(S,pos=pm)
        Mi <- getindices(S,pos=pM)
        J[ns+i,mi] <- (px-pm)/(pM-pm)
        J[ns+i,Mi] <- (pM-px)/(pM-pm)
    }
    out <- R
    samps <- subset(S,si)
    out$thedate <- rep(samps$thedate,2)
    out$irr <- rep(samps$irr,2)
    out$pos <- rep(samps$pos,2)
    out$num <- rep(samps$num,2)
    out$den <- rep(samps$den,2)
    out$nlr <- rep(samps$nlr,2)
    out$labels <- c(samps$labels,paste('J:',samps$labels,sep=''))
    out$intercepts <- J %*% S$intercepts
    out$covmat <- J %*% S$covmat %*% t(J)    
    return(out)
}

getJJ <- function(RS,ns,lambda,ts){
    J <- matrix(0,nrow=2*ns+1,ncol=2*ns+2)
    J[1:ns,1:ns] <- diag(ns)
    J[2*ns+1,2*ns+1] <- 1
    for (i in (ns+1):(2*ns)){
        J[i,i] <- (1-exp(lambda*ts))/(RS[i]^2)
        J[i,2*ns+1] <- ts*exp(lambda*ts)/RS[i]
        J[i,2*ns+2] <- lambda*exp(lambda*ts)/RS[i]
    }
    return(J)
}

#' Calculate the irradiation parameter ('J factor')
#'
#' Interpolate the irradiation parameters for the samples
#' given the 40Ar*/39ArK ratios of the samples and fluence monitors
#' 
#' @param R a vector of 40Ar*/39ArK ratios
#' @return an object of class \code{redux} containing, as
#' \code{intercepts}, the 40Ar*/39ArK ratios of the samples, the
#' interpolated J-factors, and the 40K decay constant; and as
#' \code{covmat}: the covariance matrix. All other class properties
#' are inherited from \code{R}.
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' J <- getJfactors(R)
#' plotcorr(J)
#' @export
getJfactors <- function(R){
    RS <- interpolateRJ(R)
    lambda <- R$param$l0
    ts <- R$param$ts
    covmat <- mergematrices(RS$covmat,
              matrix(c(R$param$sl0^2,0,0,R$param$sts^2),nrow=2))
    ns <- length(RS$labels)/2
    out <- RS
    out$intercepts[(ns+1):(2*ns)] <-
        (exp(lambda*ts)-1)/RS$intercepts[(ns+1):(2*ns)]
    out$intercepts[2*ns+1] <- lambda
    out$labels[2*ns+1] <- 'lambda'
    out$thedate[2*ns+1] <- NA
    J <- getJJ(RS$intercepts,ns,lambda,ts)
    out$covmat <- J %*% covmat %*% t(J)
    return(out)
}

#' Calculate 40Ar/39Ar ages
#'
#' Calculate 40Ar/39Ar ages from a vector of 40Ar/39Ar-ratios and
#' J-factors
#' 
#' @param RJ an object of class \code{Redux}
#' containing the amalgamated $^{40}$Ar$^*$/$^{39}$Ar$_K$-ratios and
#' J-factors with their covariance matrix
#' @return an object of class \code{results} containing the
#' ages and their covariance matrix
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' J <- getJfactors(R)
#' ages <- getages(J)
#' plotcorr(ages)
#' @export
getages <- function(RJ){
    out <- list()
    class(out) <- "results"
    ns <- (length(RJ$intercepts)-1)/2
    lambda <- utils::tail(RJ$intercepts,n=1)
    out$thedate <- RJ$thedate[1:ns]
    out$labels <- RJ$labels[1:ns]
    out$ages <- log(1+RJ$intercepts[1:ns]*
                    RJ$intercepts[(ns+1):(2*ns)])/lambda
    J <- matrix(0,nrow=ns,ncol=2*ns+1)
    for (i in 1:ns){
        r <- RJ$intercepts[i]
        j <- RJ$intercepts[ns+i]
        J[i,i] <- j/(lambda*(1+j*r))
        J[i,ns+i] <- r/(lambda*(1+j*r))
        J[i,2*ns+1] <- - log(1+j*r)/(lambda^2)
    }
    out$covmat <- J %*% RJ$covmat %*% t(J)
    return(out)
}


#' Calculate the weighted mean age
#'
#' Computes the error weighted mean and MSWD of some samples taking
#' into covariances.
#' 
#' @param ages an object of class \code{results}
#' @param prefix is either a string with the prefix of the
#' samples that need to be averaged, or a vector of sample names.
#' @return a list with items:
#'  
#' \code{avgt}: the weighted mean age\cr
#' \code{err}: the standard error of \code{avgt}\cr
#' \code{MSWD}: the Mean Square of the Weighted Deviates
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' weightedmean(ages,"MD2-1")
#' @export
weightedmean <- function(ages,prefix=NULL){
    if (is.null(prefix)){
        slabs <- ages$labels
    } else if (length(prefix)==1){
        slabs <- ages$labels[grep(prefix,ages$labels)]
    } else {
        slabs <- prefix
    }
    subs <- subset(ages,labels=slabs)
    n <- nruns(subs)
    W <- rep(1,n) # design matrix
    invSigma <- solve(subs$covmat)
    vart <- 1/(W %*% invSigma %*% W)
    avgt <- vart*(W %*% invSigma %*% subs$ages)
    mswd <- ((subs$ages - avgt) %*% solve(subs$covmat) %*%
             (subs$ages - avgt))/(n-1)    
    return(list(avgt=avgt,err=sqrt(vart),MSWD=mswd))
}

#' Compute the mass fractionation correction
#'
#' Compares the measured 40Ar/36Ar ratio of an air shot on a given
#' detector with the atmospheric ratio.
#' 
#' @param fname a .csv file with the air shot data
#' @param detector the name of the ion detector
#' @param MS the type of mass spectrometer
#' @param PH TRUE if the data were recorded in 'peak hopping' mode,
#' FALSE if recorded in multicollector mode.
#' @examples
#' data(Melbourne)
#' fd37file <- system.file("AirL2.csv",package="ArArRedux")
#' fd40file <- system.file("AirH1.csv",package="ArArRedux")
#' fract <- list(fractionation(fd37file,"L2",PH=TRUE),
#'               fractionation(fd40file,"H1",PH=FALSE))
#' if (isTRUE(all.equal(Melbourne$fract,fract))){
#'   print("We just re-created the fractionation correction for the Melbourne dataset")
#' }
#' @return an object of class \code{\link{logratios}}
#' @export
fractionation <- function(fname,detector,MS="ARGUS-VI",PH=FALSE){
    mf <- loaddata(fname,c("Ar40","Ar36"),MS,PH)
    lf <- fitlogratios(blankcorr(mf),"Ar40")
    f <- averagebyday(lf,detector)
    return(f)
}

#' Summary table
#'
#' Plots the ages and their standard errors
#' 
#' @param object an objct of class \code{\link{results}}
#' @param ... no other arguments
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' summary(ages)[1:5,]
#' @rdname summary
#' @method summary results
#' @export
summary.results <- function(object,...){
    tab <- cbind(object$labels,object$ages,sqrt(diag(object$covmat)))
    colnames(tab) <- c("name","age[Ma]","s.e.[Ma]")
    print(tab)
    return(tab)
}

#' Plot a matrix with correlation coefficients
#'
#' Converts the covariance matrix to a correlation matrix and plots
#' this is a coloured image for visual inspection.
#' 
#' @param X a data structure (list) containing an item called `covmat' (covariance matrix)
#' @examples
#' data(Melbourne)
#' plotcorr(Melbourne$X)
#' @export
plotcorr <- function(X){
    image.with.legend(z=stats::cov2cor(X$covmat),color.palette=grDevices::heat.colors)
}

# modified version of filled.contour with ".filled.contour" part replaced with "image"
# function. Note that the color palette is a flipped heat.colors rather than cm.colors
image.with.legend <- function (x = seq(1, nrow(z), length.out = nrow(z)), y = seq(1, 
    ncol(z), length.out=nrow(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = grDevices::heat.colors, 
    col = rev(color.palette(length(levels) - 1)), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) {
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(1, nrow(z), length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- graphics::par(c("mar", "las", "mfrow")))$mar
    on.exit(graphics::par(par.orig))
    w <- (3 + mar.orig[2L]) * graphics::par("csi") * 2.54
    graphics::layout(matrix(c(2, 1), ncol = 2L), widths = c(1, graphics::lcm(w)))
    graphics::par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    graphics::par(mar = mar)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = range(levels),
                xaxs = "i", yaxs = "i")
    graphics::rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            graphics::axis(4)
    }
    else key.axes
    graphics::box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    graphics::par(mar = mar)
    graphics::image(x,y,z,col=col,xlab="",ylab="")
    if (missing(plot.axes)) {
        if (axes) {
            graphics::title(main = "", xlab = "", ylab = "")
            graphics::Axis(x, side = 1)
            graphics::Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        graphics::box()
    if (missing(plot.title)) 
        graphics::title(...)
    else plot.title
    invisible()
}

#' Process logratio data and calculate 40Ar/39Ar ages
#'
#' Performs detector calibration, mass fractionation correction, decay
#' corrections, interference corrections, interpolates J-factors and
#' calculates ages.
#'
#' @param X an object of class \code{\link{redux}}
#' @param irr the irradiation schedule
#' @param fract list with air shot data (one item per denominator isotope)
#' @param ca an object of class \code{\link{logratios}} with
#' Ca-interference data (not necessary if interferences are included in X)
#' @param k an object of class \code{\link{logratios}} with
#' K-interference data (not necessary if interferences are included in X)
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' summary(ages)
#' @export
process <- function(X,irr,fract=NULL,ca=NULL,k=NULL){
    # apply the detector calibration (this won't affect the Ar40/Ar36 ratio)
    C <- calibration(X,"DCAL")
    # apply the mass fractionation correction
    if (is.null(fract)){
        A <- C
    } else {
        A <- massfractionation(C,fract)
    }
    # decay corrections
    D9 <- decaycorrection(A,irr,"Ar39")
    D7 <- decaycorrection(D9,irr,"Ar37")
    if (is.null(k)){
    # interference corrections
        K <- average(D7,grep("K:",A$labels),newlabel="K-glass")
    } else {
        K <- concat(list(D7,k))
    }
    if (is.null(ca)){
        Ca <- average(K,grep("Ca:",K$labels),newlabel="Ca-salt")
    } else {
        Ca <- concat(list(K,ca))
    }
    Cl <- clcorrection(Ca,irr)
    # calculate the 40Ar*/39ArK-ratios 
    R <- get4039(Cl,irr)
    # calculate J factors
    J <- getJfactors(R)
    # calculate ages
    ages <- getages(J)
    return(ages)
}

#' Read mass spectrometer data
#'
#' Reads raw mass spectrometer data and parses it into a
#' \code{\link{redux}} format for further processing.
#' 
#' @param xfile a .csv file with samples and fluence monitor data
#' @param masses a list which specifies the order in which the isotopes
#' are recorded by the mass spectrometer
#' @param blabel a prefix string denoting the blanks
#' @param Jpos a vector of integers denoting the positions of the
#' fluence monitors in the irradiation stack
#' @param kfile a .csv file with the K-interference monitor data
#' (optional)
#' @param cafile a .csv file with the Ca-interference monitor data
#' (optional)
#' @param dfile a .csv file with the detector calibration data
#' (optional)
#' @param dlabels a list which specifies the names of the detectors
#' and the order in which they were recorded by the mass spectrometer
#' @param MS a string denoting the type of mass spectrometer. At the
#' moment only the ARGUS-IV is listed. Please contact the author to
#' add other file formats to Ar-Ar_Redux.
#' @return an object of class \code{\link{redux}}.
#' @examples
#' samplefile <-  system.file("Samples.csv",package="ArArRedux")
#' kfile <- system.file("K-glass.csv",package="ArArRedux")
#' cafile <- system.file("Ca-salt.csv",package="ArArRedux")
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' dlabels <- c("H1","AX","L1","L2")
#' X <- read(samplefile,masses,"EXB#",c(3,15),kfile,cafile,dfile,dlabels)
#' plotcorr(X)
#' @export
read <- function(xfile,masses,blabel,Jpos,kfile=NULL,cafile=NULL,
                  dfile=NULL,dlabels=NULL,MS="ARGUS-VI"){
    # load the .csv files
    m <- loaddata(xfile,masses,MS) # samples and J-standards
    x <- fitlogratios(blankcorr(m,blabel),"Ar40")    
    if (!is.null(kfile)){ # K-interference data
        # subset of the relevant isotopes
        mk <- loaddata(kfile,masses,MS)
        lk <- fitlogratios(blankcorr(mk,blabel,"K:"),"Ar40")
        k <- getmasses(lk,"Ar39","Ar40")
        x <- concat(list(x,k))
    }
    if (!is.null(cafile)){ # Ca interference data
        mca <- loaddata(cafile,masses,MS)
        lca <- fitlogratios(blankcorr(mca,blabel,"Ca:"),"Ar37")
        ca <- getmasses(lca,c("Ar36","Ar39"),c("Ar37","Ar37"))
        x <- concat(list(x,ca))
    }
    if (!is.null(dfile)){
        md <- loaddata(dfile,dlabels,MS,PH=TRUE)
        ld <- fitlogratios(blankcorr(md))
        d <- averagebyday(ld,"DCAL")
        x <- concat(list(x,d))
    }
    return(newredux(x,Jpos))
}

test <- function(builddata=FALSE,option=TRUE){

    samplefile <- "../inst/Samples.csv"
    kfile <- "../inst/K-glass.csv"
    cafile <- "../inst/Ca-salt.csv"
    fd37file <- "../inst/AirL2.csv"
    fd40file <- "../inst/AirH1.csv"
    irrfile <- "../inst/irradiations.csv"
    dfile <- "../inst/Calibration.csv"
    masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
    blanklabel <- "EXB#"
    dlabels <- c("H1","AX","L1","L2")
    Jpos <- c(3,15)

    irr <- loadirradiations(irrfile)
    
    fract <- list(fractionation(fd37file,"L2",PH=TRUE),
                  fractionation(fd40file,"H1",PH=FALSE))
    
    if (option){ # full propagation
        X <- read(samplefile,masses,blanklabel,Jpos,
                  kfile,cafile,dfile,dlabels)
        ages <- process(X,irr,fract)
    } else {

    }

    if (builddata){
        Melbourne <- list(X=X,irr=irr,fract=fract)
        save(Melbourne,file="../data/Melbourne.rda")
    }

    return(ages)

}

#graphics.off()
#rm(list=ls())
#setwd("/home/pvermees/Dropbox/Ar-Ar_Redux/ArArRedux/R")
#ages <- test(option=TRUE)
#plotcorr(ages)
#summary(ages)
#weightedmean(ages,"MD")

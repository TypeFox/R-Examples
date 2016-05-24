###
tgcd <-
function(Sigdata, npeak, inis=NULL, mdt=3,
         nstart=200, model=c("f","g"),
         elim=NULL, logy=FALSE, hr=NULL, 
         outfile=NULL, plot=TRUE)  {
    UseMethod("tgcd")
} #
### 2015.06.01, revised in 2015.09.12.
tgcd.default <- 
function(Sigdata, npeak, inis=NULL, mdt=3,
         nstart=200, model=c("f","g"),
         elim=NULL, logy=FALSE, hr=NULL, 
         outfile=NULL, plot=TRUE)  {
        ### Stop if not.
        stopifnot(ncol(Sigdata)==2L, 
                  all(Sigdata[,1L,drop=TRUE]>0),
                  all(Sigdata[,2L,drop=TRUE]>=0),
                  length(npeak)==1L, is.numeric(npeak), 
                  npeak %in% seq(13L), 3L*npeak<nrow(Sigdata), 
                  is.null(inis) || is.matrix(inis),
                  length(mdt)==1L, is.numeric(mdt), mdt>0,
                  length(nstart)==1L, is.numeric(nstart), nstart>=1L, nstart<=50000L,
                  length(model) %in% c(1L, 2L), all(model %in% c("f","g")),
                  is.null(elim) || is.numeric(elim),
                  is.logical(logy), length(logy)==1L,
                  is.null(hr) || is.numeric(hr),
                  is.null(outfile) || is.character(outfile),
                  is.logical(plot), length(plot)==1L)
        if (!is.null(inis))  {
            if(!is.numeric(inis))  stop("Error: inis should be a numeric matrix!")
            if(dim(inis)[1L]!=npeak) stop("Error: incorrect dimensions of inis!")
            if(model[1L]=="f" && dim(inis)[2L]!=3L) stop("Error: incorrect dimensions of inis!")
            if(model[1L]=="g" && dim(inis)[2L]!=4L) stop("Error: incorrect dimensions of inis!")
        } # end if.
        if (!is.null(elim)) {
            if(length(elim)!=2L) stop("Error: elim should be a two-element vector!")
            if(elim[1L]>=elim[2L]) stop("Error: invalid elim!")
            if(elim[1L]>=1.8) stop("Error: lower limit of elim is too large!")
            if(elim[2L]<=2.2) stop("Error: upper limit of elim is too small!")
        } # end if.
        if (!is.null(hr)) {
            if(length(hr)!=1L) stop("Error: hr should be a one-element vector!")
            if(hr<=0) stop("Error: hr should be larger than zero!")
        } # end if.
        if (!is.null(outfile)) {
            if(length(outfile)!=1L) stop("Error: outfile should be an one-element vector!")
        } # end if.
        ###
        temp <- as.numeric(Sigdata[,1L,drop=TRUE])
        signal <- as.numeric(Sigdata[,2L,drop=TRUE])
        ###
        if(is.null(inis))  {
            abzero <- which(signal>(.Machine$double.eps)^0.3)
            plot(temp[abzero], signal[abzero], type="l", col="skyblue3", 
                 lwd=5, xlab="Temperature (K)", log=ifelse(logy,"y",""), 
                 ylab="TL intensity (counts)", main=paste("Click the mouse to select ", 
                 npeak, " peak maxima:", sep=""))
            grid(col="orangered2", lwd=1)
            sldxy <- try(locator(n=npeak), silent=TRUE)
            if(class(sldxy)=="try-error") {
                stop("Error: no available starting values!")
            } else {
                sldxy_index <- order(sldxy$x, decreasing=FALSE)
                sldxy$x <- sldxy$x[sldxy_index]
                sldxy$y <- sldxy$y[sldxy_index]
            } # end if.
        } # end if.
        ###
        minTEMPER <- min(temp)
        maxTEMPER <- max(temp)
        minINTENS <- min(signal)
        maxINTENS <- max(signal)
        ###
        ### Function used for setting initial
        ### parameters manually (interactively).
        setpars <-function(npeak, sldx, sldy)  {
            mat1 <- mat2 <- mat3 <- mat4 <-
            as.data.frame(matrix(nrow=npeak+1L, ncol=5L))
            ###
            ### Default TL growth peak intensity. 
            mat1[1L,] <-  c("Peak", "INTENS(min)", "INTENS(max)",  
                            "INTEN(ini)", "INTENS(fix)")
            mat1[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat1[-1L,2L] <-  round(minINTENS*0.8, 3L)
            mat1[-1L,3L] <-  round(maxINTENS*1.2, 3L) 
            mat1[-1L,4L] <- if (is.null(inis)) {
                round(sldy, 3L)
            } else {
                round(inis[,1L,drop=TRUE], 3L)
            } # end if.
            mat1[-1L,5L] <- FALSE
            ###
            ### Default activation energy.
            mat2[1L,]<- c("Peak", "ENERGY(min)", "ENERGY(max)", 
                          "ENERGY(ini)", "ENERGY(fix)")
            mat2[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat2[-1L,2L] <-  if (is.null(elim)) 0.5 else round(elim[1L], 3L)
            mat2[-1L,3L] <-  if (is.null(elim)) 5.0 else round(elim[2L], 3L)
            mat2[-1L,4L] <- if (is.null(inis)) {
                round(runif(n=npeak,min=1.8, max=2.2), 3L)
            } else {
                round(inis[,2L,drop=TRUE], 3L)
            } # end if.
            mat2[-1L,5L] <- FALSE
            ###
            ### Default temperature at the peak maximum.
            mat3[1L,] <- c("Peak", "TEMPER(min)", "TEMPER(max)",  
                           "TEMPER(ini)", "TEMPER(fix)")
            mat3[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat3[-1L,2L] <-  round(minTEMPER, 3L)
            mat3[-1L,3L] <-  round(maxTEMPER, 3L)
            mat3[-1L,4L] <- if (is.null(inis)) {
                round(sldx, 3L)
            } else {
                round(inis[,3L,drop=TRUE], 3L)
            } # end if.
            mat3[-1L,5L] <- FALSE
            ###
            ### Default bValue for a glow peak.
            if (model[1L]=="g")  {
                mat4[1L,] <- c("Peak", "bValue(min)", "bValue(max)",  
                               "bValue(ini)", "bValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0
                mat4[-1L,3L] <-  2.0
                mat4[-1L,4L] <- if (is.null(inis)) {
                    round(runif(n=npeak,min=1.1, max=1.3), 3L)
                } else {
                    round(inis[,4L,drop=TRUE], 3L)
                } # end if.
                mat4[-1L,5L] <- FALSE
            } # end if.
            ###
            ###
            if (model[1L]=="f")  {
                mat <- rbind(mat1, rep("    ", 5L),
                             mat2, rep("    ", 5L), 
                             mat3)
            } else {
                mat <- rbind(mat1, rep("    ", 5L),
                             mat2, rep("    ", 5L), 
                             mat3, rep("    ", 5L),
                             mat4)
            } # end if.
            ###
            if(is.null(inis)) {
                pars <- try(edit(name=mat), silent=TRUE)
                if (class(pars)=="try-error") {
                    stop("Error: incorrect parameter modification!")
                } # end if.
            } else {
                pars <- mat
            } # end if.
            ###
            return(pars)
        } # end function setpars.
       ###
       ###
       ### cat("Set parameter constraints:\n")
       pars <- setpars(npeak=npeak, sldx=sldxy$x, sldy=sldxy$y)
       indx <- seq(from=2L, to=npeak+1L, by=1L)
       ###
       ###
       ### TL growth peak intensity.
       ### Check non-finite vlaues.
       intensity1 <- as.numeric(pars[indx,2L,drop=TRUE])
       whichloc <- which(!is.finite(intensity1))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of INTENS")
       } # end if.
       intensity2 <- as.numeric(pars[indx,3L,drop=TRUE])
       whichloc <- which(!is.finite(intensity2))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of INTENS")
       } # end if.
       intensity3 <- as.numeric(pars[indx,4L,drop=TRUE])
       whichloc <- which(!is.finite(intensity3))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of INTENS")
       } # end if.
       ###
       ### Check bounds.
       whichloc <- which(intensity3<intensity1 | 
                         intensity3>intensity2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of INTENS")
       } # end if.
       ### 
       ### Check logical values.
       fix_intensity <- pars[indx,5L,drop=TRUE]
       if (!all(fix_intensity %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of INTENS!")
       } # end if.
       fix_intensity <- as.logical(fix_intensity)
       intensity1[fix_intensity==TRUE] <- 
       intensity3[fix_intensity==TRUE]
       intensity2[fix_intensity==TRUE] <- 
       intensity3[fix_intensity==TRUE]
       ###
       ###
       ### Activation energy. 
       ### Check non-finite values.
       energy1 <- as.numeric(pars[indx+(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(energy1)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of ENERGY!")
       } # end if.
       energy2 <- as.numeric(pars[indx+(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(energy2)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of ENERGY!")
       } # end if.
       energy3 <- as.numeric(pars[indx+(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(energy3)) 
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of ENERGY!")
       } # end if.
       ###
       ### Check bounds.
       whichloc <- which(energy3<energy1 | 
                         energy3>energy2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of ENERGY")
       } # end if.
       ### Check logical values.
       fix_energy <- pars[indx+(npeak+2L),5L,drop=TRUE]
       if (!all(fix_energy %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of ENERGY!")
       } # end if.
       fix_energy <- as.logical(fix_energy)
       energy1[fix_energy==TRUE] <- 
       energy3[fix_energy==TRUE]
       energy2[fix_energy==TRUE] <- 
       energy3[fix_energy==TRUE]
       ###
       ###
       ### Temperature at the peak maximum.
       ### Check non-finite values.
       temperature1 <- as.numeric(pars[indx+2L*(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(temperature1))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite lower bound of TEMPER")
       } # end if.
       temperature2 <- as.numeric(pars[indx+2L*(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(temperature2))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite upper bound of TEMPER")
       } # end if.
       temperature3 <- as.numeric(pars[indx+2L*(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(temperature3))
       if (length(whichloc)>=1L)  {
           stop("Error: non-finite initial of TEMPER")
       } # end if.
       ### Check bounds.
       whichloc <- which(temperature3<temperature1 | 
                         temperature3>temperature2)
       if (length(whichloc)>=1L)  {
           stop("Error: unbounded initial of TEMPER")
       } # end if.
       ### Check logical values.
       fix_temperature <- pars[indx+2L*(npeak+2L),5L,drop=TRUE]
       if (!all(fix_temperature %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of TEMPER!")
       } # end if.
       fix_temperature <- as.logical(fix_temperature)
       temperature1[fix_temperature==TRUE] <- 
       temperature3[fix_temperature==TRUE]
       temperature2[fix_temperature==TRUE] <- 
       temperature3[fix_temperature==TRUE]
       ###
       ###
       ### bValue for the the glow peak.
       ### Check non-finite values.
       if (model[1L]=="g") {
           bValue1 <- as.numeric(pars[indx+3L*(npeak+2L),2L,drop=TRUE])
           whichloc <- which(!is.finite(bValue1))
           if (length(whichloc)>=1L)  {
               stop("Error: non-finite lower bound of bValue")
           } # end if.
           bValue2 <- as.numeric(pars[indx+3L*(npeak+2L),3L,drop=TRUE])
           whichloc <- which(!is.finite(bValue2))
           if (length(whichloc)>=1L)  {
               stop("Error: non-finite upper bound of bValue")
           } # end if.
           bValue3 <- as.numeric(pars[indx+3L*(npeak+2L),4L,drop=TRUE])
           whichloc <- which(!is.finite(bValue3))
           if (length(whichloc)>=1L)  {
               stop("Error: non-finite initial of bValue")
           } # end if.
           ### Check bounds.
           whichloc <- which(bValue3<bValue1 | 
                             bValue3>bValue2)
           if (length(whichloc)>=1L)  {
               stop("Error: unbounded initial of bValue")
           } # end if.
           ### Check logical values.
           fix_bValue <- pars[indx+3L*(npeak+2L),5L,drop=TRUE]
           if (!all(fix_bValue %in% c("TRUE", "FALSE", "T", "F")))  {
               stop("Error: non-logical variable in the 5th column of bValue!")
           } # end if.
           fix_bValue <- as.logical(fix_bValue)
           bValue1[fix_bValue==TRUE] <- 
           bValue3[fix_bValue==TRUE]
           bValue2[fix_bValue==TRUE] <- 
           bValue3[fix_bValue==TRUE]
       } # end if.
       ###
       ###
       nd <- length(temp)
       n2 <- ifelse(model[1L]=="f", 3L*npeak, 4L*npeak)
       fmin <- message <- 0
       ###
       if(model[1L]=="f")  {
           lower <- c(intensity1,energy1, temperature1)
           upper <- c(intensity2,energy2, temperature2)
           pars <- c(intensity3,energy3, temperature3)
       } else {
           lower <- c(intensity1, energy1, temperature1, bValue1)
           upper <- c(intensity2, energy2, temperature2, bValue2)
           pars <- c(intensity3, energy3, temperature3, bValue3)
       } # end if.
       ###
       subroutine <- ifelse(model[1L]=="f", "tgcd", "tgcd1")
       res <- .Fortran(subroutine, as.numeric(temp), as.numeric(signal),
                       as.integer(nd), pars=as.double(pars), as.integer(n2), 
                       fmin=as.double(fmin), message=as.integer(message), 
                       as.double(lower), as.double(upper), as.integer(nstart),
                       as.double(mdt), PACKAGE="tgcd")
       if (res$message!=0) {
           stop("Error: fail in glow curve deconvolution!")
       } # end if.
       ###
       pars <- matrix(res$pars, ncol=ifelse(model[1L]=="f",3L,4L))
       index <- order(pars[,3L,drop=TRUE], decreasing=FALSE)
       pars <- pars[index,,drop=FALSE]
       colnames(pars) <- if (model[1L]=="f") {
           c("INTENS", "ENERGY", "TEMPER")
       } else {
           c("INTENS", "ENERGY", "TEMPER", "bValue")
       } # end if.
       rownames(pars) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       kbz <- 8.617385e-5
       ###
       CompSig <- matrix(nrow=nd, ncol=npeak)
       if (model[1L]=="f")  {
           alpha <- function (x)  {
               a0 <- 0.267773734; a1 <- 8.6347608925
               a2 <- 18.059016973; a3 <- 8.5733287401
               b0 <- 3.9584969228; b1 <- 21.0996530827
               b2 <- 25.6329561486; b3 <- 9.5733223454
               return(1.0-(a0+a1*x+a2*x^2L+a3*x^3L+x^4L)/
                     (b0+b1*x+b2*x^2L+b3*x^3L+x^4L) )
           } # end function alpha.
           ###
           for(i in seq(npeak)) {
               maxi <- pars[i,1L]
               engy <- pars[i,2L]
               maxt <- pars[i,3L]
               xa <- engy/kbz/maxt
               xb <- engy/kbz/temp
               CompSig[,i] <- maxi*exp(xa-xb)*exp(xa*(alpha(xa)-
                              temp/maxt*alpha(xb)*exp(xa-xb)))
           } # end for. 
       } else {
           for(i in seq(npeak)) {
               maxi <- pars[i,1L]
               engy <- pars[i,2L]
               maxt <- pars[i,3L]
               bv <- pars[i,4L]
               xa <- 2.0*kbz*temp/engy
               xb <- 2.0*kbz*maxt/engy 
               expv <- exp(engy/kbz/temp*(temp-maxt)/maxt)
               CompSig[,i] <- maxi*(bv^(bv/(bv-1.0)))*expv*
                 ((bv-1.0)*(1.0-xa)*((temp/maxt)^2L)*expv+1.0+(bv-1.0)*xb)^(-bv/(bv-1.0))
           } # end for. 
       } # end if.
       rowsumSig <- rowSums(CompSig)
       FOM <- res$fmin/sum(rowsumSig)*100
       ###
       calShape <- function(y, x)  {
           ny <- length(y)
           maxloc <- which.max(y)
           hmaxval <- max(y)/2.0
           Tm <- x[maxloc]
           T1 <- approx(x=y[1L:maxloc], y=x[1L:maxloc], xout=hmaxval)$y
           T2 <- approx(x=y[maxloc:ny], y=x[maxloc:ny], xout=hmaxval)$y
           d1 <- Tm-T1
           d2 <- T2-Tm
           thw <- T2-T1
           sf <- d2/thw
           return( c("T1"=T1, "T2"=T2, "Tm"=Tm, "d1"=d1, 
                     "d2"=d2, "thw"=thw, "sf"=sf) )
       } # end function calShape.
       ###
       sp <- t(apply(CompSig, MARGIN=2L, calShape, temp))
       rownames(sp) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       if (!is.null(hr))  {
           calff1 <- function(et)  {
               energy <- et[1L]
               temper <- et[2L]
               return(hr*energy/kbz/temper^2L*exp(energy/kbz/temper))             
           } # end function calff1.
           ###
           calff2 <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               bv <- et[3L]
               return(hr*energy/kbz/temper^2L*exp(energy/kbz/temper)/
                      (1.0+(bv-1.0)*(2.0*kbz*temper)/energy))
           } # end function calff2.
           ###
           if (model[1L]=="f")  {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff1)
           } else {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff2)
           } # end if.
           ###
           output <-list("pars"=pars, "ff"=ff, "sp"=sp, "FOM"=FOM)
       } else {
           output <-list("pars"=pars, "sp"=sp, "FOM"=FOM)
       } # end if
       ###
       if (plot==TRUE) {
           layout(cbind(c(1L,1L,1L,2L),c(1L,1L,1L,2L)))
           par(mar=c(0,5.1,3.1,1.1))
           lineCol <- c("deepskyblue", "orangered", "purple", 
                        "violetred", "yellowgreen", "lightblue", 
                        "goldenrod", "forestgreen", "blue", 
                        "plum", "tan", "violet", "grey50")
           plot(temp, signal, type="p", pch=21, bg="white", cex=0.8,
                ylab="TL intensity (counts)", las=0, lab=c(7,7,9), xaxt="n", 
                xaxs="r", yaxs="i", cex.lab=1.5)
           box(lwd=2L)
           XaxisCentral <- median(axTicks(side=1L))
           for (i in seq(npeak)) {
               points(temp,CompSig[,i,drop=TRUE], type="l", 
                      lwd=2, col=lineCol[i])
           } # end for.
           points(temp, rowsumSig, 
                  type="l", lwd=2, col="black")
           legend(ifelse(temp[which.max(signal)]>XaxisCentral,"topleft","topright"),
                  legend=c("Fitted.Curve", paste(seq(npeak),"th-Peak",sep="")), 
                  col=c("black", lineCol[seq(npeak)]), pch=c(21, rep(NA,npeak)),
                  lty=rep("solid",npeak), yjust=2, ncol=1, cex=1.3, bty="o", 
                  lwd=2, pt.bg="white")
           ###
           par(mar=c(5.1,5.1,0,1.1))
           plot(temp, signal-rowsumSig, type="o", 
                xlab="Temperature (K)", ylab="Residuals",
                las=0, lab=c(7,7,9), xaxs="r", yaxs="i", 
                pch=21, bg="black", cex=0.5, cex.lab=1.3)
           abline(h=0)
           box(lwd=2L)
           ###
           par(mar=c(5,4,4,2)+0.1)
           layout(1L)
       } # end if.
       ###
       if (!is.null(outfile)) {
           CompSig <- cbind(temp, signal, rowsumSig,  CompSig)
           colnames(CompSig) <- c("Temperature", "Obs.Signal", 
               "Fit.Signal", paste("Comp.", seq(npeak), sep = ""))
           write.csv(CompSig, file=paste(outfile, ".csv", sep = ""))
       } # end if.
       ###
       return(output)                 
} # end fucntion tgcd.
###

# $Id: boot.fear.R 128 2014-06-14 16:19:18Z B002961 $

# Bootstrap DEA functions, a wrapper for FEAR::boot.sw98


# boot.sw98(XOBS, YOBS, NREP = 2000, DHAT = NULL, 
#    RTS = 1, ORIENTATION = 1, alpha = 0.05, CI.TYPE=2,
# 	  XREF = NULL, YREF = NULL, DREF = NULL, 
# 	  OUTPUT.FARRELL = FALSE, NOPRINT = FALSE, errchk = TRUE)





boot.fear <- function(X,Y, NREP=200, EFF=NULL, RTS="vrs", ORIENTATION="in",
       alpha=0.05, XREF=NULL, YREF=NULL, EREF=NULL)  
{

   if ( !("FEAR" %in% .packages(T)) )
      stop("Note: boot.fear only works if the package FEAR is installed")


   rts_ <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( is.numeric(RTS) )  {
      RTStemp <- rts_[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
   }
   RTS <- tolower(RTS)
   rts <- which(RTS == rts_) -1
   if ( rts < 1 || 3 < rts )
      stop("Invalid value of RTS in call to boot.fear")


   orientation_ <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   orientation <- which(ORIENTATION == orientation_) -1
   if ( orientation < 1 || 3 < orientation )
      stop("Invalid value of ORIENTATION in call to boot.fear")

   if ( !is.null(EREF) && (is.null(XREF) || is.null(YREF)) )
      stop("When EREF is present then XREF and YREF must also be present") 

   if (ORIENTATION=="out") farrell<-TRUE else farrell<-FALSE

   # print(paste("rts =",rts,"   orientation=",orientation))

   b <- NA
   if ( !is.null(XREF) )  XREF <- t(XREF)
   if ( !is.null(YREF) )  YREF <- t(YREF)
   if ( !is.null(EFF)  )  EFF  <- 1/EFF   # Konverter til Shephard
   if ( !is.null(EREF) )  EREF <- 1/EREF  # Konverter til Shephard
##   tryCatch( b <- FEAR::boot.sw98(t(X), t(Y), NREP, EFF, rts, 
##                    orientation, alpha,,XREF, YREF, EREF,
##                    OUTPUT.FARREL=farrell) , 
##      warning = function(w) print(w),
##      error = function(e) {
##         print(e)
##         stop("boot.fear aborted:  Could be that FEAR is not installed")
##      } # ,
##      # finaly=print("FEAR::boot.sw98 finished with bootstrap",quote=FALSE)
##   )

   # boot.sw98(t(x), t(y), NREP=NREP)
   # print("FEAR done; now for calculating the aggregate statistics")

   if ( farrell )  {
      # bb <- b
      bb <- list(eff=b$dhat, eff.bc=b$dhat.bc, bias=b$bias, var=b$var, 
                 conf.int=b$conf.int, eref=b$dref, boot=b$boot, fear=b)
  } else {
      # print("Omregner til Farrell")
      # Omregn til Farrell efficiencer naar det nu ikke er Farrell
      eff <- 1/b$dhat
      bias <- 1/b$dhat - 1/b$dhat.bc
      # bias0 <- rowMeans(1/b$boot) - eff
      # eff.bc <- eff - bias
      # eff.bc1 <- 2*b$eff - 1/rowMeans(b$boot)

      # boot_ <- b$boot - b$dhat
      # ci_ <- 
      #   t(apply(boot_,1, quantile, 
      #     probs=c(0.5*alpha, 1-0.5*alpha), type=9, na.rm=TRUE))
      # ci <-  1/(b$dhat- ci_)
      # rm(boot_, ci_)

      # Forste ordens tilnaermelse til variansen, Taylor, delta
      var <- b$var/b$dhat^2
      var0 <- apply(1/b$boot,1,var)
      if ( is.null(b$dref) )  eref <- NULL  else  eref <- 1/b$dref

      sboot <- t(apply(1/b$boot,1,sort))
 
      bb <- list(eff=eff, 
                 eff.bc=1/b$dhat.bc, # eff.bc0=eff.bc, eff.bc1=eff.bc1,
                 bias=bias, # bias0=bias0,
                 var=var, var0=var0,
                 conf.int=1/b$conf.int[,c(2,1)], 
                 # conf.int0=ci, 
                 eref=eref, boot=sboot, fear=b)

   }

   return( bb )
}


# $Id: addModel.R 125 2013-01-20 16:54:54Z Lars $

# Additive model, corresponds to eqs. 4.34-4.38 in Cooper et al., 2007 
dea.add <- function(X, Y, RTS="vrs", XREF=NULL, YREF=NULL, 
             FRONT.IDX=NULL, param=NULL, TRANSPOSE=FALSE, LP=FALSE)  {

   rts <- c("fdh","vrs","drs","crs","irs","irs","add","fdh+")
   if ( is.numeric(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) ) stop(paste("Unknown scale of returns:", RTS))

	if ( TRANSPOSE ) {
		K <- dim(X)[2]
	} else {
		K <- dim(X)[1]
	}

   if ( RTS == "fdh+" )  {
      # Saet parametrene low og high
      if ( is.null(param) )  {
         param <- .15
      }
      if ( length(param) == 1 )  {
         low <- 1-param
         high <- 1+param
      } else {
         low <- param[1]
         high <- param[2]
      }
      param <- c(low=low, high=high)
   }

	e <- list(eff=rep(1,K), objval=NULL, RTS=RTS, 
        ORIENTATION="in", TRANSPOSE=TRANSPOSE, param=param)
   class(e) <- "Farrell"
	sl <- slack(X,Y,e, XREF=XREF, YREF=YREF, FRONT.IDX=FRONT.IDX, LP=LP)

	# Make slack the sum of slacks, not a logical variable 
	# wether there is slack or not
# 	if (TRANSPOSE)
#      sl$sum2 <- colSums(sl$sx) + colSums(sl$sy)
# 	else 
#      sl$sum2 <- rowSums(sl$sx) + rowSums(sl$sy)

	return(sl)
}  # dea.aff


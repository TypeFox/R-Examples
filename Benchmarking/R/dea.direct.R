# $Id: dea.direct.R 156 2015-07-08 13:34:15Z b002961 $


dea.direct <- function(X,Y, DIRECT, RTS="vrs", ORIENTATION="in", 
                  XREF=NULL, YREF=NULL, FRONT.IDX=NULL, SLACK=FALSE, 
                  param=NULL, TRANSPOSE=FALSE)  {


   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)


   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }
   
   transpose <- FALSE
   if ( TRANSPOSE )  {
      transpose <- TRUE
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }

   m <- dim(X)[2]  # number of inputs
   n <- dim(Y)[2]  # number of outputs
   K <- dim(X)[1]  # number of units, firms, DMUs

   ee <- dea(X,Y, RTS=RTS, ORIENTATION=ORIENTATION, XREF=XREF, YREF=YREF,
             FRONT.IDX=FRONT.IDX, SLACK=SLACK, DUAL=FALSE, 
             DIRECT=DIRECT, param=param, TRANSPOSE=FALSE)  

   mmd <- switch(ORIENTATION, "in"=m, "out"=n, "in-out"=m+n) 
   if ( is.null(ee$objval) )  ee$objval <- ee$eff
   ob <- matrix(ee$objval,nrow=K, ncol=mmd)
   if ( class(DIRECT)=="matrix" && dim(DIRECT)[1] > 1 )  {
       dir <- DIRECT
   } else {
       dir <- matrix(DIRECT,nrow=K, ncol=mmd, byrow=TRUE)
   }
   if ( ORIENTATION=="in" )  {
      e <- 1 - ob*dir/X
   } else if ( ORIENTATION=="out" )  {
      e <- 1 + ob*dir/Y
   } else if ( ORIENTATION=="in-out" )  {
      e <- cbind(1 - dir[,1:m,drop=FALSE]*ob[,1:m,drop=FALSE]/X, 
        1 + dir[,(m+1):(m+n),drop=FALSE]*ob[,(m+1):(m+n),drop=FALSE]/Y)
   } else {
      warning("Illegal ORIENTATION for argument DIRECT") 
   }
   if ( class(e)=="matrix" && dim(e)[2]==1 )
      e <- c(e) 

   if ( transpose )  {
      transpose <- FALSE
      TRANSPOSE <- FALSE
   } 


   if ( TRANSPOSE ) {
      if ( class(e)=="matrix" )
         e <- t(e)
      ee$lambda <- t(ee$lambda)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }

   ee$eff <- e
   ee$direct <- DIRECT
   ee$TRANSPOSE <- TRANSPOSE
   return(ee)
}

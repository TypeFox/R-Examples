# $Id: rev.R 125 2013-01-20 16:54:54Z Lars $

# Function to calculate maximun revenue for given input and given
# output prices.

# Calculations are done with matrices compared to R standards, but
# according to LP practice

# The function cannot be called rev.opt because then print.rev.opt is
# confused with the method rev from base R.

revenue.opt <- function(XREF, YREF, P, XOBS=NULL, RTS="vrs", param=NULL,
            TRANSPOSE=FALSE, LP=FALSE, LPK=NULL)  {
   if ( missing(XOBS) )  {
   	XOBS <- XREF
   }
   if (!TRANSPOSE) {
      XREF <- t(XREF)
      YREF <- t(YREF)
      P <- t(P)
      XOBS <- t(XOBS)
   }
   m = dim(XREF)[1]  # number of inputs
   n = dim(YREF)[1]  # number of outputs
   K = dim(XOBS)[2]  # number of units, firms, DMUs
   Kr = dim(XREF)[2]  # number of units in tevhnology, firms, DMUs

   if ( dim(P)[2] > 1 && dim(P)[2] != K )
      stop("Dimensions for P and XOBS are different")
   if ( Kr != dim(YREF)[2] )
      stop("Number of firms in XREF and YREF differ")
   if ( n != dim(P)[1] )
      stop("Number of inputs in P and YREF differ")
   if ( m != dim(XOBS)[1] )
      stop("Number of outputs in XREF and XOBS differ")

   rts <- c("fdh","vrs","drs","crs","irs","irs","add","fdh+")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.numeric(RTS) )  {
      if (LP) cat(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) cat(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  {
      print(paste("Unknown scale of returns:", RTS))
      print("continuees asssuming RTS = \"vrs\"\n")
      RTS <- "vrs"
   } 

   if ( RTS != "crs" && RTS != "add" )  {
      rlamb <- 2
   } else 
      rlamb <- 0

   lps <- make.lp(m+n +rlamb,n+Kr)
   name.lp(lps, paste("DEA rev,",RTS,"technology"))

   # saet raekker i matrix med restriktioner, saet 0'er for den foerste
   # soejle for den skal alligevel aendres for hver firm.
   dia <- diag(1,nrow=n)
   for ( h in 1:m )
       set.row(lps,h, c(rep(0,n), XREF[h,]))
   for ( h in 1:n)
       set.row(lps,m+h, c(dia[h,], -YREF[h,]))
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(rep(0,n),rep(1,Kr)))
      set.row(lps, m+n+2, c(rep(0,n),rep( -1,Kr)))
   }

   if ( RTS == "fdh" ) {
      set.type(lps,(n+1):(n+Kr),"binary")
      set.rhs(lps,1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
   } else if ( RTS == "vrs" )  {
      set.rhs(lps, c(1,-1), (m+n+1):(m+n+2))
   } else if ( RTS == "drs" )  {
      set.rhs(lps, 1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
   } else if ( RTS == "irs" )  {
      set.rhs(lps, -1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
   } else if ( RTS == "add" )  {
      set.type(lps,(n+1):(n+Kr),"integer")
   } else if ( RTS == "fdh+" )  {
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
      set.rhs(lps, c(high, -low), (m+n+1):(m+n+2))
      add.SOS(lps,"lambda", 1,1, (n+1):(n+Kr), rep(1, Kr))
   }

   set.objfn(lps, c(P[,1],rep(0,Kr)))
   set.constr.type(lps, rep("<=",m+n+rlamb))
   lp.control(lps, sense="max")

   yopt <- matrix(NA,n,K)
   lambda <- matrix(NA,nrow=Kr,ncol=K)
   rev <- rep(NA,K)

   for ( k in 1:K )  {
      if (LP) print(paste("===> firm",k),quote=FALSE)
      if ( dim(P)[2] != 1 && k > 1 ) { 
         set.objfn(lps, c(P[,k],rep(0,K)))
	   }
      set.rhs(lps, XOBS[,k], 1:m)

      if (LP) print(lps)

      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if ( status != 0 ) {
	      print(paste("Error in solving for firm",k,":  Status =",status),
               quote=FALSE)
      }  else {
         rev[k] <- get.objective(lps)
         sol <- get.variables(lps)
         yopt[,k] <- sol[1:n]
         lambda[,k] <- sol[(n+1):(n+Kr)]
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }

   }  # for ( k in 1:K )
   # delete.lp(lps)
   
   rownames(lambda) <- paste("L",1:Kr,sep="")
   names(rev) <- colnames(XOBS)

   if (!TRANSPOSE) {
      yopt <- t(yopt)
      # rev <- t(rev)
      lambda <- t(lambda)
   }

   svar <- list("yopt"=yopt, "rev"=rev, "lambda"=lambda, RTS=RTS, 
        TRANSPOSE=TRANSPOSE)
   class(svar) <- "revenue.opt"
   return (svar)
} # revenue.opt



print.revenue.opt  <- function(x, ...)  {
   a <- cbind("Optimalt output"=x$yopt)
   print(a,...)
   invisible(a)
} ## print.revenue.opt


summary.revenue.opt <- function(object, ...)  {
   cat("Optimal output:\n")
   print.revenue.opt(object)
   cat("Revenue:\n")
   print(object$rev,...)
   cat("Weights (lambda):\n")
   x <- object$lambda
   xx <- format(unclass(x), digits=4)
   if (any(ina <- is.na(x))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(x) < 1e-9) ) 
      xx[i0] <- sub("0.0000", ".", xx[i0])
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(object)
   # printSpMatrix(Matrix(object$lambda),digits=4, col.names=T,...)
   # print(object$lambda,digits=4)
}  ## summary.revenue.opt



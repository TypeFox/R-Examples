# $Id: profit.R 140 2015-05-15 21:48:02Z B002961 $

# Calculates optimal input and output to maximize profit for given
# input and output prices.

profit.opt <- function(XREF, YREF, W, P, RTS="vrs", param=NULL,
                  TRANSPOSE=FALSE, LP=FALSE, LPK=NULL)
{
   if (!TRANSPOSE) {
      XREF <- t(XREF)
      YREF <- t(YREF)
      W <- t(W)
      P <- t(P)
   }
   m = dim(XREF)[1]  # number of inputs
   n = dim(YREF)[1]  # number of outputs
   K = dim(W)[2]  # number of price sets
   Kr = dim(XREF)[2]  # number of units, firms, DMUs

   if ( dim(P)[2] > 1 && dim(P)[2] != K )
      stop("Dimensions for W and P are different")
   if ( Kr != dim(YREF)[2] )
      stop("Number of firms in XREF and YREF differ")
   if ( m != dim(W)[1] )
      stop("Number of inputs in W and XREF differ")
   if ( n != dim(P)[1] )
      stop("Number of outputs in P and YREF differ")

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

   lps <- make.lp(m+n +rlamb,m+n+Kr)
   name.lp(lps, paste("DEA profit,",RTS,"technology"))

   set.objfn(lps, c(-W[,1],P[,1], rep(0,Kr)))
   # saet raekker i matrix med restriktioner, 
   # foerst for input
   dia <- diag(1,nrow=m)
   for ( h in 1:m )
       set.row(lps,h, c(-dia[h,], rep(0,n), XREF[h,]))
   # saa for output
   dia <- diag(1,nrow=n)
   for ( h in 1:n)
       set.row(lps,m+h, c(rep(0,m), dia[h,], -YREF[h,]))
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(rep(0,m+n),rep(-1,Kr)))
      set.row(lps, m+n+2, c(rep(0,m+n),rep( 1,Kr)))
   }

   if ( RTS == "fdh" ) {
      set.type(lps,(m+n+1):(m+n+Kr),"binary")
      set.rhs(lps,1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
   } else if ( RTS == "vrs" )  {
      set.rhs(lps, c(-1,1), (m+n+1):(m+n+2))
   } else if ( RTS == "drs" )  {
      set.rhs(lps, 1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
#   } else if ( RTS == "crs" )  {
#     # En mystisk restriktion for at tvinge loesning til eksisterende firm
#     add.constraint(lps, rep(1,Kr),">=", 1, (m+n+1):(m+n+Kr))
   } else if ( RTS == "irs" )  {
      set.rhs(lps, -1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
   } else if ( RTS == "add" )  {
      set.type(lps,(m+n+1):(m+n+Kr),"integer")
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
      set.rhs(lps, c(-low,high), (m+n+1):(m+n+2))
      add.SOS(lps,"lambda", 1,1, (m+n+1):(m+n+Kr), rep(1, Kr))
   }

   set.constr.type(lps, rep("<=",m+n+rlamb),  1:(m+n+rlamb))
   lp.control(lps, sense="max")

   xopt <- matrix(NA,m,K)
   yopt <- matrix(NA,n,K)
   lambda <- matrix(NA,nrow=Kr,ncol=K)
   profit <- rep(NA,K)

   for ( k in 1:K )  {
      if ( dim(W)[2] != 1 && k > 1 ) { 
         set.objfn(lps, c(-W[,k],P[,k], rep(0,Kr)))
      }

      if (LP) print(lps)

      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if ( status == 3 )  {
         print(paste("Profit is unbounded for firm ",k), quote=FALSE)
         profit[k] <- Inf
         sol <- get.variables(lps)
         xopt[,k] <- sol[1:m]
         yopt[,k] <- sol[(m+1):(m+n)]
         lambda[,k] <- sol[(m+n+1):(m+n+Kr)]
      } else if ( status != 0 ) {
	      print(paste("Error in solving for firm ",k,":  Status =",status),
               quote=FALSE)
      }  else {
         profit[k] <- get.objective(lps)
         sol <- get.variables(lps)
         xopt[,k] <- sol[1:m]
         yopt[,k] <- sol[(m+1):(m+n)]
         lambda[,k] <- sol[(m+n+1):(m+n+Kr)]
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }

   }  # for ( k in 1:K )
   # delete.lp(lps)

   rownames(lambda) <- paste("L",1:Kr,sep="")
   names(profit) <- colnames(W)

   if (!TRANSPOSE) {
      xopt <- t(xopt)
      yopt <- t(yopt)
      # profit <- t(profit)
      lambda <- t(lambda)
   }

   svar <- list("xopt"=xopt, yopt=yopt, profit=profit, "lambda"=lambda, 
        RTS=RTS, TRANSPOSE=TRANSPOSE)
   class(svar) <- "profit.opt"
   return (svar)

}  # profit.opt



print.profit.opt  <- function(x, ...)  {
   a <- cbind("Optimalt input"=x$xopt, "Optimal output"=x$yopt)
   print(a,...)
   invisible(a)
} ## print.profit.opt


summary.profit.opt <- function(object, ...)  {
   cat("Optimal input and output:\n")
   print.profit.opt(object)
   cat("Profit:\n")
   print(object$profit,...)
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
}  ## summary.profit.opt



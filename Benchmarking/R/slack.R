# $Id: slack.R 153 2015-07-07 14:01:23Z b002961 $

# Calculate slack at the efficient points.

# Hvis data ikke er transponeret bliver x og y transponeret og alle
# beregninger laves transponeres. Det betyder at resultat skal
# transponeres inden afslut og retur.

slack <- function(X, Y, e, XREF=NULL, YREF=NULL, FRONT.IDX=NULL, 
                  LP=FALSE)  {       # , CONTROL=NULL, LPK=NULL
   if ( !methods::is(e,"Farrell") )  {
       stop("In call of slack: argument 'e' must be of class 'Farrell'")
   }

   RTS <- e$RTS
   if ( RTS == "fdh0" ) RTS <- "fdh"
   if (LP) print(paste("slack:  RTS =",RTS),quote=F)
   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }
   if ( e$TRANSPOSE ) {
      if (LP) print("X and Y are transposed")
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
   }

   okr <- dim(XREF)[1]
   if (LP) cat("okr =",okr,"\n")

   if ( FALSE && length(FRONT.IDX) == 0 )  {
      # Kun units med eff==1 bruges i referenceteknologi
      FRONT.IDX <- which( abs(e$eff -1) < 1e-5 )
   }
   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in method 'slack'")
      XREF <- XREF[FRONT.IDX,, drop=FALSE]
      YREF <- YREF[FRONT.IDX,, drop=FALSE]
      if (LP) cat("FRONT.IDX =",FRONT.IDX,"\n")
   }

   K = dim(X)[1]  # number of units, firms, DMUs
   Kr = dim(XREF)[1]  # number of units in reference set
   m = dim(X)[2]  # number of inputs
   n = dim(Y)[2]  # number of outputs
   if (LP) cat("In slack: m n K Kr = ",m,n,K,Kr,"\n")

   if ( m != dim(XREF)[2] )
      stop("Number of inputs must be the same in X and XREF")
   if ( n != dim(YREF)[2] )
      stop("Number of outputs must be the same in Y and YREF")
   if ( K !=  dim(Y)[1] )
      stop("Number of units must be the same in X and Y")
   if ( Kr != dim(YREF)[1] )
      stop("Number of units must be the same in XREF and YREF")

   # For at undgaa afrundingsproblemer i forbindelse med slacks skal
   # in- og ouput skaleres til at vaere i omegnen af 1
   mmm <- (colMeans(X))
   nnn <- (colMeans(Y))
   if ( min(mmm) < 1e-4 || max(mmm) > 1e4 || 
       min(nnn) < 1e-4 || max(nnn) > 1e4 )  {
      SKALERING <- TRUE
      X <- X / matrix(mmm, nrow=K, ncol=m, byrow=TRUE)
      XREF <- XREF / matrix(mmm, nrow=Kr, ncol=m, byrow=TRUE)
      Y <- Y / matrix(nnn, nrow=K, ncol=n, byrow=TRUE)
      YREF <- YREF /  matrix(nnn, nrow=Kr, ncol=n, byrow=TRUE)
   } else {
     SKALERING <- FALSE
   }

   if ( RTS != "crs" && RTS != "add" )  {
      rlamb <- 2
   } else { 
      rlamb <- 0
   }

   if ( is.null(e$objval) )
      eff <- e$eff
   else
      eff <- e$objval

   if ( is.null(e$direct) )  {
      if ( e$ORIENTATION == "in" ) {
         E <- eff
         FF <- rep(1,K)  # Naturligt at bruge F, men F er FALSE i R
      } else if ( e$ORIENTATION == "out" )  {
         FF <- eff
         E <- rep(1,K)
      } else if ( e$ORIENTATION == "graph" )  {
         E <- eff
         FF <- 1/eff
      } else {
        stop(paste("Unknown orientation in slack:", e$ORIENTATION))
      }
   } else {
      mmd <- switch(e$ORIENTATION, "in"=m, "out"=n, "in-out"=m+n) 
      ob <- matrix(e$objval,nrow=K, ncol=mmd)
      if ( class(e$direct)=="matrix" && dim(e$direct)[1] > 1 )  {
          dir <- e$direct
      } else {
          dir <- matrix(e$direct,nrow=K, ncol=mmd, byrow=TRUE)
      }

      if ( e$ORIENTATION=="in" )  {
         dirRhs <- cbind(X - ob*dir, -Y)
      } else if ( e$ORIENTATION=="out" )  {
         dirRhs <- cbind(X, -Y - ob*dir)
      } else if ( e$ORIENTATION=="in-out" )  {
         dirRhs <- cbind(X -ob[,1:m,drop=FALSE]*dir[,1:m,drop=FALSE], 
           -Y -ob[,(m+1):(m+n),drop=FALSE]*dir[,(m+1):(m+n),drop=FALSE])
      } else {
         warning("Illegal ORIENTATION for argument DIRECT in slacks",
                 immediate. = TRUE) 
      }

   }  # if ( is.null(e$direct) )


   # Initialiser LP objekt
   lps <- make.lp(m+n +rlamb, m+n+Kr)
   name.lp(lps, 
          paste("DEA-slack",RTS,",",e$ORIENTATION,"orientated",sep="-"))
   
   # saet raekker i matrix med restriktioner
   for ( h in 1:m )
       set.row(lps,h, XREF[,h], (m+n+1):(m+n+Kr) )
   for ( h in 1:n)
       set.row(lps,m+h, -YREF[,h], (m+n+1):(m+n+Kr) )
   for ( h in 1:(m+n) )
       set.mat(lps,h,h,1)
   # restriktioner paa lambda
   if ( RTS != "crs" && RTS != "add" )  {
      set.row(lps, m+n+1, c(rep(0,m+n),rep(-1,Kr)))
      set.row(lps, m+n+2, c(rep(0,m+n),rep( 1,Kr)))
   }
   set.constr.type(lps, c(rep("=",m+n), rep("<=",rlamb)))
   set.objfn(lps, rep(1, m+n), 1:(m+n))
   lp.control(lps, sense="max")

   if ( RTS == "fdh" ) {
      set.type(lps,(m+n+1):(m+n+Kr),"binary")
      delete.constraint(lps, m+n+1)
      set.rhs(lps,1, m+n+1)
      rlamb <- rlamb -1
   } else if ( RTS == "vrs" )  {
      set.rhs(lps, c(-1,1), (m+n+1):(m+n+2))
   } else if ( RTS == "drs" )  {
      set.rhs(lps, -1, m+n+1)
      delete.constraint(lps, m+n+2)
      rlamb <- rlamb -1
      set.constr.type(lps, ">=", m+n+1)
   } else if ( RTS == "irs" )  {
      set.rhs(lps, 1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
      set.constr.type(lps, ">=", m+n+1)
   } else if ( RTS == "irs2" )  {
      set.rhs(lps, 1, m+n+2)
      delete.constraint(lps, m+n+1)
      rlamb <- rlamb -1
      set.constr.type(lps, ">=", m+n+1)
      set.semicont(lps, 2:(1+Kr))
      set.bounds(lps, lower=rep(1,Kr), columns=(m+n+1):(m+n+Kr))
   } else if ( RTS == "add" )  {
      set.type(lps, (m+n+1):(m+n+Kr),"integer")
   } else if ( RTS == "fdh+" )  {
      param <- e$param
      low <- param["low"]
      high <- param["high"]
      set.rhs(lps, c(-low,high), (m+n+1):(m+n+2))
      add.SOS(lps,"lambda", 1,1, (m+n+1):(m+n+Kr), rep(1, Kr))
   }
 
   if (LP) {
      print("Right hand side not sat yet")
      print(lps)
   }

   objval <- rep(NA,K)   # vector with sums of slacks
   sx <- matrix(NA,K,m)
   sy <- matrix(NA,K,n)
   lambda <- matrix(NA,K,Kr)

   cont <- lp.control(lps)
   eps <- cont$epsilon["epsint"]
   # if ( !is.null(CONTROL) )  {
   #    lp.control(lps,CONTROL)
   # }
   for ( k in 1:K )  { # beregn for hver enhed
      if (LP) cat("\n---\nFirm",k,"\n")

      # her er der problem med afrunding, sammen med venstrsiden er
      # det ikke sikkert at restriktionen holder naar der er
      # sket afrunding i mellemregningerner.
      if ( is.null(e$direct) )  {
         rhs <- c(E[k] * X[k,], -FF[k] * Y[k,])
      } else {
         rhs <- dirRhs[k,]
      }
      set.rhs(lps, rhs, 1:(m+n))
      if (LP)  {
         print(paste("Hoejresiden for firm",k))
         if ( is.null(e$direct) )  {
            print(E[k] * X[k,])
            print( -FF[k] * Y[k,])
         } else {
            print(rhs)
         }
         print(lps)
      }
      if ( !is.null(LP) && k %in% LP )  {
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }

      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if ( status != 0 ) {
         if ( status == 2 || status == 3 ) {
	        # print(paste("Firm",k,"not in the technology set"), quote=F)
           # print(paste("Status =",status))
           objval[k] <- 0
           sol <- rep(0,m+n+Kr)
           sol[m+n+k] <- 1
         } else {
	        print(paste("Error in solving for firm",k,":  Status =",status), 
              quote=F)
           objval[k] <- NA
           sol <- NA
         }
      }  else {
         objval[k] <- get.objective(lps)
         sol <- get.variables(lps)
      }
      sx[k,] <- sol[1:m]       
      sy[k,] <- sol[(m+1):(m+n)]
      lambda[k,] <- sol[(m+n+1):(m+n+Kr)]

      if (LP)  {
         print(paste("Obj.value =",get.objective(lps)))
         cat(paste("Solutions for",k,": "))
         print(sol)
      }

   }  # loop for each firm

   # Drop soejler i lambda der alene er nuller, dvs lambda soejler skal 
   # alene vaere reference firms

   colnames(sx) <- paste("sx",1:m,sep="")
   colnames(sy) <- paste("sy",1:n,sep="")

   # For at undgaa at regneunoejagtigheder giver ikke-nul slack bliver
   # slack sat til nul hvis de er mindre end den relative noejagitghed
   # i input og output, noejagtighed der er i beregningerne.  Det
   # betyder at sum af slacks bliver sum mens objval er sum plus
   # regneunoejagtighed som giver en forskel hvis X er meget stoerre
   # eller mindre end 1.

   sx[abs(sx) < eps ] <- 0
   sy[abs(sy) < eps ] <- 0
   if ( SKALERING )  {
      sx <- sx * matrix(mmm, nrow=K, ncol=m, byrow=TRUE)
      sy <- sy * matrix(nnn, nrow=K, ncol=n, byrow=TRUE)
   } 
   sum <- rowSums(sx) + rowSums(sy)

   if ( length(FRONT.IDX)>0 )  {
      colnames(lambda) <- paste("L",(1:okr)[FRONT.IDX],sep="")
   } else {
      colnames(lambda) <- paste("L",1:Kr,sep="")
   }


   if (FALSE && LP)  {
   	print("Done with loop\nColumn names for lambda:")
	   print(colnames(lambda))
	   print(lambda)
   }

   if ( e$TRANSPOSE ) {
      sx <- t(sx)
      sy <- t(sy)
      lambda <- t(lambda)
   }

   if ( FALSE && LP ) {
      print("Faerdig med slack: lambda")
	   print(colnames(lambda))
      print(lambda)
      print("slack:")
      print(paste("laengden af objval:",length(objval)))
      print(objval==0)
      # print(!slackZero)
   }

   oe <- list(eff=eff, slack=sum>eps, sum=sum, objval=objval, sx=sx, sy=sy,
        lambda=lambda, RTS=e$RTS, ORIENTATION=e$ORIENTATION, 
        TRANSPOSE=e$TRANSPOSE)

   class(oe) <- "slack"

	return(oe)
} # slack






print.slack  <- function(x, digits=4, ...)  {
   sum <- x$sum
   print(sum, digits=digits, ...)
   invisible(sum)
} ## print.slack



summary.slack <- function(object, digits=4, ...)  {
   eps <- 1e-6
   cat("Efficiency and slacks:\n")

   # if ( object$ORIENTATION != "out" ) 
   #    a <- cbind(E=round(object$eff, digits=2),
   #             "sum"=object$sum, object$sx, object$sy)
   # else
   #    a <- cbind(F=round(object$eff,2),
   #             "sum"=object$sum, object$sx, object$sy)
   # print(a,...)

   cat("Number of firms with efficiency==1 and positive slacks:",
         sum(abs(object$eff-1) < eps & object$slack ),"\n" )

   if ( dim(object$sx)[2]==1 )  {
     SX <- object$sx > eps
   } else {
     SX <- rowSums(object$sx) > eps  
   }
   if ( dim(object$sy)[2]==1 )  {
     SY <- object$sy > eps
   } else {
     SY <- rowSums(object$sy) > eps
   }
   cat("Number of firms with:\n")
   cat("  only x slacks: ", sum(SX & !SY), "\n")
   cat("  only y slacks: ", sum(!SX & SY), "\n")
   cat("  x and y slacks:", sum(SY & SX), "\n")

   cat("\n  x slacks: ", sum(SX), "\n")
   cat(   "  y slacks: ", sum(SY), "\n")
   cat(   "all slacks: ", sum(SX | SY), "\n")


   # invisible(a)
   if ( FALSE && dim(lambda(object))[2] <= 10 )  {
      cat("Weights (lambda):\n")
      x <- lambda(object)
      xx <- round(x, digits)
      # xx <- format(unclass(x), digits=digits)
      if (any(ina <- is.na(x))) 
         xx[ina] <- ""
      if ( any(i0 <- !ina & abs(x) < 1e-9) ) 
         xx[i0] <- sub("0", ".", xx[i0])
      print(xx, quote=FALSE, rigth=TRUE, ...)
      invisible(xx)
   }
}


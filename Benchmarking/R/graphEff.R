# $Id: graphEff.R 140 2015-05-15 21:48:02Z B002961 $

# Funktion til beregning af graf efficiens.  Beregning sker via
# bisection hvor der itereres mellem mulige og ikke-mulige loesninger
# et LP problem hvor venstreside er som in- og output orienteret
# efficiens. blot er foerste soejle erstattet af rene nuller, og G*X og
# (1/G)*Y optraeder paa hoejresiden. Minimering af 0 saa der blot soeges om
# der er en mulig loesning. Da soejlen for efficiens er bar 0'er vil
# justering af efficiens ud over G ikke ske, dvs. det er kun lambdaer
# der tilpasses for at se om der er en mulig loesning.

graphEff <- function(lps, X, Y, XREF, YREF, RTS, FRONT.IDX, rlamb, oKr, 
         param=param, TRANSPOSE=FALSE, SLACK=FALSE, FAST=FALSE, LP=FALSE) 
{
   m = dim(X)[2]  # number of inputs
   n = dim(Y)[2]  # number of outputs
   K = dim(X)[1]  # number of units, firms, DMUs
   Kr = dim(YREF)[1]  # number of units, firms, DMUs

   objval <- rep(NA,K)   # vector for the final efficiencies
   if ( FAST ) {
     lambda <- NULL
   } else {
      lambda <- matrix(NA, nrow=K, ncol=Kr) # lambdas one column per unit
   }
   set.column(lps, 1, rep(0,dim(lps)[1]))
   lpcontr <- lp.control(lps)
   tol <- lpcontr$epsilon["epsint"]
   for ( k in 1:K)  {
      if ( LP )  print(paste("Firm",k), quote=FALSE)
      # Lav bisection
      a <- 0
      b <- 2  # medfoerer start med G=1
      nIter <- 0
      gFundet <- FALSE
      xIset <- TRUE

      # Er G==1?
      G <- 1
      set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if (LP) print(paste("For G=1, status =",status))
      nIter <- nIter + 1
      if ( status == 0 )  {
         G <- 1 - sqrt(tol)
         set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
         if (RTS=="add") set.basis(lps, default=TRUE)
         status <- solve(lps)
         if (LP) print(paste("For G=1-eps, status =",status))
         nIter <- nIter + 1
         if ( status != 0 )  {
            # G=1 er mulig og G=1-eps er ikke-mulig ==> G==1
            G <- b <- 1
            gFundet <- TRUE
         }
      } else {
         # warning("Firm outside technology set, firm =",k)
         # G=1 er ikke mulig; firm uden for teknology set saa G > 1
         xIset <- FALSE
         # Bestem oevre graense
         b <- 2
         while (status != 0 && nIter < 50)  {
            set.rhs(lps, c(-b*X[k,], Y[k,]/b), 1:(m+n))
            status <- solve(lps)
            if ( status==5 )  {
               set.basis(lps, default=TRUE)
               status <- solve(lps)
            }
            nIter <- nIter + 1
            b <- b^2
         }
         # nedre graense
         if ( b > 2 )  { a <- sqrt(b) 
	      }  else  { a <- 1 }
      }

      status <- 0 # status kunne godt have en anden vaerdi og saa 
                  # ville naeste loekke ikke blive gennemloebet
      if ( !gFundet && xIset && status == 0 )  {
         # G==1 er mulig er G er ikke fundet endnu
         dif <- .1
         i <- 1
         while ( status==0 && i < 10 )  {
            G <- 1 - i*dif
            set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
            if (RTS=="add") set.basis(lps, default=TRUE)
            status <- solve(lps)
            if ( status==5 )  {
               set.basis(lps, default=TRUE)
               status <- solve(lps)
            }
            if (LP) print(paste("G = ",G,"; status = ",status))
            nIter <- nIter + 1
            i <- i+1
         }
         # enten er i==10 eller ogsaa er status!=0
         if ( i==10 )  {
            a <- 0
            b <- dif
         } else {
            a <- 1 - (i-1)*dif
            b <- 1 - (i-2)*dif
         }
      }


      # bisection loekke
      if (LP) print(paste("Bisection interval: [",a,",",b,"]"))
      while ( !gFundet && b-a > tol && nIter < 50 )  {
         G <- (a+b)/2
         set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
         if (RTS=="add") set.basis(lps, default=TRUE)
         # if ( k==1 ) print(lps)
         status <- solve(lps)
         if (LP) print(paste("G = ",G,"(",k,"); status =",status))
         if ( status == 0 ) {
            # loesning findes
            b <- G
         } else {
            a <- G
         }
         nIter <- nIter + 1
      }
      if (LP) print(paste("nIter =",nIter,"; status =",status))

      if ( status != 0 )  {
         # Hvis den sidste vaerdi af G ikke var mulig bruger vi den
         # oevre graense. Det er noedvendigt med en mulig loesning for at
         # kunne faa lambdaer og duale vaerdier.
         G <- b
	      set.rhs(lps, c(-G*X[k,],Y[k,]/G), 1:(m+n))
         status <- solve(lps)
	   }
      if (LP)  {
         print(paste("G = ",G,"(",k,"); status =",status))
         print(rlamb)
         print("Solution")
         print(get.variables(lps))
         print(lps)
      }
      objval[k] <- G
      if ( LP && k == 1 )  print(lps)
      if ( !FAST )  {
         sol <- get.variables(lps)
         lambda[k,] <- sol[2:(1+Kr)]
      }
   	if (LP && status==0) {
         print(paste("Objval, firm",k))
         print(get.objective(lps))
         print("Solution/varaibles")
         print(get.variables(lps))
         print("Primal solution")
         print(get.primal.solution(lps))
         print("Dual solution:")
         print(get.dual.solution(lps))
      }
   }  # loop for each firm

   e <- objval
   e[abs(e-1) < sqrt(tol)] <- 1

   if ( FAST ) { 
      return(e)
      stop("Her skulle vi ikke kunne komme i 'dea'")
   }

   if ( length(FRONT.IDX)>0 )  {
      colnames(lambda) <- paste("L",(1:oKr)[FRONT.IDX],sep="")
   } else {
      colnames(lambda) <- paste("L",1:Kr,sep="")
   }

   primal <- dual <- NULL
   ux <- vy <- NULL

   if ( TRANSPOSE ) {
      lambda <- t(lambda)
   }

   oe <- list(eff=e, lambda=lambda, objval=objval, RTS=RTS,
              primal=primal, dual=dual, ux=ux, vy=vy, gamma=gamma,
              ORIENTATION="graph", TRANSPOSE=TRANSPOSE
              # ,slack=slack_, sx=sx, sy=sy
              , param=param 
              )
   class(oe) <- "Farrell"

   if ( SLACK ) {
      if ( TRANSPOSE )  { # Transponer tilbage hvis de blev transponeret
         X <- t(X)
         Y <- t(Y)
         XREF <- t(XREF)
         YREF <- t(YREF)
      }
      sl <- slack(X, Y, oe, XREF, YREF, FRONT.IDX, LP=LP)
      oe$slack <- sl$slack
      oe$sx <- sl$sx
      oe$sy <- sl$sy
      oe$lambda <- sl$lambda
      if (LP)  {
         print("slack fra slack:")
         print(sl$slack)
         print("slack efter slack:")
         print(oe$slack)
      }
   }

   return(oe)
}


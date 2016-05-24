# $Id: dea.dual.R 156 2015-07-08 13:34:15Z b002961 $

# In the calculation in the method input/output matrices X and Y are
# of the order good x firms.  Ie. X, Y etc must be transformed as
# default in R is firm x good.

# DUALIN og DUALOUT eller blot DUAL.  Her er valgt DUAL saa DUAL skal
# vaere en matrix saaledes at der en raekke i DUAL for hvert input og
# hvert output, bortset fra det foerste input og det foerste output,
# dvs m-1+n-1=m+n-2 raekker.  Der skal vaere 2 soejler, foerste soejle
# er den nedre graense, og anden soejle er den oevre graense.


dea.dual <- function(X,Y, RTS="vrs", ORIENTATION="in", 
            XREF=NULL,YREF=NULL,
            FRONT.IDX=NULL, DUAL=NULL, DIRECT=NULL,
            TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  {
   # XREF, YREF determines the technology
   # FRONT.IDX index for units that determine the technology

   rts <- c("fdh","vrs","drs","crs","irs","irs","add","fdh+")
   if ( missing(RTS) ) RTS <- "vrs" 
   if (LP)  print(paste("Vaerdi af 'RTS' er ",RTS),quote=FALSE)
   if ( is.numeric(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=FALSE)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=FALSE)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  stop(paste("Unknown scale of returns:", RTS))

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop(paste("Unknown value for ORIENTATION:",ORIENTATION),quote=F)
   }

   if ( RTS %in% c("fdh","add", "fdh+") )
      stop("dea.dual does not work for \"fdh\", \"fdh+\" or \"add\"")
   if ( !ORIENTATION %in% c("in","out","in-out") )
      stop("dea.dual does not work for \"graph\"")

   if ( class(X)=="data.frame" && data.kontrol(X) || is.numeric(X) ) 
      { X <- as.matrix(X) }
   if ( class(Y)=="data.frame" && data.kontrol(Y) || is.numeric(Y) ) 
      { Y <- as.matrix(Y) }
   if ( class(XREF)=="data.frame" && data.kontrol(XREF)||is.numeric(XREF))
      { XREF <- as.matrix(XREF) }
   if ( class(YREF)=="data.frame" && data.kontrol(YREF)||is.numeric(YREF)) 
      { YREF <- as.matrix(YREF) }

   if ( class(X)!="matrix" || !is.numeric(X) )
      stop("X is not a numeric matrix (or data.frame)")
   if ( class(Y)!="matrix" || !is.numeric(X) )
      stop("Y is not a numeric matrix (or data.frame)")
   if ( !is.null(XREF) && (class(XREF)!="matrix" || !is.numeric(XREF)) )
      stop("XREF is not a numeric matrix (or data.frame)")
   if ( !is.null(YREF) && (class(YREF)!="matrix" || !is.numeric(YREF)) )
      stop("YREF is not a numeric matrix (or data.frame)")


   if ( missing(XREF) || is.null(XREF) )  {
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      YREF <- Y
   }
   
   if ( TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & class(DIRECT)=="matrix" )
         DIRECT <- t(DIRECT)
   }
   orgKr <- dim(XREF)

   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX) )
         stop("FRONT.IDX is not a vector in 'eff'")
      XREF <- XREF[FRONT.IDX,, drop=FALSE]
      YREF <- YREF[FRONT.IDX,, drop=FALSE]
   }

   m = dim(X)[2]  # number of inputs
   n = dim(Y)[2]  # number of outputs
   K = dim(X)[1]  # number of units, firms, DMUs
   Kr = dim(XREF)[1]  # number of units, firms, DMUs
   if ( !is.null(DIRECT) )  {
      if ( class(DIRECT)=="matrix" ) {
         md <- dim(DIRECT)[2]
         Kd <- dim(DIRECT)[1]
      } else {
         md <- length(DIRECT)
         Kd <- 0
      }
   } else {
      Kd <- 0
   }

   if (LP) cat("m n k kr = ",m,n,K,Kr,"\n")
   if (LP & !is.null(DIRECT) ) cat("md, Kd =",md,Kd,"\n") 

   if ( m != dim(XREF)[2] )
      stop("Number of inputs must be the same in X and XREF")
   if ( n != dim(YREF)[2] )
      stop("Number of outputs must be the same in Y and YREF")
   if ( K != dim(Y)[1] )
      stop("Number of units must be the same in X and Y")
   if ( Kr != dim(YREF)[1] )  
      stop("Number of units must be the same in XREF and YREF")

   if ( !is.null(DUAL) && ( !is.matrix(DUAL) ||
               dim(DUAL)[1] != (m+n-2) || dim(DUAL)[2] != 2 ) ) {
      stop("DUAL must be a (m+n-2) x 2 matrix with lower and upper bounds for restrictions")
   }


   if ( !is.null(DIRECT) & ORIENTATION=="graph" )
         stop("DIRECT cannot not be used with ORIENTATION=\"graph\"")
 
   if ( !is.null(DIRECT) & length(DIRECT) > 1 )  {
      if ( ORIENTATION=="in" & md!=m )
         stop("Length of DIRECT must be the number of inputs")
      else if ( ORIENTATION=="out" & md!=n )
         stop("Length of DIRECT must be the number of outputs")
      else if ( ORIENTATION=="in-out" & md!=m+n )
         stop("Length of DIRECT must be the number of inputs plus outputs")
      if ( class(DIRECT)=="matrix" & (Kd>0 & Kd!=K) )
         stop("Number of firms in DIRECT must equal firms in X and Y") 
   }
   if ( !is.null(DIRECT) & length(DIRECT) == 1 )  {
      if ( ORIENTATION=="in" & length(DIRECT)!=m )
         DIRECT <- rep(DIRECT,m)
      else if ( ORIENTATION=="out" & length(DIRECT)!=n )
         DIRECT <- rep(DIRECT,n)
      else if ( ORIENTATION=="in-out" & length(DIRECT)!=m+n )
         DIRECT <- rep(DIRECT,m+n)
   }


   if ( RTS == "vrs" )  { 
      rlamb <- 2
   } else if ( RTS == "drs" || RTS == "irs" )  {
      rlamb <- 1
   } else if ( RTS == "crs" )  {
      rlamb <- 0
   } else {
      stop(paste("Unknown value for RTS in 'dea.dual':", RTS))
   }



if ( !missing(DUAL) && !is.null(DUAL) )  {
   # Make matrix for dual restrictions.
   # Foerst for input hvis der er mere end 1 input
   if ( m > 1 )  {
      Udiag <- diag(1,m-1)
      DL <- rbind(DUAL[1:(m-1),1], -Udiag)
      DU <- rbind(-DUAL[1:(m-1),2], Udiag)
      if (LP)  {
         print("DL"); print(DL)
         print("DU"); print(DU)
      }
      ADin <- cbind(DL,DU)
   } else ADin <- NULL
   if (LP)  {
      print("ADin"); print(ADin)
   }

   # og saa restriktioner for output hvis der er mere end 1 output
   if ( n > 1 )  {
      Udiag <- diag(1,n-1)
      DL <- rbind(DUAL[m:(m+n-2),1], -Udiag)
      DU <- rbind(-DUAL[m:(m+n-2),2], Udiag)
      ADout <- cbind(DL,DU)
   } else ADout <- NULL

   nulln <- matrix(0,nrow=n,ncol=2*(m-1))
   nullm <- matrix(0,nrow=m,ncol=2*(n-1))
   AD <- cbind(rbind(ADin,nulln),rbind(nullm,ADout))
   # AD <- t(AD)
   if (LP) {
      print("AD"); print(AD)
   }
} else {
  # AD <- matrix(0, nrow=m+n, ncol=2*(m+n-2))
  AD <- NULL
}

if ( is.null(AD) )  {
   restr <- 0
} else {
   restr <- 2*(m-1 +n-1)
}


   # Initialiser LP objekt
   lps <- make.lp(1+Kr+restr, m+n+rlamb )
   name.lp(lps, paste(ifelse(is.null(AD),"Dual","DualAC"),ORIENTATION,
             RTS,sep="-"))
   # if ( LP==TRUE ) print(lps)

   # saet soejler i matrix med restriktioner, saet 0'er for den foerste
   # raekke for den skal alligevel aendres for hver firm.
   for ( h in 1:m ) 
      set.column(lps,h, c( 0, -XREF[,h], AD[h,]))
   # if ( LP || !is.null(LPK) ) print(lps)
   for ( h in 1:n)
       set.column(lps,m+h, c( 0, YREF[,h], AD[m+h,]))
   # restriktioner paa lambda
   objgamma <- NULL
   if ( rlamb > 0 )  {
      set.column(lps, m+n+1, c(0,rep(-1,Kr)), 1:(Kr+1))
      objgamma <- -1
   }
   if ( rlamb > 1 )  {
      set.column(lps, m+n+2, c(0,rep( 1,Kr)), 1:(Kr+1))
      objgamma <- c(-1,1)
   }
   if ( rlamb == 1 && RTS == "irs" )  {
      set.column(lps, m+n+1, c(0,rep(1,Kr)), 1:(Kr+1))
      objgamma <- 1
   }

   if ( !is.null(DIRECT) & Kd==0 )  {
      # print(Kd)
      # print(DIRECT)
      # Samme retning for alle enheder
      if ( ORIENTATION=="in" )
         set.row(lps, 1, c(-DIRECT),1:m)
      else if ( ORIENTATION=="out" )
         set.row(lps, 1, c(-DIRECT),(m+1):(m+n))
      else if ( ORIENTATION=="in-out" )
         set.row(lps, 1, c(-DIRECT),1:(m+n))
   }

   set.constr.type(lps, rep("<=", 1+Kr+restr))
 
   if ( ORIENTATION == "in" )  {
      lp.control(lps, sense="max")
      set.rhs(lps, 1, 1)
   } else if ( ORIENTATION == "out" )  {
      lp.control(lps, sense="min")
      set.rhs(lps, -1, 1)
      if ( !is.null(objgamma) ) objgamma <- -objgamma
   } 

   if ( !is.null(DIRECT) )  {
      lp.control(lps, sense="min")
      set.rhs(lps, c(-1,rep(0,Kr)))
      if ( ORIENTATION == "in" | ORIENTATION == "in-out" )      
         if ( !is.null(objgamma) ) objgamma <- -objgamma
   }

   if ( !is.null(CONTROL) )  {
      if( !is.list(CONTROL)) {
         stop( "argument 'control' must be a 'list' object")
      }
      do.call( lp.control, c( list( lprec = lps ), CONTROL ) )
   }


# if ( LP || !is.null(LPK) ) print(lps)

   u <- matrix(NA,K,m)   # vector for the final efficiencies
   v <- matrix(NA,K,n)   # vector for the final efficiencies
   objval <- rep(NA,K)   # vector for the final efficiencies
   sol <- matrix(NA, K, 1 + sum(dim(lps)) )
   gamma <- NULL
   if (rlamb > 0)  {
      gamma <- matrix(NA,K,rlamb)
   }



   for ( k in 1:K )  { # Finds the efficiencies for each unit
       if ( is.null(DIRECT) )  {
          # object function and first collumn
          if ( ORIENTATION == "in" )  {
             objrow <- c(rep(0,m),Y[k,], objgamma)
             if (LP) { 
                print(lps); 
                print("objgamma:") 
                print(objgamma) 
                print("objrow:")
                print(objrow)
             }
             set.objfn(lps, objrow)
             set.row(lps, 1, c(X[k,]),1:m)
          }  else  {
             objrow <- c(X[k,], rep(0,n), objgamma)
             set.objfn(lps, objrow)
             set.row(lps, 1, c(-Y[k,]),(m+1):(m+n))
          }
       } else {
          objrow <- c(X[k,],-Y[k,], objgamma)
          set.objfn(lps, objrow)
          # print(Kd)
          # print(DIRECT)
          if ( Kd > 1 )  {
             # retning for enheden
             if ( ORIENTATION=="in" )
                set.row(lps, 1, -DIRECT[k,], 1:m)
             else if ( ORIENTATION=="out" )
                set.row(lps, 1, -DIRECT[k,], (m+1):(m+n))
             else if ( ORIENTATION=="in-out" )
                set.row(lps, 1, -DIRECT[k,], 1:(m+n))
          }
       }


      if ( LP )  print(paste("Firm",k), quote=FALSE)
      if ( LP && k == 1 )  print(lps)
      set.basis(lps, default=TRUE)
      status <- solve(lps)

      if ( status != 0 ) {
        if (status == 2 || status == 3) {
	        # print(paste("Firm",k,"not in the technology set"), quote=F)
           objval[k] <- ifelse(ORIENTATION=="in",Inf,-Inf)
        } else {
	        print(paste("Error in solving for firm",k,":  Status =",status), 
             quote=F)
           objval[k] <- NA
        }
      }  else {
         objval[k] <- get.objective(lps)

         losning <- get.variables(lps)
         # sol[k,] <- get.variables(lps)
         sol[k,] <- get.primal.solution(lps)
         u[k,] <- losning[1:m]
         v[k,] <- losning[(m+1):(m+n)]
         if ( rlamb > 0 )  {
            gamma[k,] <- losning[(m+n+1):(m+n+rlamb)]
         }
      }


   	if (LP && status==0) {
         print(paste("Objval, firm",k))
         print(objval[k])
         print("Solution")
         print(sol)
         print("Dual values:")
         print(u[k,])
         print(v[k,])
         print("get.variables")
         print(get.variables(lps))
         print("Primal solution")
         print(get.primal.solution(lps))
         print("Dual solution:")
         print(get.dual.solution(lps))
      }

      if ( !is.null(LPK) && k %in% LPK )  {
         print(paste("Model",k,"(",name.lp(lps),")"))
         print(lps)
         write.lp(lps, paste(name.lp(lps),k,".mps",sep=""),
                type="mps",use.names=TRUE)
      }
   } #    for ( k in 1:K )

   if ( is.null(DIRECT) )  {
      eff <- objval
   } else { 
      mmd <- switch(ORIENTATION, "in"=m, "out"=n, "in-out"=m+n) 
      ob <- matrix(objval, nrow=K, ncol=mmd)
      if ( class(DIRECT)=="matrix" && dim(DIRECT)[1] > 1 )  {
          dir <- DIRECT
      } else {
          dir <- matrix(DIRECT, nrow=K, ncol=mmd, byrow=TRUE)
      }
      if ( ORIENTATION=="in" )  {
         eff <- 1 - ob*dir/X
      } else if ( ORIENTATION=="out" )  {
         eff <- 1 + ob*dir/Y
      } else if ( ORIENTATION=="in-out" )  {
         eff <- cbind(1 - ob[,1:m,drop=FALSE]*dir[,1:m,drop=FALSE]/X, 
           1 + ob[,(m+1):(m+n),drop=FALSE]*dir[,(m+1):(m+n),drop=FALSE]/Y)
      } else {
         warning("Illegal ORIENTATION for argument DIRECT") 
      }
 # print("eff")
 # print(eff)
      if ( class(eff)=="matrix" && ( dim(eff)[1]==1 || dim(eff)==1 ) )
         eff <- c(eff) 
   }

   # undgaa afrundingsfejl i e naar den er taet ved 1.
   lpcontr <- lp.control(lps)
   eps <- sqrt(lpcontr$epsilon["epsint"])
   eff[abs(eff-1) < eps] <- 1

   # Slut med at bruge lps
   # delete.lp(lps)


   colnames(u) <- paste("u",1:m,sep="")
   colnames(v) <- paste("v",1:n,sep="")

   if ( TRANSPOSE )  {
      u <- t(u)
      v <- t(v)
      if ( class(eff)=="matrix" )
         eff <- t(eff)
   }

   oe <- list(eff=eff, objval=objval, RTS=RTS,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE,
              u=u, v=v, gamma=gamma, sol=sol)
   class(oe) <- "Farrell"


   return (oe)
} ## dea.dual



########################

# X = x
# Y = y
# RTS="crs"
# ORIENTATION="in"
# XREF=NULL
# YREF=NULL
# FRONT.IDX=NULL
# DUAL=NULL
# TRANSPOSE=FALSE
# LP=F
# CONTROL=NULL
# LPK=NULL


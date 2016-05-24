tgetcal <-
function(tv, # an object of class tverify
                    quad=FALSE # choose linear or quadratic calibration
                  ) 
{
#
# first define some utility functions
#

getopt <- function(p,o)
{
 pmat <- p
 omat <- o
 initpars <- rep(0,12)
 T <- oldgetT(p,o)       # might want to change this
 initpars[2] <- T[1,1] # I think this is right
 initpars[3] <- T[1,3]
 initpars[8] <- T[3,1]
 initpars[9] <- T[3,3]   
 opt <- nlm(f=costf,p=initpars,pmat,omat,iterlim=1000)
 if (opt$code > 2)
 { 
   print("WARNING!! abnormal exit from nlm")
   print(opt)
 }
 opt
}

costf <- function(pars=rep(0,12),
                  pmat,
		          omat,
                  L = diag(c(1,1,1))/sqrt(2)
		 )
{
   phat <- quadf(pars,pmat)
   costf <- tscore(p=phat,o=omat,L=L) 
   costf   
}

getT <- function(p=cbind(.8,.1,.1),
                 o=cbind(1,0,0)
		 )
{
T <- t(solve(t(p)%*%p) %*% (t(p)%*%o))
return(T)		     
}


oldgetT <- function(p=cbind(.8,.1,.1),
                 o=cbind(1,0,0),
		 L = diag(c(1,1,1))/sqrt(2) # Brier score by default
		 )
{
#print("oldgetT")
#print(L)
n <- nrow(p)                  # number of points in sample
T <- array(NA,dim=c(3,3))     # initialise transformation matrix		     

t1 <- t(p %*% t(L) ) %*% (p %*% t(L))  #
t2 <- t(p %*% t(L) ) %*% (o %*% t(L))

#
# now perform "stacking" compatible with using 9-by-1 vector b
# in place of 3-by-3 matrix T
#

D    <- rbind(
              cbind(1*t1,0*t1,0*t1),
              cbind(0*t1,1*t1,0*t1),
              cbind(0*t1,0*t1,1*t1)
             )
	     
d    <- as.vector(t2)

#
# solve.QP finds vector b that minimises:
#
# -t(d) %*% b + 0.5 * t(b) %*% D %*% b
#
# subject to constraints
#
# t(A) %*% b >= b_0
#
# where the first meq of these constraints are strict equalities
#
#
if (!is.na( det(D)))
{    
      temp <- solve.QP(Dmat = D,
                       dvec = d,
		       Amat = t(rbind(
 		                      cbind(-1, 0, 0,-1, 0, 0,-1, 0, 0),
 				      cbind( 0,-1, 0, 0,-1, 0, 0,-1, 0),
 				      cbind( 0, 0,-1, 0, 0,-1, 0, 0,-1),
				      diag(rep(1,9)),
		                      diag(rep(-1,9))
				      )
				      ),
		       bvec = c(
 		                 rep(-1,3),   # sum to 1
				 rep( 0,9),   # all >= 0
				 rep(-1,9)    # all <= 1
				),
		       meq  = 3)              # first 3 are equalities
#
# now "unstack" by getting 3-by-3 matrix T from 9-by-1 vector b
#		       
      T <- matrix(temp$solution,nrow=3,ncol=3,byrow=TRUE)

   for (i in 1:3)
   {
      T[,i] <- tscale(T[,i]) # make sure T is a proper transition matrix
   }
}
return(T)		     
}

linf <- function(pars,pmat)
{
   p <- pmat
   phat   <-  p %*% t(pars)
   phat # transformed p vector
}

quadf <- function(pars,pmat)
{
   p <- pmat
   if (length(pars) != 12) stop("wrong dimensions for param vector pars in quadf")
   phat   <- NA * p
   
#
# this is really a 2d problem so we transform to input space (pi1,pi3) and output
# space (pi1hat,pi3hat)
#      
   pi1    <- rep(NA,nrow(p)) # work vectors in the nonlinear function
   pi3    <- rep(NA,nrow(p))
   pi1hat <- rep(NA,nrow(p))
   pi3hat <- rep(NA,nrow(p))

   eps <- 1.0e-12 # a small number
   for (i in 1:nrow(p))
   {
      for (j in 1:3)
      {
        if (p[i,j] <= 0) p[i,j] <- eps
	if (p[i,j] >= 1) p[i,j] <- 1-eps
      }
      p[i,] <- tscale(p[i,])
   } 
   pi1 <- p[,1]
   pi3 <- p[,3]
#
# chose quadratic as functional form
# but could choose anything here really
#  
   pi1hat <- pars[1] +
             pars[2] * pi1 +
	     pars[3] * pi3 +
	     pars[4] * pi1 * pi1 + 
	     pars[5] * pi1 * pi3 +
	     pars[6] * pi3 * pi3
	     
   pi3hat <- pars[7]  +
             pars[8]  * pi1 +
	         pars[9]  * pi3 +
	         pars[10] * pi1 * pi1 +
	         pars[11] * pi1 * pi3 +
	         pars[12] * pi3 * pi3
	    
   phat[,1] <- pi1hat
   phat[,3] <- pi3hat
   phat[,2] <- 1 - phat[,1] - phat[,3]	
   
   for (i in 1:nrow(phat))
   {
     phat[i,] <- tscale(phat[i,])
   }       	    
   phat # transformed p vector
}
   calhexc <- NA * tv$hexc # initialise calibrated hexagon corners
   rel     <- NA * tv$rel
   score   <- NA * tv$score
  
   if (quad)
   { 
      opt     <- getopt(p=tv$pk,o=tv$ok)  # calibrate on raw binned forecasts
      pars    <- opt$estimate               # get parameters 
      phat    <- quadf(pars,tv$pbin) # transform the binned forecasts
      outf    <- quadf      
   }
   else
   {
      pars <- oldgetT(p=tv$pk,o=tv$ok) 
      phat <- linf(pars,tv$pbin) # transform the binned forecasts
      outf <- linf
      opt  <- "linear fit"
   }
   
   

   

   for (c in 1:nrow(tv$pbin))
   {      
      rel[c] <- tscore(p=matrix(phat[c,],     nrow=1,ncol=3),
                       o=matrix(tv$obar[c,],nrow=1,ncol=3),
   		               L=tv$L)  		    
      score[c]     <- tv$unc[c]+rel[c]-tv$res[c]
 		                     
      if (quad) {calhexc[c,,] <- quadf(pars,tv$hexc[c,,])	}
         else   {calhexc[c,,] <- linf(pars,tv$hexc[c,,])}	                     
   }

unc <- tv$unc
res <- tv$res
Nobs <- tv$Nobs

relbar   <- weighted.mean(rel,Nobs,na.rm=TRUE)    # mean reliability
resbar   <- weighted.mean(res,Nobs,na.rm=TRUE)    # mean resolution 
scorebar <- weighted.mean(score,Nobs,na.rm=TRUE)  # mean score
uncbar   <- weighted.mean(unc,Nobs,na.rm=TRUE)    # mean score

if (quad) {f <- quadf}
else      {f <- linf }

   out <- list(pbin=phat,
               Nobs=Nobs,
	           obar=tv$obar,
	           score=score,        # score
	           unc=unc,            # uncertainty
	           rel=rel,            # reliability
	           res=res,            # resolution
	           scorebar=scorebar,        # score
	           uncbar=uncbar,            # uncertainty
	           relbar=relbar,            # reliability
	           resbar=resbar,            # resolution
	           ncirc=tv$ncirc,        # number of points along side
               p=tv$pbin,
	           o=tv$obar,
	           assigned=tv$assigned,  # index of pbin assigned to p
               L=tv$L,
		       hexc=calhexc,
               q=tv$q,
               pk = NA,
               ok = NA,
		       pars = pars,
		       opt = opt,
               f = f
	          )
   class(out) <- "tverify"
   return(out)
}

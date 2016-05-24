tverify.default <-
function(p=cbind(1,1,1)/3,
                            o=cbind(0,0,1),
                            q=cbind(1,1,1)/3,
                            forceq = FALSE,
		                    ncirc = 11,
		                    L = diag(c(1,1,1))/sqrt(2)    # score matrix is Brier by default
		                   )
{
n <- nrow(p)
if (!forceq) q <- t(tscale(colSums(o)))
pbin     <- tgrid(ncirc)                  # binned locations for p
obar     <- matrix(0,nrow=nrow(pbin),ncol=3) # initialise vector of mean obs
dens     <- matrix(0,nrow=nrow(pbin),ncol=1) # initialise vector of prediction densities
score    <- matrix(0,nrow=nrow(pbin),ncol=1) # initialise vector of scores using p
res      <- matrix(0,nrow=nrow(pbin),ncol=1) # initialise vector of resolution
rel      <- matrix(0,nrow=nrow(pbin),ncol=1) # initialise vector of reliability
unc      <- matrix(0,nrow=nrow(pbin),ncol=1) # initialise vector of uncertainties
assigned <- rep(NA,n)

#
# assign predictions to bins
#

for (i in 1:n)
{
  p[i,]    <- tscale(p[i,]) # ensure all forecasts are scaled
  mindiff  <- 1
  jbest    <- NA
  cand     <- matrix(0,ncol=3,nrow=6) # the three nirc*(candidate vectors) at corners
  n1       <- floor(ncirc*p[i,1])
  n2       <- floor(ncirc*p[i,2])
  cand[1,] <- c(n1,n2,ncirc-n1-n2) 

  n1 <- floor(ncirc*p[i,1])
  n3 <- floor(ncirc*p[i,3])
    cand[2,] <- c(n1,ncirc-n1-n3,n3) 
  n2 <- floor(ncirc*p[i,2])
  n3 <- floor(ncirc*p[i,3])	      
    cand[3,] <- c(ncirc-n2-n3,n2,n3) 	      
  n1 <- ceiling(ncirc*p[i,1])
  n2 <- min(ceiling(ncirc*p[i,2]),ncirc-n1)
    cand[4,] <- c(n1,n2,ncirc-n1-n2)
  n1 <- ceiling(ncirc*p[i,1])
  n3 <- min(ceiling(ncirc*p[i,3]),ncirc-n3)
    cand[5,] <- c(n1,ncirc-n1-n3,n3)
  n2 <- ceiling(ncirc*p[i,2])
  n3 <- min(ceiling(ncirc*p[i,3]),ncirc-n2)
    cand[6,] <- c(ncirc-n2-n3,n2,n3)
  index <- rep(NA,6)
  for (j in 1:6)
  {                    # find the index of the candidate vectors in the tgrid 
    index[j] <- 1 + 
                ((2*ncirc+3)/2)*cand[j,1] -
	        (1/2) * cand[j,1]*cand[j,1] +
	        cand[j,2]
    if (index[j] < 1 || index[j] > 0.5*(ncirc+1)*(ncirc+2)) 
       {print("");print(cand);print(index)}
  }	      
	      
  for (j in 1:6)
  { 
    diff <- max(c(p[i,1]-pbin[index[j],1],
                  p[i,2]-pbin[index[j],2],
                  p[i,3]-pbin[index[j],3]
		 )
	        )
				 
    if (diff < mindiff) {mindiff <- diff; jbest <- index[j]}		 
  }

  assigned[i]  <- jbest 
  dens[jbest]  <- dens[jbest]  + 1
  obar[jbest,] <- obar[jbest,] + o[i,]
  
  score[jbest] <- score[jbest] + tscore(p=matrix(pbin[jbest,],nrow=1,ncol=3),
                                        o=matrix(o[i,],       nrow=1,ncol=3),L=L) 
					
  unc[jbest]   <-   unc[jbest] + tscore(p=q,
                                        o=matrix(o[i,],nrow=1,ncol=3),L=L)
}      

#
# normalise counts of observations to get obar
#
  for (c in 1:nrow(pbin))
  {
    if (dens[c] > 0)
    {
     obar[c,] <- obar[c,]  / dens[c]   
     score[c] <- score[c]  / dens[c]  # get mean score
       unc[c] <-   unc[c]  / dens[c]  # get mean uncertainty
    }
    else
    {
     obar[c,] <- NA  
     score[c] <- NA
     unc[c]   <- NA
    }
  } 

for (c in 1:nrow(pbin))
{
   obarl=matrix(obar[c,],nrow=1,ncol=3)
   pbinl=matrix(pbin[c,],nrow=1,ncol=3)
   
   qloc       <- 0 * pbinl #initialise
   qloc[,1:3] <- q
      
   rel[c] <- tscore(p=pbinl,o=obarl,L=L)                   
   res[c] <- tscore(p=qloc, o=obarl,L=L)
}

eps <- 1/(ncirc)
eps <- 2*eps/3
bigN   <- (ncirc+1)*(ncirc+2)/2
d1 <- cbind(eps,-eps/2,-eps/2)
d2 <- cbind(-eps/2,eps,-eps/2)
d3 <- cbind(-eps/2,-eps/2,eps)

hexc <- array(NA,dim=c(bigN,6,3))   # array containing corners of hexagons centred on pbin
for (i in 1:bigN)
{
   hexc[i,1,] <- pbin[i,] + d2
   hexc[i,2,] <- pbin[i,] - d3
   hexc[i,3,] <- pbin[i,] + d1
   hexc[i,4,] <- pbin[i,] - d2
   hexc[i,5,] <- pbin[i,] + d3
   hexc[i,6,] <- pbin[i,] - d1
}
   eps <- 1.0e-12 # a small number
   for (i in 1:bigN)
   {  
     for (k in 1:6)
     {
      for (j in 1:3)
      {
        if (hexc[i,k,j] <= 0) hexc[i,k,j] <- eps
	if (hexc[i,k,j] >= 1) hexc[i,k,j] <- 1-eps
      }
      hexc[i,k,] <- tscale(hexc[i,k,])
     }
   }

score    <- matrix(score,nrow=nrow(pbin),ncol=1)
rel      <- matrix(rel,nrow=nrow(pbin),ncol=1)
res      <- matrix(res,nrow=nrow(pbin),ncol=1)
Nobs     <- matrix(dens,nrow=nrow(pbin),ncol=1)

relbar   <- weighted.mean(rel,Nobs,na.rm=TRUE)    # mean reliability
resbar   <- weighted.mean(res,Nobs,na.rm=TRUE)    # mean resolution 
scorebar <- weighted.mean(score,Nobs,na.rm=TRUE)  # mean score
uncbar   <- weighted.mean(unc,Nobs,na.rm=TRUE)    # mean score

pk <- NA * p
ok <- NA * p

for (i in 1:nrow(p))
{
  pk[i,] <- pbin[assigned[i],]
  ok[i,] <- obar[assigned[i],]
}

f <- function(pars=NA,p=cbind(1,1,1)/3)
{
   print("identity transformation: no calibration has been performed.")
   out <- p
   out
}

out <- list(pbin=pbin,
           Nobs=Nobs,
	       obar=obar,
	       score=score,        # score
	       unc=unc,            # uncertainty
	       rel=rel,            # reliability
	       res=res,            # resolution
	       scorebar=scorebar,  # score
	       uncbar=uncbar,      # uncertainty
	       relbar=relbar,      # reliability
	       resbar=resbar,      # resolution
	       ncirc=ncirc,                                 # number of points along side
           p=p,
	       o=o,
	       assigned=assigned,                            # index of pbin assigned to p
           L=L,
	       hexc = hexc,
           q = q,
           pk = pk,
           ok = ok,
           pars = NA,
           opt = NA,
           f = f
	      )
class(out) <- "tverify"
out
}

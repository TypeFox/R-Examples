CARrampsOcl.fit <-
function( alpha, beta, Q, y,  nsamp, seed, fixed=FALSE,
coefs = FALSE, randeffs = FALSE,  
designMat = NULL, mult=20, filename="params.txt")
{
# alpha = c(alphae, alphaz_1 to alphaz_{F-1})
# beta = c(betae, betaz_1 to betaz_{F-1})
# Q = list of structure matrices 
#      data have to be ordered correctly; spell this out
# coefs:  if true, compute regressions coefficients   (design mat must be
#                           provided)
# randeffs:  if true, compute random effects (the phis)
# preds:  if true, compute predicted values (not implemented yet)
# designMat: regression design mat; required if coefficients are requested
#              must be in null space of kronecker sum
###############################################################################
multicoreflag = 0
gpuflag = 2
preds = FALSE
nodenames = "localhost"
gputoolsflag = 0

# initiate OpenCL and set up contexts and kernels

#require(OpenCL)
plat <- oclPlatforms()
dev <- oclDevices( plat[[1]] )

if(!any(grepl("cl_khr_fp64", oclInfo(dev[[1]])$exts)))
    stop("GPU with double precision and Open CL capabilities required.")


k.kronyD <- setup1(dev=dev)
k.kronyD3 <- setup2(dev=dev)
k.sampling <- setup3(dev=dev)


nQ <- length(Q)   # moved earlier
#n <- sapply( Q, function(l) nrow(l) )   # dim of each Q matrix; moved earlier
N <- length(y)

set.seed(seed)    # single random number seed

# diagonalize the kronecker sum
# first, eigendecompose each matrix in Q 

eigenQ <- list()
valid_types <- c("CAR1","RW1","Gen")
n <- rep(0,nQ)
for( i in 1:nQ )
{
    if( grep( Q[[i]]$type, valid_types ) ==1 && is.vector(Q[[i]]$content)
          && length(Q[[i]]$content)==2) 
        {
          eigenQ[[i]] <- eigenCAR1( Q[[i]]$content[1], Q[[i]]$content[2] )
          n[i] <- prod(Q[[i]]$content)
        }
    else  if( grep( Q[[i]]$type, valid_types ) == 2 && 
          is.numeric(Q[[i]]$content) && length(Q[[i]]$content)== 1)
        {
          eigenQ[[i]] <- eigenRW1( Q[[i]]$content )
          n[i] <- Q[[i]]$content
        }
    else  if( grep( Q[[i]]$type, valid_types ) == 3 && 
          is.matrix(Q[[i]]$content))
          {
          eigenQ[[i]] <- eigen( Q[[i]]$content )
          n[i] <- dim(Q[[i]]$content)[1]
          }
    else stop("invalid Q specification")
}

epseigen <- 0.00000000001

rm(Q)

# construct the kronecker sum, its eigenval matrix D,
#     and its eigenvect matrix BT

By <- rep(0,N)
 
    D <- matrix( eigenQ[[1]]$values, ncol = 1)
    sofar <- 1

i <- 2
while (i <= nQ )
{
    sofar <- sofar * n[i-1]
    D <- cbind( kronecker( rep(1, n[i]), D),
                kronecker( eigenQ[[i]]$values, rep(1, sofar ) ))
    i <- i + 1
}
#print("done with computing D")

F <- ncol(D) + 1  
sums <- apply(D,1,sum)


k <- length(sums[abs(sums) < epseigen ])   # rank deficiency of sum of Qs

nullSpaceIndex <- (1:nrow(D))[abs(sums) < epseigen ] 



if(coefs)
{
# change 11/02/11
   # if( !( length(nullSpaceIndex) == ncol(X) ))
    if( !( length(nullSpaceIndex) == ncol(designMat) ))
        stop("Dimension of covariate matrix does not match dimension of
             null space of structure matrices")
	nullSpaceEigenvects <- matrix(0, nrow = N, 
                ncol = length(nullSpaceIndex))
#   print(dim(nullSpaceEigenvects))
}

sumalpha <- sum(alpha)

# multiply B %*% y
# and, if coefs is true, build matrix of eigenvectors in null space of
#           kronecker sum of Q's

if( nQ == 1)
   {
 
    By <- t(eigenQ[[1]]$vectors) %*% y
    if(coefs)
        {
        #print(nullSpaceIndex)
        nullSpaceEigenvects <- eigenQ[[1]]$vectors[,nullSpaceIndex]
        if(is.vector(nullSpaceEigenvects) )
                   nullSpaceEigenvects <- matrix( nullSpaceEigenvects,
               ncol=1) 
        }
    } 

else if (nQ == 2)
{

    if(gpuflag==2)
        {
        # changed for OpenCL MKC 09/06/12
        By <- as.vector(oclKronVectMult1col( k.kronyD, t(eigenQ[[2]]$vectors), 
              t(eigenQ[[1]]$vectors),  y))
        }
    else
    {
        for( ind in 1:N)
             {
               tmp <- ind %% n[1]
               j <- tmp + as.numeric(!(tmp)) * n[1]
               i <- (ind-j) / n[1] + 1
              vind <- kronecker(eigenQ[[2]]$vectors[,i],
                      eigenQ[[1]]$vectors[,j]) 
              if( coefs & ind %in% nullSpaceIndex)
                  # change 03/25/11
                  nullSpaceEigenvects[,which(nullSpaceIndex==ind)] <-  vind
              By[ind] <- inprod( y, vind)
             }
    }
}

else if (nQ == 3)
{
    if(gpuflag==2)
        {
        # changed for OpenCL MKC 09/06/12
        By <- as.vector(oclKronVectMult1col3Q( k.kronyD3, t(eigenQ[[3]]$vectors), 
              t(eigenQ[[2]]$vectors), t(eigenQ[[1]]$vectors), 
               y))
        }

  else
  {
   for( j in 1:n[1])
      for( i in 1:n[2])
      {
         vij <- kronecker(eigenQ[[2]]$vectors[,i], eigenQ[[1]]$vectors[,j])
         for( h in 1:n[3])
         {
          ind <- (h-1) * n[1] * n[2] + (i-1) * n[1] + j
          vind <- kronecker( eigenQ[[3]]$vectors[,h], vij)
          if( coefs & ind %in% nullSpaceIndex)
              # change 03/25/11
              nullSpaceEigenvects[,which(nullSpaceIndex==ind)] <-  vind
              #nullSpaceEigenvects <- cbind(nullSpaceEigenvects, vind)
          By[ind] <- inprod(y, vind)
         }
      }
   }
#print("done with computing By")
}

else if (nQ == 4)
{

   for( j in 1:n[1])
      for( i in 1:n[2])
      {
         vij <- kronecker(eigenQ[[2]]$vectors[,i], eigenQ[[1]]$vectors[,j])
         for( h in 1:n[3])
         {

            vhij <- kronecker( eigenQ[[3]]$vectors[,h], vij)
            for( g in 1:n[4])
            {
             ind <- (g-1) * n[1] * n[2] * n[3] +
                      (h-1) * n[1] * n[2] + (i-1) * n[1] + j
             vind <- kronecker( eigenQ[[4]]$vectors[,g], vhij)
             if( coefs & ind %in% nullSpaceIndex)
                  # change 03/25/11
                  nullSpaceEigenvects[,which(nullSpaceIndex==ind)] <-  vind
                  # nullSpaceEigenvects <- cbind(nullSpaceEigenvects, vind)
             By[ind] <- inprod(y, vind)
            }
         }
      }

}

if(coefs & gpuflag==2)   # if we haven't built nullSpaceEigenvects yet
{
if(nQ==2)
   {
    for( ind in nullSpaceIndex)
       {
               tmp <- ind %% n[1]
               j <- tmp + as.numeric(!(tmp)) * n[1]
               i <- (ind-j) / n[1] + 1
               vind <- kronecker(eigenQ[[2]]$vectors[,i],
                      eigenQ[[1]]$vectors[,j])
              nullSpaceEigenvects[,which(nullSpaceIndex==ind)] <-  vind
        }
   }
else if(nQ==3)
  {
    for( ind in nullSpaceIndex)
       {
       for( j in 1:n[1])
         for( i in 1:n[2])
         {
            vij <- kronecker(eigenQ[[2]]$vectors[,i], eigenQ[[1]]$vectors[,j])
            for( h in 1:n[3])
            {  
             ind <- (h-1) * n[1] * n[2] + (i-1) * n[1] + j
             vind <- kronecker( eigenQ[[3]]$vectors[,h], vij)
             nullSpaceEigenvects[,which(nullSpaceIndex==ind)] <-  vind
            }
         }

    
     }
  }
}
#print("done with computing By")
#print(system("date"))
#print(c("coefs",coefs))

if(coefs)
     {
     #print(designMat)
     #print(nullSpaceEigenvects)

     #print("dim nullSpaceEigenvects")
     #print(dim(nullSpaceEigenvects))

     A <- coefficients(lm( designMat ~ -1 + nullSpaceEigenvects ))

     same <- all.equal( designMat, nullSpaceEigenvects %*% A)
     if(!isTRUE(same))
             {
               print("regressions coefs invalid")
             }

     # so designMatrix = nullSpaceEigenvects %*% A
     # and regression coefs for designMat will be solve(A) 
     #            %*% Bphi[nullSpaceIndex]
     }



# find mode of log posterior
logmode <- optimizelogpost( alpha=alpha,beta=beta, D=D, y=y, By=By,  k=k,
func = multivspost2blog  )

print("done optimizing ")
#print(system("date"))

# rejection from uniform envelope

accepted <- matrix( nrow=0, ncol= F + 1)
batch = 0
nacpt = 0
while( nrow(accepted) < nsamp)
{

cands <- rdirichlet( mult * nsamp, rep(1,F)  )

unifs <- runif( mult * nsamp )

results <- oclSampling( k.sampling, smat = matrix(cands[,1:(F-1)], ncol=F-1), 
     alpha=alpha, beta=beta, D=D, By=By, k = k)


# results has logposterior for each s in first column, newbeta in 2nd col

logpost <- results[,1]

keep <- ( logpost > log(unifs) + logmode )  # check the log(F-1) part

# return matrix of accepted values with s in first col and corresponding tausqtot in 2nd


accepted <- rbind( accepted, cbind(cands, results[,2])[keep,])

if(nrow(accepted) > nacpt)
{
   write.table(accepted, file=filename )
}


batch <- batch + 1
print(paste("accepted so far", nrow(accepted)) )
nacpt <- nrow(accepted)

if(fixed)     # new 09/14/10
     break
}


print("done generating s")
print(system("date"))


# if the number of candidates isn't fixed, strip off extra accepted results


acptrate <- nrow(accepted) / (mult * nsamp * batch)  # first calc acptrate
if(!fixed)
    accepted <- accepted[1:nsamp, ]

newalpha <- sum(alpha) + (N-k)/2 
newbeta <- accepted[,F+1]
naccpt <- nrow(accepted)
accepted[,F+1] <- rgamma( naccpt, rep(newalpha, naccpt), newbeta )

tausqy <- (accepted[,F]) * accepted[,F+1]
tausqphi <- matrix(accepted[,1:(F-1)] * accepted[,F+1],ncol=F-1)


output <- as.matrix(cbind(accepted, tausqy, tausqphi))

output <- output[ !is.na(output[,1]), ]

tausqy <- tausqy[!is.na(output[,1])]
tausqphi <- tausqphi[!is.na(output[,1]), ]


# -----------------------------------------------------------
if(randeffs | preds )
{
#if(gpuflag==1)
 #  phi <- combo1colFromC(eigenQ[[2]]$vectors, eigenQ[[1]]$vectors, 
 #     D, tausqy, tausqphi, By)
#else
   if(gpuflag==2)
      {
      if(nQ==1)
         {
 #        print("executing 1Q on GPU")
         phi <- oclCombo1col1( eigenQ[[1]]$vectors,
            D, tausqy, tausqphi, By)
         }
      if(nQ==2)
         {
 #        print("executing 2Q on GPU")
         phi <- oclCombo1col(eigenQ[[2]]$vectors, eigenQ[[1]]$vectors,
            D, tausqy, tausqphi, By)
         }
      else
         if(nQ==3)
         {
        # print("executing 3Q on GPU")
         phi <- oclCombo1col3(eigenQ[[3]]$vectors, eigenQ[[2]]$vectors, 
            eigenQ[[1]]$vectors,
            D, tausqy, tausqphi, By)
         }
      }


}
else
    phi <- NULL

print("done generating phi")
print(system("date"))

# generate coefficients
# change 03/05/11 MKC

if(coefs)
    {
     ncoefs <- length( nullSpaceIndex )
     coefs1 <- matrix( 0, nrow=ncoefs, ncol=nrow(output))   # this is transpose of orig

     for( i in 1:nrow(output) )  {

        neweigendenom <- D  %*% 
               (output[i, (F+3):(2 * F + 1)]) + tausqy[i]

        coefs1[,i ] <-  ( rnorm( N, tausqy[i] * By / neweigendenom,
           1 /sqrt(neweigendenom) ))[nullSpaceIndex]    # transpose


        }

    regcoefs <- t(solve(A) %*% coefs1)
    }
else
    regcoefs <- NULL


if(preds)
{

}
else
   preds <- NULL



# -------------------------------------------------------------

colnames(output) <- c(paste("s",1:(F-1),sep=""), "s0", "tausqtot", "tausqy", 
paste("tausqphi",1:(F-1),sep=""))

numparms <- ncol(output)

#detach(kernels)
#rm(kernels)

list(params = output[ , (F+2):numparms], phi = phi, preds=preds, 
regcoefs = regcoefs, 
D=D, y=y, acptrate = acptrate, n=n )

}


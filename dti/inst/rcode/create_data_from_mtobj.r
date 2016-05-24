datamixtens <- function(mtobj,grad,sigma,bf=1,ns0=NULL,smoothS0=TRUE,maxs0=30000){
#
#  this is a mmodification of the function in
#  source(system.file("rcode/gen_mixtens.r",package="dti"))
#  that can handle gradients with norm <> 1
#
#  w  - matrix of wghts (m+1,n1,n2,n3)
#  th - ev-parameters (2,n1,n2,n3) lambda1 = th[1]+th[2]
#                                  lambda2 = th[2]
#  alpha, beta (m,n1,n2,n3)  spherical coordinates of directions
#
args <- list(sys.call())
m <- dim(mtobj@mix)[1]
ddim <- dim(mtobj@mix)[2:4]
w <- mtobj@mix
w0 <- 1-apply(w,2:4,sum)
ngrad <- dim(grad)[2]
d <- extract(mtobj,"andir")$andir
th <- mtobj@ev
order <- mtobj@order
s0 <- extract(mtobj,"s0")$s0
if(length(dim(s0))==4) s0 <- apply(s0,1:3,mean)
ms0 <- median(s0[mtobj@mask])
if(ms0 < 10000) s0 <- s0/ms0*10000
ds0 <- dim(s0)
if(is.null(ns0)) ns0 <- length(mtobj@s0ind)
if(!require(aws)) smoothS0 <- FALSE
if(smoothS0){
s0 <- aws(s0,hmax=4)@theta
}
grad <- cbind(matrix(0,3,ns0),grad)
ngrad <- ngrad + ns0
dim(s0) <- ds0
s0[s0>maxs0] <- maxs0
sb <- array(0,c(dim(s0),ngrad))
for(i in 1:ns0) sb[,,,i] <- s0
for(i in (ns0+1):ngrad){
#  initialize with isotropic part
   sb[,,,i] <- w0*exp(-th[1,,,]-th[2,,,])
   for(j in 1:m){
      dg <- d[1,j,,,]*grad[1,i]+d[2,j,,,]*grad[2,i]+d[3,j,,,]*grad[3,i]
      sb[,,,i] <- sb[,,,i]+w[j,,,]*exp(-th[2,,,]-th[1,,,]*dg^2)
   }
   sb[,,,i] <- sb[,,,i]*s0
}
sb <- as.integer(sqrt(rnorm(sb,sb,sigma)^2+rnorm(sb,0,sigma)^2))
dim(sb) <- c(prod(ddim),ngrad)
sb[!mtobj@mask,] <- 0
dim(sb) <- c(ddim,ngrad)
ddim0 <- ddim
invisible(new("dtiData",
                call = args,
                si     = sb,
                gradient = grad,
                btb    = create.designmatrix.dti(grad),
                ngrad  = as.integer(ngrad), # = dim(btb)[2]
                s0ind  = as.integer(1:ns0), # indices of S_0 images
                replind = as.integer(c(rep(1,ns0),(ns0+1):ngrad)),
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = 1:ddim[1],
                yind   = 1:ddim[2],
                zind   = 1:ddim[3],
                level  = 1,
                sdcoef = c(5.e-01, 0, 1, 2.e+04),
                voxelext = mtobj@voxelext,
                orientation = mtobj@orientation,
                rotation = diag(3),
                source = "artificial")
            )
}
create.designmatrix.dti <- function(gradient, bvalue=1) {
  dgrad <- dim(gradient)
  if (dgrad[2]==3) gradient <- t(gradient)
  if (dgrad[1]!=3) stop("Not a valid gradient matrix")

  btb <- matrix(0,6,dgrad[2])
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]

  btb * bvalue
}

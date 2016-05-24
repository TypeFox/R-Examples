tdatamixtens <- function(w,th,alpha,beta,grad,bvalue,sigma,ns0){
#  
#  w  - matrix of wghts (m+1,n1,n2,n3)
#  th - ev-parameters (2,n1,n2,n3) lambda1 = th[1]+th[2]
#                                  lambda2 = th[2]
#  alpha, beta (m,n1,n2,n3)  spherical coordinates of directions
#
args <- list(sys.call())
ddim <- dim(w)[2:4]
m <- dim(w)[1]-1
ngrad <- dim(grad)[2]
#  check array dimensions
if(any(dim(th)!=c(2,ddim))) return("wrong dimension of th")
if(any(dim(alpha)!=c(m,ddim))) return("wrong dimension of alpha")
if(any(dim(beta)!=c(m,ddim))) return("wrong dimension of beta")
if(any(dim(grad)!=c(3,ngrad))) return("wrong dimension of grad")
sb <- array(0,c(ddim,ngrad))
d <- array(0,c(3,m,ddim))
sal <- sin(alpha)
cal <- cos(alpha)
sbe <- sin(beta)
cbe <- cos(beta)
d[1,,,,] <- sal*cbe
d[2,,,,] <- sal*sbe
d[3,,,,] <- cal
s0 <- array(10000,ddim)
for(i in 1:ns0) sb[,,,i] <- s0
for(i in (ns0+1):ngrad){
   sb[,,,i] <- w[1,,,]*exp(-bvalue[i]*(th[1,,,]+th[2,,,]))
   for(j in 1:m){
      dg <- (d[1,j,,,]*grad[1,i]+d[2,j,,,]*grad[2,i]+d[3,j,,,]*grad[3,i])
      sb[,,,i] <- sb[,,,i]+w[j+1,,,]*exp(-bvalue[i]*(th[2,,,]+th[1,,,]*dg^2))
   }
   sb[,,,i] <- sb[,,,i]*s0
}
sb <- sqrt(rnorm(sb,sb,sigma*s0)^2+rnorm(sb,0,sigma*s0)^2)
dim(sb) <- c(ddim,ngrad)
ddim0 <- ddim
invisible(new("dtiData",
                call = args,
                si     = sb,
                gradient = grad,
                btb    = sweep( create.designmatrix.dti(grad), 2, bvalue, "*"),
                bvalue = bvalue,
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
                voxelext = c(1,1,1),
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                source = "artificial")
            )
}
truemixtens <- function(w,th,alpha,beta,grad,bvalue,sigma,ns0){
#  
#  w  - matrix of wghts (m+1,n1,n2,n3)
#  th - ev-parameters (2,n1,n2,n3) lambda1 = th[1]+th[2]
#                                  lambda2 = th[2]
#  alpha, beta (m,n1,n2,n3)  spherical coordinates of directions
#
args <- list(sys.call())
ddim <- dim(w)[2:4]
m <- dim(w)[1]-1
ngrad <- dim(grad)[2]
#  check array dimensions
if(any(dim(th)!=c(2,ddim))) return("wrong dimension of th")
if(any(dim(alpha)!=c(m,ddim))) return("wrong dimension of alpha")
if(any(dim(beta)!=c(m,ddim))) return("wrong dimension of beta")
if(any(dim(grad)!=c(3,ngrad))) return("wrong dimension of grad")
ddim0 <- ddim
orient <- array(0,c(2,m,ddim))
orient[1,,,,] <- alpha
orient[2,,,,] <- beta

  invisible(new("dwiMixtensor",
                model  = "homogeneous_prolate",
                call   = args,
                ev     = th,
                mix    = w[-1,,,,drop=FALSE],
                orient = orient,
                order  = array(m,ddim),
                p      = 1,
                th0    = array(10000,ddim),
                sigma  = array(10000*sigma,ddim),
                scorr  = array(0,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = array(TRUE,ddim),
                hmax   = 1,
                gradient = grad,
                btb    = sweep( create.designmatrix.dti(grad), 2, bvalue, "*"),
                bvalue = bvalue,
                ngrad  = as.integer(ngrad), # = dim(btb)[2]
                s0ind  = as.integer(1:ns0),
                replind = as.integer(c(rep(1,ns0),(ns0+1):ngrad)),
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = 1:ddim[1],
                yind   = 1:ddim[2],
                zind   = 1:ddim[3],
                voxelext = c(1,1,1),
                level = 1,
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                source = "artificial",
                outlier = numeric(0),
                scale = 1,
                method = "mixtensor")
            )
}
create.designmatrix.dti <- function(gradient) {
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
  btb
}
odfdist <- function(obj1,obj2,poly=4,mask=NULL){
  if(any(obj1@ddim!=obj2@ddim)) return(warning("incompatible dimensions"))
  if(!class(obj1)%in%c("dwiMixtensor","dtiTensor","dwiQball")) return(warning("incompatible class of obj1"))
  if(!class(obj2)%in%c("dwiMixtensor","dtiTensor","dwiQball")) return(warning("incompatible class of obj2"))
  if(!exists("icosa0")) data("polyeders")
  n <- prod(obj1@ddim)
  if(n>419103) poly <- min(poly,3)
  if(n>1672495) poly <- min(poly,2)
  if(n>6628036) poly <- min(poly,1)
  polyeder <- switch(poly+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  if(is.null(mask)) mask <- rep(TRUE,n)
  if(class(obj1)=="dwiQball"){
     sphdesign <- design.spheven(obj1@order,polyeder$vertices,obj1@lambda)$design
     sphcoef <- obj1@sphcoef
     dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
  }
  radii1 <- if(class(obj1)=="dwiMixtensor") {
               .Fortran("mixtradi",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj1@ev),
                    as.double(obj1@orient),
                    as.double(obj1@mix),
                    as.integer(obj1@order),
                    as.integer(dim(obj1@mix)[1]),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii 
                } else if(class(obj1)=="dtiTensor"){
                .Fortran("odfradii",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj1@D),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
                } else t(sphdesign)%*%sphcoef
  dim(radii1) <- c(polyeder$nv,prod(obj1@ddim))
  if(class(obj2)=="dwiQball") {
     sphdesign <- design.spheven(obj2@order,polyeder$vertices,obj2@lambda)$design
     sphcoef <- obj2@sphcoef
     dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
  }
  radii2 <- if(class(obj2)=="dwiMixtensor") {
               .Fortran("mixtradi",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj2@ev),
                    as.double(obj2@orient),
                    as.double(obj2@mix),
                    as.integer(obj2@order),
                    as.integer(dim(obj2@mix)[1]),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
                } else if(class(obj2)=="dtiTensor"){
                .Fortran("odfradii",
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(obj2@D),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    DUP=FALSE,
                    PACKAGE="dti")$radii
                } else t(sphdesign)%*%sphcoef
  dim(radii2) <- c(polyeder$nv,prod(obj1@ddim))
  t((radii2-radii1)[,mask]^2)%*%rep(1/polyeder$nv,polyeder$nv)
  }

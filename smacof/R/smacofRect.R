`smacofRect` <- function(delta, ndim = 2, circle = c("none", "row", "column"), weightmat = NULL, init = NULL, verbose = FALSE,
                       itmax = 1000, reg = 1e-6, eps = 1e-6)

# init ... either a list of 2 matrices of dimension n \times p and m \times p with 
# starting values. if NULL, svd is used.
{
  
  circle <- match.arg(circle, c("none", "row", "column"), several.ok = FALSE)
  diss <- delta
  rnames <- rownames(delta)
  if (is.data.frame(diss)) diss <- as.matrix(diss)
  checkdiss(diss)
  
  n <- dim(diss)[1]                       #number of individuals
  m <- dim(diss)[2]                       #number of objects
  p <- ndim

  if (is.null(weightmat)) {
    w <- matrix(1,n,m)                    #initialize weights (as 1)
  } else w <- weightmat
  

  itel <- 1
  delta <- ifelse(is.na(diss),0,diss)     #replace NA's by 0
  #delta <- delta/sqrt(sum(w*delta^2))*sqrt(n*m)       #normalize dissimilarities
  
  delta_plus <- ifelse(delta>=0,delta,0)  #delta decomposition (+)
  delta_min <- ifelse(delta<=0,-delta,0)  #delta decomposition (-) (if all >0 --> complete 0)

  if (is.list(init)) {
    x <-init[[1]]                         #list as input structure
    y <-init[[2]]
  } else {
    e <- delta_plus^2
    e <- -0.5*(e-outer(rowSums(e)/m,colSums(e)/n,"+")+(sum(e)/(n*m)))
    #e <- e/sqrt(sum(e^2))*sqrt(n*m) 
    
    z <- svd(e,nu=p,nv=0)                 #SVD for e (pos. distances)
    x<-z$u                                #starting value for x
    y<-crossprod(e,x)                     #starting value for y
  }
  
  if (circle != "none"){
    r <- projCircle(x,y,x,y,circle=circle)
    x <- r$x
    y <- r$y
    wr <- rowSums(w)
    wc <- colSums(w)
    lambda <- 2*max(c(wr,wc))    
  }
  
  d <- distRect(x,y,reg)                  #n times m of reproduced diss
  lold <- sum(w*(delta-d)^2)              #stress value
  

  
  #------------------- begin majorization -----------------------------
  repeat {
    if (circle == "none") {
	      ww <- w*(1+(delta_min/d)) 
        wr <- rowSums(ww)
        wc <- colSums(ww)

        v <- solve(diag(wc)+(1/m) - crossprod(ww,ww/wr)) -(1/m)

        b <- w*delta_plus/d                #B matrix
        br <- rowSums(b)                   #rows B
        bc <- colSums(b)                   #columns W
       
        xraw <- (br*x)-(b%*%y)
        yraw <- (bc*y)-crossprod(b,x)
        
        y <- v%*%(yraw+crossprod(ww,xraw/wr)) #x update 
        x <- (xraw+(ww%*%y))/wr               #y update
        
     } else {
       b  <- w*(1-delta/d)               #B matrix
       br <- rowSums(b)                   #rows B
       bc <- colSums(b)                   #columns W
       xunc <- x - outer(br,rep(1/lambda,p),"*")*x + b%*%(y/lambda)
       yunc <- y - outer(bc,rep(1/lambda,p),"*")*y + t(b)%*%(x/lambda)
       r <- projCircle(xunc,yunc,x,y,circle=circle)
       x <- r$x
       y <- r$y
     }  

    d <- distRect(x,y,reg)             #compute distances (update)
    
    lnew <- sum(w*(delta-d)^2)         #compute stress

    if (verbose) cat("Iteration: ",formatC(itel,digits=6,width=6), "   Stress (not normalized):",formatC(lnew,digits=6,width=12,format="f"),"   Dif:",formatC(lold-lnew,digits=6,width=12,format="f"),"\n")

    if (((lold-lnew) < eps) || (itel==itmax)) break() 
        
	  lold <- lnew                       #update stress
    itel <- itel+1
  }
  #-------------------- end majorization --------------------------

colnames(y) <- colnames(x) <- paste("D",1:(dim(y)[2]),sep="")
rownames(x) <- rownames(diss) <- rownames(d) <- rnames

# point stress 
resmat <- as.matrix(d - diss)^2    #point stress
spp.col <- colMeans(resmat, na.rm = TRUE)
spp.col <- spp.col/sum(spp.col)*100
spp.row <- rowMeans(resmat, na.rm = TRUE)
spp.row <- spp.row/sum(spp.row)*100

if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")

## stress normalization
lnew <- sqrt(sum(w*(diss-d)^2, na.rm = TRUE)/sum(d^2, na.rm = TRUE))
  
## congruence coefficients
diss0 <- diss
diss0[is.na(diss0)] <- 0
congnum <- diag(diss0 %*% t(d))
congdenom <- sqrt(diag(diss0 %*% t(diss0)) * diag(d %*% t(d)))
congvec <- congnum/congdenom
  
#return configuration distances, row and column configurations, stress 
result <- list(obsdiss = diss, confdiss = d, conf.row = x, conf.col = y, stress = lnew, 
               spp.row = spp.row, spp.col = spp.col, congvec = congvec,
               ndim = p, model = "Rectangular smacof", niter = itel, nind = n, nobj = m, call = match.call()) 
class(result) <- "smacofR"
result 
}

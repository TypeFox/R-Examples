qskellam <- function(p, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE){
 # inverse CDF of Skellam distriubition (difference of Poissons)
    if (missing(p)|missing(lambda1)) stop("first 2 arguments are required")
   # make all args the same length (for subsetting)
    lens <- c(length(p),length(lambda1),length(lambda2))
    len <- max(lens)
    if(len>min(lens)) {
        if (all(len/lens!=len%/%lens)) warning("longer object length is not a multiple of shorter object length", domain=NA)
        p <- rep(p,length.out=len)
        lambda1 <- rep(lambda1,length.out=len)
        lambda2 <- rep(lambda2,length.out=len)
    }
    ret <- rep(NaN,length.out=len)
   # handle a zero lambda separately (quicker than interpreted search)
    if (getOption("verbose")){  # set via  options(verbose=TRUE)  # could create a new option called validate, but would would have to test for existence before value
        nz<-rep(TRUE,length.out=len)    # verify search by using it for Poisson too
    } else {
        ret[lambda2==0] <- qpois(p[lambda2==0],lambda1[lambda2==0],lower.tail=lower.tail,log.p=log.p)
        ret[lambda1==0] <- -qpois(p[lambda1==0],lambda2[lambda1==0],lower.tail=!lower.tail,log.p=log.p)
        nz <- (lambda1!=0)&(lambda2!=0)
    }
   # handle boundaries correctly
    bdry <- nz & ((p==0) | (p+1.01*.Machine$double.eps>=1))     # match qpois in assuming that p with 2.25e-16 of 1 are actually 1
    if(any(bdry)){
        if(lower.tail){
            ret[p[bdry]==0] <- ifelse(lambda2[p[bdry]==0]==0,0,-Inf)
            ret[p[bdry] >0] <- ifelse(lambda1[p[bdry] >0]==0,0, Inf)
        } else {
            ret[p[bdry] >0] <- ifelse(lambda2[p[bdry]==0]==0,0,-Inf)
            ret[p[bdry]==0] <- ifelse(lambda1[p[bdry] >0]==0,0, Inf)
        }
    }
    if(any(bdry)|(!getOption("verbose")&any(!nz))){
       # avoid repeated susetting later
        nz <- nz & !bdry
        p <- p[nz]
        lambda1 <- lambda1[nz]
        lambda2 <- lambda2[nz]
    }
   # Cornish-Fisher approximations
    z <- qnorm(p,lower.tail=lower.tail,log.p=log.p)
    mu <- lambda1-lambda2
    vr <- lambda1+lambda2
    sg <- sqrt(vr)
    c0 <- mu+z*sg
    c1 <- (z^2-1)*mu/vr/6
    c2 <- -( c1*mu - 2*lambda1*lambda2*(z^2-3)/vr )*z/12/vr/sg
   # test and linear search (slow if p extreme or lambda1+lambda2 small)
    q0 <- round(c0+c1+c2)
    p0 <- pskellam(q0,lambda1,lambda2,lower.tail=lower.tail,log.p=log.p)
    p <- p*(1-64*.Machine$double.eps)       # match qpois in assuming that a value within 1.4e-14 of p actually equals p
    if (lower.tail){    # smallest x such that F(x) >= p (consider floor for greater efficiency?)
        up <- up1 <- p0 < p
        while(any(up1)){
            q0[up1] <- q0[up1]+1
            up1[up1] <- pskellam(q0[up1],lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) < p[up1]
        }
        down1 <- (!up)&(p0>p)   # oops: \qskellam(p,lambda,0) == qpois(p,lambda)+1 usually
        down1[down1] <- (pskellam(q0[down1]-1,lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1])
        while(any(down1)){
            q0[down1] <- q0[down1]-1
            down1[down1] <- (pskellam(q0[down1]-1,lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1])
        }
    } else {        # largest x such that F(x,lower.tail=FALSE) <= p (consider ceiling for greater efficiency?)
        down <- down1 <- p0 > p
        while(any(down1)){
            q0[down1] <- q0[down1]+1
            down1[down1] <- pskellam(q0[down1],lambda1[down1],lambda2[down1],lower.tail=lower.tail,log.p=log.p) > p[down1]
        }
        up1 <- (!down)&(p0<p)
        up1[up1] <- !(pskellam(q0[up1]-1,lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) > p[up1])
        while(any(up1)){
            q0[up1] <- q0[up1]-1
            up1[up1] <- !(pskellam(q0[up1]-1,lambda1[up1],lambda2[up1],lower.tail=lower.tail,log.p=log.p) > p[up1])
        }
    }
    ret[nz] <- q0
    ret
}

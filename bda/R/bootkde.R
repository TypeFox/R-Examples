
.bootkde <- function(x, from, to, gridsize=512L,
                    method='quantile', conf.level)
{
    stopifnot(class(x)=='bdata')
    if(x$top.coded!=0)
        stop("Not applicable to top-coded data")
    
    method <- match.arg(tolower(method),c("z.score","quantile","cdf"))

    X <- x$mids;
    N <- length(X)
    F <- x$counts;
    B <- diff(x$breaks)/2;
    A <- -B;
    a <- x$breaks[1];
    b <- rev(x$breaks)[1]
    ucb=NULL; lcb=NULL;
    ## Set default bandwidth
    h0 = bw.nrd(rep(X,F)+runif(sum(F),-B,B))
    h = .Fortran(.F_hbmise,as.double(X), as.double(F),as.double(2*B),
        as.integer(N), hopt=as.double(h0))$hopt

    bw <- diff(x$breaks)
    
    a <- ifelse(missing(from), x$breaks[1] - bw[1], from)
    b <- ifelse(missing(to), x$breaks[N+1] + bw[N], to)
    gpoints <- seq(a, b, length = gridsize)
    
    y <- .Fortran(.F_ofcpdf,
                  as.double(X), as.double(F),as.double(-B),as.double(B),
                  as.integer(N), y=as.double(gpoints), as.integer(gridsize),
                  para=as.double(h))$y
    y <- cumsum(y)

    if(missing(conf.level)){
        conf.level <- NULL
        y <- .bootemp(sum(F),x=gpoints,y=X,f=F,lb=A,ub=B,Fx=y,from=a,to=b)
    }else{
        stopifnot(is.numeric(conf.level))
        stopifnot(conf.level>0&&conf.level<1)
        alpha <- 1.0 - conf.level
        iter <- 1000
        tmp <- as.matrix(rep(sum(F),iter), ncol=1);
        ## simulate the rounding errors by a pilot estimate of F(x) based on a BME
        out <- apply(tmp,1, .bootemp,x=gpoints,y=X,f=F,lb=A,ub=B,
                     Fx=y,from=from,to=to)
        
        if(method=="z.score"){
            sigs <- apply(out,1,sd);
            z0 <- qnorm(1-alpha/2.);
            ym <- apply(out,1,median);
            lcb <- ym-sigs*z0;
            lcb[lcb<0] <- 0;
            ucb <- ym+sigs*z0;
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Density"
        }else if(method=="quantile"){
            ym <- apply(out,1,median);
            lcb <- as.numeric(apply(out,1,quantile, alpha/2));
            lcb[lcb<0] <- 0;
            ucb <- as.numeric(apply(out,1,quantile, 1-alpha/2));
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Density"
        }else{
            out <- apply(out,2,cumsum)*(b-a)/gridsize
            ym <- apply(out,1,median);
            lcb <- as.numeric(apply(out,1,quantile, alpha/2));
            lcb[lcb<0] <- 0;
            ucb <- as.numeric(apply(out,1,quantile, 1-alpha/2));
            ucb[ucb>1] <- 1
            y <- out[,1];
            type <- "Probability"
        }
    }
    
    f0 <- y
    x0 <- gpoints;

    pars <- list(npar=NULL)
    return(structure(list(y=y,x=gpoints,
                          ucb=ucb,lcb=lcb, conf.level=conf.level,
                          call = match.call(),
                          name = x$xname
                          ), class = "histosmooth"))
}

.bootemp <- function(n,x,y,f,lb,ub,Fx,from,to){
    M=length(x);u=runif(n);N=length(f);

    smpl = .Fortran(.F_remp, as.integer(N), as.double(y), as.double(f),
        as.double(lb), as.double(ub), as.integer(M), as.double(Fx),
        as.double(x), smpl=as.double(u), as.integer(n))$smpl

    density(smpl,n=M,from=from,to=to)$y
}


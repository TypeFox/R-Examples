FastS <- function(x, y, Scontrol=list(int=FALSE, N=100, kk=2, tt=5, b=.5, cc=1.547,seed=NULL), beta_gamma) {
    # Fast-S algorithm for linear regression
    # This is a slightly changed version of the R-function fast.s published in
    # Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
    # Journal of Computational and Graphical Statistics, 15 (2), 414-427
    #
    # int = 1 -> add a column of ones to x
    # N = cant de sub-samples
    # kk = number of refining iterations in ea.
    #     subsample
    # kk = 0 means "raw-subsampling"
    # b = right hand side of the equation
    # cc = corresponding tuning constant
    # tt = number of "best betas" to remember
    #          from the subsamples. These will be later
    #          iterated until convergence

    if (is.null(Scontrol$int)) int <- FALSE else int <- Scontrol$int
    if (is.null(Scontrol$N)) N <- 100 else N <- Scontrol$N
    if (is.null(Scontrol$kk)) kk <- 2 else kk <- Scontrol$kk
    if (is.null(Scontrol$tt)) tt <- 5 else tt <- Scontrol$tt
    if (is.null(Scontrol$b)) b <- .5 else b <- Scontrol$b
    if (is.null(Scontrol$cc)) cc <- 1.547 else cc <- Scontrol$cc

    if (!is.null(Scontrol$seed))  seed  <- Scontrol$seed

    n <- dim(x)[1]

    if( int == 1) x <- cbind(rep(1,n), x)
       
    p <- dim(x)[2]       
               
    x <- as.matrix(x)
    x.red<- unique(x) # no double rows

    n.red<- nrow(x.red)
    if(n.red< n) {
        ff<- function(x) {
            x<-as.factor(x)
            levels(x)<- (1: length(levels(x)))
            return(x)
        }
        index<-apply(apply(x,2,ff),1,paste, collapse="_")
        index.red<- unique(index)  # index for different regressors
    }
	
    isStepmodel <- all(apply(x!=0,1,sum)==1)&(p>1)
	
    if(!is.null(Scontrol$seed)) set.seed(Scontrol$seed) 

    best.betas <- matrix(0, tt, p)
    best.scales <- rep(1e20, tt)
    s.worst <- 1e20
    n.ref <- 1

    if(!missing(beta_gamma)) N<- N+1 # beta_gamma is the N+1-th candidate

    rho <- function(u, cc=1.56) {
        w <- abs(u)<=cc
        v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
        v <- v*6/cc^2
        return(v)
    }
    loss.S <- function(u,s,cc) mean(rho(u/s,cc) )
    our.solve <- function(a,b) {
        a <- qr(a)
        da <- dim(a$qr)
        if(a$rank < (p <- da[2])) return(NA) else qr.coef(a, b)
    }

    stepsampler <- function(hnr, nichtnull) sample(rep(nichtnull[which(nichtnull[,2]==hnr),1],2),1)

    scale1 <- function(u, b, cc, initial.sc=median(abs(u))/.6745) {
        # find the scale, full iterations
        max.it <- 200
        # magic number alert
        #sc <- median(abs(u))/.6745
        sc <- initial.sc
        i <- 0
        eps <- 1e-20
        # magic number alert
        err <- 1
        while( ( (i <- i+1) < max.it ) && (err > eps) ) {
            sc2 <- sqrt( sc^2 * mean( rho( u / sc, cc ) ) / b   )
            err <- abs(sc2/sc - 1)
            sc <- sc2
        }
        return(sc)
    }
	
    for(i in 1:N) {
        #Generate candidates:
        if(!missing(beta_gamma)& i==1) {
            beta<- beta_gamma
        } else {
            if(isStepmodel) { # Get a subsample in general position, which means in this case: One Point from every step
                indices<- apply(cbind(1:p), 1, stepsampler, nichtnull=which(x!=0, arr.ind=TRUE))
                xs <- x[indices,]
                ys <- y[indices]
                beta <- our.solve(xs,ys)
            } else {
                # get a subsample
                singular <- TRUE
                itertest <- 1
                while(singular&& itertest<100) {
                    if(n.red<n) {
                        ranind <- sample(index.red, p, replace=FALSE) # choose p different regressor-label
                        zahlen <- which(index%in% ranind) # which indices belong to these?
                        i.red <- index[which(index%in% ranind)]  #which regressor-labels do they have
                        indices <- tapply(zahlen, INDEX=i.red, sample,1) # only one per regressor-label
                    } else {
                        indices <- sample(n,p)
                    }
                    xs <- x[indices,]
                    ys <- y[indices]
                    beta <- our.solve(xs,ys)
                    singular <- any(is.na(beta))
                    itertest <- itertest + 1
                }
                if (itertest==100) return(list(scale=NA))
            }
        }
        if(kk>0) {
            # do the refining
            tmp <- re.s(x=x,y=y,initial.beta=beta,kk=kk,conv=0,b=b,cc=cc)
            beta.rw <- tmp$beta.rw
            scale.rw <- tmp$scale.rw
            res.rw <- y - x %*% beta.rw
        } else { #kk = 0 means "no refining"
            beta.rw <- beta
            res.rw <- y - x %*% beta.rw
            scale.rw <- median(abs(res.rw))/.6745
        }
        if( i > 1 ) {
            # if this isn't the first iteration....
            scale.test <- loss.S(res.rw,s.worst,cc)
            if( scale.test < b ) {
                s.best <- scale1(res.rw,b,cc,scale.rw)
                ind <- order(best.scales)[tt]
                best.scales[ind] <- s.best
                best.betas[ind,] <- beta.rw
                s.worst <- max(best.scales)
            }
        } else { # if this is the first iteration, then this is the best beta...
            best.scales[tt]  <- scale1(res.rw,b,cc,scale.rw)
            best.betas[tt,] <- beta.rw
        }
    }
    # do the complete refining step until convergence (conv=1) starting
    # from the best subsampling candidate (possibly refined)
    super.best.scale <- 1e20
    # magic number alert
    for(i in tt:1) {
        tmp <- re.s(x=x,y=y,initial.beta=best.betas[i,],initial.scale=best.scales[i],kk=0,conv=1,b=b,cc=cc)
        if(tmp$scale.rw < super.best.scale) {
            super.best.scale <- tmp$scale.rw
            super.best.beta <- tmp$beta.rw
        }
    }
    return(list(beta=as.vector(super.best.beta), scale=super.best.scale))
}

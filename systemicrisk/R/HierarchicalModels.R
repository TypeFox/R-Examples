#########################################
## Independent Models for p and lambda ##
#########################################


#' Model for a Constant p
#'
#' This model assumes that the link existence probabilities of the
#' matrix are known.
#'
#' @param n dimension of matrix.
#' @param p existence probability of a link. Either a matrix of
#' dimension n or a single numeric value. A single numeric value leads
#' to a matrix of existence probabilities that has 0 on the diagonal.
#'
#' @return the resulting model.
#'
#' @examples
#' m <- Model.p.constant(5,0.25)
#' m$matr(m$rtheta())
#'
#' p <- matrix(c(0,0.99,0.99,0.5,0.5,0,0.01,0.01,0),nrow=3)
#' m <- Model.p.constant(5,p)
#' m$matr(m$rtheta())

#' @export
Model.p.constant <- function(n,p){
    list(dim=0,
         update= function(L,theta){
             NULL
         },
         matr=function(theta){
             if (is.matrix(p)){
                 p
             }else{
                 m <- matrix(p,nrow=n,ncol=n)
                 diag(m) <- 0
                 m
             }
         },
         rtheta=function() NULL,
         inittheta=function() NULL
         )
}


#' Model for a Random One-dimensional p
#'
#' Assumes a Beta prior on the one-dimensional link existence
#' probabilities p. This model has a one-dimensional parameter.
#'
#' @param n dimension of matrix.
#' @param shape1 first parameter of Beta prior. Default 1.
#' @param shape2 second  parameter of Beta prior. Default 1.
#'
#' @return the resulting model.
#'
#' @examples
#' m <- Model.p.BetaPrior(5)
#' m$matr(m$rtheta())
#' @export
Model.p.BetaPrior <- function(n,shape1=1,shape2=1){
    list(dim=1,
         update=function(L,theta){
             rbeta(1,shape1=shape1+sum(L>0),
                   shape2=shape2+sum(L==0)-dim(L)[1])
         },
         matr=function(theta){
             m <- matrix(theta,nrow=n,ncol=n)
             diag(m) <- 0
             m
         },
         rtheta=function() rbeta(1,shape1=shape1,shape2=shape2),
         inittheta=function() 0.5
         )
}


internal_genmatrfromvec <- function(I,vec){
    matr <- matrix(0,nrow=dim(I)[1],ncol=dim(I)[2])
    for (j in 1:length(vec)){
        matr <- ifelse(I==j,vec[j],matr)
    }
    matr
}

#' Model Using Multiple Independent Components
#'
#' Assumes a multivariate hyperparameter \eqn{\theta}{theta} with each
#' component following an independent Beta distribution. A matrix
#' indicates which component \eqn{\theta}{theta} is used for what
#' component of p.
#'
#' @param Ip matrix consisting of integers that describe which
#' component of \eqn{theta}{theta} is used for a given position in the
#' matrix. Must consist of nonnegative integers (0 encoding forced 0s
#' in the matrix), using all integers in the range.
#' @param shape1 first parameter of Beta prior on
#' \eqn{theta}{theta}. Default 1.
#' @param shape2 second parameter of Beta prior
#' \eqn{theta}{theta}. Default 1.
#'
#' @return the resulting model.
#'
#'
#' @export
Model.p.Betaprior_mult <- function(Ip, shape1=1,shape2=1){
    allids <- unique(c(Ip))
    dimp <- length(allids[allids!=0])
    if (any(1:dimp!=sort(allids[allids!=0]))) stop("Ip must consist of nonnegative integers, using all integers in the range")
    list(dim=dimp,
         update=function(L,theta){
             sum1 <- sapply(1:dimp,function(j)sum(L>0&Ip==j))
             npar <- sapply(1:dimp,function(j)sum(Ip==j))
             rbeta(dimp,shape1=shape1+sum1,
                   shape2=shape2+(npar-sum1))
         },
         matr=function(theta){
             internal_genmatrfromvec(Ip,theta)
         },
         rtheta=function(){
             rbeta(dimp,shape1=shape1,shape2=shape2)
         },
         inittheta=function() rep(0.5,dimp)
         )
}

#' Multiplicative Fitness Model for Power Law
#'
#' This model has a power law of the degree distribution with a
#' parameter {\eqn{\alpha}{alpha}} and is tuned to a desired link
#' existence probability. It is based on a fitness model.
#'
#' Every node \eqn{i} has a fitness \eqn{\theta_i}{theta_i} being an
#' independent realisation of a U[0,1] distribution.  The probability
#' of a link between a node with fitness x and a node with fitness y
#' is g(x)g(y) where g is as follows.  If \eqn{\alpha=-1}{alph=-1a}
#' then \deqn{g(x)=g0*\exp(-\log(g0)*x)}{g(x)=g0*exp(-log(g0)*x)}
#' Otherwise,
#' \deqn{g(x)=(g0^(\alpha+1)+(1-g0^(\alpha+1))*x)^(1/(\alpha+1))}{g(x)=(g0^(\alpha+1)+(1-g0^(\alpha+1))*x)^(1/(\alpha+1))}
#' where \eqn{g0}{g0} is tuned numerically to achieve the desired
#' overall mean degree.
#'
#' Updating of the model parameters in the MCMC setup is done via a
#' Metropolis-Hastings step, adding independent centered normal random
#' variables to each node fitness in \eqn{\theta}{theta}.
#'
#' @param n dimension of matrix.
#' @param alpha exponent for power law. Must be <=-1.
#' @param meandegree overall mean degree (expected degree divided by number of nodes). Must be in (0,1).
#' @param sdprop standard deviation of updated steps.
#'
#'
#' @examples
#' n <- 5
#' mf <- Model.p.Fitness.Servedio(n=n,alpha=-2.5,meandegree=0.5)
#' m <- Model.Indep.p.lambda(model.p=mf,
#'                           model.lambda=Model.lambda.GammaPrior(n,scale=1e-1))
#' x <- genL(m)
#' l <- rowSums(x$L)
#' a <- colSums(x$L)
#' res <- sample_HierarchicalModel(l,a,model=m,nsamples=10,thin=10)
#'
#'
#' @references
#' Servedio V. D. P. and Caldarelli G. and Butta P. (2004)
#' Vertex intrinsic fitness: How to produce arbitrary scale-free networks.
#' \emph{Physical Review E} 70, 056126.
#' @export
Model.p.Fitness.Servedio <- function(n, alpha, meandegree,sdprop=0.1){
   ### Mean degree of network
    getmeandegree <- function(ming){
        if (alpha==-1){
            (log(1/ming)^(-1)*(1-ming))^2
        }else{
            if (alpha==-2){
                ((alpha+1)/(1^(alpha+1)-ming^(alpha+1))*(log(1)-log(ming)))^2
            }else{
                ((alpha+1)/(alpha+2)/(1^(alpha+1)-ming^(alpha+1))*(1^(alpha+2)-ming^(alpha+2)))^2
            }
        }
    }

    ##calibrate minp to get a given mean degree
    if (meandegree<=0|meandegree>=1)
        stop("Parameter meandegree must be in (0,1).")
    ming <- uniroot(function(ming) getmeandegree(ming)-meandegree,interval=c(1e-10,1-1e-10))$root
    if (alpha==-1){
        g <- function(x) ming*exp(-log(ming)*x)
    } else
        g <- function(x) (ming^(alpha+1)+(1-ming^(alpha+1))*x)^(1/(alpha+1))
    if (g(1)>1) stop("Choice of gamma and beta not consistent - no such probability distribution exists")
    genmatr <- function(theta){
        m <- outer(g(theta),g(theta),FUN="*")
        diag(m) <- 0
        m
    }
    list(dim=n,
         update=function(L,theta){
             loglikelihood <- function(theta){
                 pakt <- genmatr(theta)
                 prob <- ifelse(L>0,log(pakt),log(1-pakt))
                 diag(prob) <- 0
                 sum(prob)
             }
             thetanew <- theta+rnorm(n,sd=sdprop)
             if (any(thetanew<=0)||any(thetanew>=1))
                 theta
             else{
                 if (exp(loglikelihood(thetanew)-loglikelihood(theta))>runif(1))
                     thetanew
                 else
                     theta
             }
         },
         matr=genmatr,
         rtheta=function() runif(n),
         inittheta=function() rep(0.5,n)
         )
}

#######################
## Models for lambda ##
#######################

#' Model for a Constant lambda
#'
#' This model assumes that the parameter lambda is known.
#'
#' @param n dimension of matrix.
#' @param lambda paramer for the size of the liabilities. Either a matrix of
#' dimension n or a single numeric value.
#'
#' @return the resulting model.
#'
#' @examples
#' m <- Model.lambda.constant(n=5,lambda=0.25)
#' m$matr(m$rtheta())
#' lambda<-matrix(c(NA,1,1,1e-4,NA,1e-4,1e4,1e4,NA),nrow=3)
#' m <- Model.lambda.constant(n=3,lambda=lambda)
#' m$matr(m$rtheta())
#' @export
Model.lambda.constant <- function(lambda,n){
    list(dim=0,
         update= function(L,theta){
             NULL
         },
         matr=function(theta){
             matrix(lambda,nrow=n,ncol=n)
         },
         rtheta=function() NULL,
         inittheta=function() NULL
         )
}

#' Model with Gamma Prior on Lambda
#'
#' Assumes that all elements of lambda are equal to a parameter
#' \eqn{\theta}{theta}, which has a Gamma prior.
#'
#' @param n dimension of matrix
#' @param shape shape paramer for prior on
#' \eqn{\theta}{theta}. Default 1.
#' @param scale scale paramer for prior on
#' \eqn{\theta}{theta}. Default 1.
#'
#' @export
Model.lambda.GammaPrior <- function(n,shape=1,scale=1){
    list(dim=1,
         update=function(L,theta){
             rgamma(1,shape=shape+sum(L>0),rate=1/scale+sum(L))
         },
         matr=function(theta){
             matrix(theta,nrow=n,ncol=n)
         },
         rtheta=function() rgamma(1,shape=shape,scale=scale),
         inittheta=function() 1
         )
}

#' Model Using Multiple Independent Components
#'
#' Assumes a multivariate hyperparameter \eqn{\theta}{theta} with each
#' component following an independent Beta distribution. A matrix
#' indicates which component \eqn{\theta}{theta} is used for what
#' component of lambda.
#'
#' @param Ilambda matrix consisting of integers that describe which
#' component of \eqn{theta}{theta} is used for a given position in the
#' matrix. Must consist of nonnegative integers using all integers in
#' the range.
#' @param shape shape paramer for prior on
#' \eqn{\theta}{theta}. Default 1.
#' @param scale scale paramer for prior on
#' \eqn{\theta}{theta}. Default 1.
#'
#' @return the resulting model.
#'
#' @export
Model.lambda.Gammaprior_mult <- function(Ilambda, shape=1,scale=1){
    allids <- unique(c(Ilambda))
    dimlambda <- length(allids[allids!=0])
    if (any(1:dimlambda!=sort(allids[allids!=0]))) stop("Ilambda must consist of nonnegative integers, using all integers in the range")
    list(dim=dimlambda,
         update=function(L,theta){
             npar <- sapply(1:dimlambda,function(j)sum(Ilambda==j&L>0))
             sumL <- sapply(1:dimlambda,function(j)sum(L*(Ilambda==j)))
             rgamma(dimlambda,shape=shape+npar,rate=1/scale+sumL)      #########TODO:  check
         },
         matr=function(theta){
             internal_genmatrfromvec(Ilambda,theta)
         },
         rtheta=function(){
             rgamma(dimlambda,shape=shape,scale=scale)
         },
         inittheta=function() rep(shape*scale,dimlambda)
         )
}


#' Combination of Independent Models for p and lambda
#'
#' @param model.p model for p.
#' @param model.lambda model for lambda.
#'
#' @examples
#' n <- 5
#' m <- Model.Indep.p.lambda(Model.p.BetaPrior(n),
#'                           Model.lambda.GammaPrior(n,scale=1e-1))
#' genL(m)
#'
#' @export
Model.Indep.p.lambda <- function(model.p,model.lambda){
    get.theta.p <- function(theta){
        if (model.p$dim>0){
            theta[1:model.p$dim]
        }else{
            NULL
        }
    }
    get.theta.lambda <- function(theta){
        if (model.lambda$dim>0){
            theta[(model.p$dim+1):length(theta)]
        }else{
            NULL
        }
    }

    list(dim=model.p$dim+model.lambda$dim,
         update=
         function(L,theta){
             c(p=model.p$update(L,get.theta.p(theta)),
               lambda=model.lambda$update(L,get.theta.lambda(theta))
               )
         },
         matr=function(theta){
             list(p=model.p$matr(get.theta.p(theta)),
                  lambda=model.lambda$matr(get.theta.lambda(theta)))
         },
         rtheta=function()c(p=model.p$rtheta(),lambda=model.lambda$rtheta()),
         inittheta=function()c(p=model.p$inittheta(),lambda=model.lambda$inittheta())
         )
}




#' Sample from Hierarchical Model with given Row and Column Sums
#'
#'
#' @param model Underlying model for p and lambda.
#' @param l observed row sum
#' @param a observerd column sum
#' @param L_fixed Matrix containing known values of L, where NA
#' signifies that an element is not known. If \code{L_fixed} equates
#' to \code{NA} (the default) then no values are assumed to be known.
#' @param burnin number of iterations for the burnin. Defaults to 5%
#' of the steps in the sampling part.
#' @param silent (default FALSE) suppress all output (including progress bars).
#' @param nsamples number of samples to return.
#' @param thin how many updates of theta to perform before outputting a sample.
#' @param matrpertheta number of matrix updates per update of theta.
#'
#' @examples
#' n <- 10
#' m <- Model.Indep.p.lambda(Model.p.BetaPrior(n),
#'                           Model.lambda.GammaPrior(n,scale=1e-1))
#' x <- genL(m)
#' l <- rowSums(x$L)
#' a <- colSums(x$L)
#' \dontrun{
#' res <- sample_HierarchicalModel(l,a,model=m)
#' }
#' # fixing one values
#' L_fixed <- matrix(NA,ncol=n,nrow=n)
#' L_fixed[1,2:5] <- x$L[1,2:5]
#' \dontrun{
#' res <- sample_HierarchicalModel(l,a,model=m,L_fixed=L_fixed,
#'                                 nsamples=1e2)
#' sapply(res$L,function(x)x[1,2:5])
#' }
#'
#' @export
sample_HierarchicalModel <- function (l, a, L_fixed=NA,
                                      model,
                                      nsamples = 10000,
                                      thin = choosethin(l=l, a=a, L_fixed=L_fixed,
                                          model=model,
                                          matrpertheta=matrpertheta,
                                          silent=silent),
                                      burnin = NA,
                                      matrpertheta=length(l)^2,
                                      silent=FALSE) {
    if (thin<1) stop("thin needs to be a positive integer")
    n <- length(l)
    res <- list()
    if (model$dim>0)
        restheta <- matrix(nrow=model$dim,ncol=nsamples)
    theta <- model$inittheta()
    u <- model$matr(theta)

    if (is.matrix(L_fixed)){ ## find starting values for fixed values.
        if (any(dim(L_fixed)!=c(length(l),length(a))))
            stop("Dimenstions of L_fixed, l and a do not match")
        if (any(L_fixed<0,na.rm=TRUE)) stop("L_fixed has negative entries")
        l_fixed <- rowSums(L_fixed,na.rm=TRUE)
        a_fixed <- colSums(L_fixed,na.rm=TRUE)
        if (any(l_fixed>l)) stop("A row sums of fixed entries is greater than desired row sum.")
        if (any(a_fixed>a)) stop("A column sum of fixed entries is greater than desired row sum.")
        L_fixed_add <- ifelse(is.na(L_fixed),0,L_fixed)
        L <- findFeasibleMatrix_targetmean(l-l_fixed,a-a_fixed,p=ifelse(is.na(L_fixed),u$p,0),targetmean=mean(genL(model)$L>0))
        samplestep <- expression({
            theta <- model$update(L+L_fixed_add,theta)
            u <- model$matr(theta)
            GibbsSteps_kcycle(L = L, p = ifelse(is.na(L_fixed),u$p,0),
                              lambda = u$lambda, it = matrpertheta)
        })
    }else{
        L <- findFeasibleMatrix_targetmean(l, a,p=u$p, targetmean=mean(genL(model)$L>0))
        samplestep <- expression({
            theta <- model$update(L,theta)
            u <- model$matr(theta)
            GibbsSteps_kcycle(L = L, p = u$p, lambda = u$lambda, it = matrpertheta)
        })
    }
    if (!silent) cat("Burn-in\n")
    if (is.na(burnin))
        burnin <- round(0.05*nsamples*thin)
    if (burnin>0){
        if (!silent)  pb <- txtProgressBar(min=0,max=burnin,style=3)
        for (i in 1:burnin){
            eval(samplestep);
            if (!silent) setTxtProgressBar(pb, i)
        }
    }

    if (!silent) cat("\nSampling\n")
    if (!silent) pb <- txtProgressBar(min=0,max=nsamples,style=3)
    for (i in 1:nsamples) {
        for (j in 1:thin){
            eval(samplestep)
        }
        if (is.matrix(L_fixed)){
            res[[i]] <- L+L_fixed_add
        }else{
            res[[i]] <- cloneMatrix(L)
        }
        if (model$dim>0)
            restheta[,i] <- theta
        if (!silent) setTxtProgressBar(pb, i)
    }
    if (!silent) cat("\n")

    if (model$dim>0)
        list(L=res,theta=restheta)
    else
        list(L=res)

}


#' Calibrate Thinning
#'
#' Attempts to automatically choose a thinning paramter to achieve an
#' overall relative effective sample size (defined as the effective
#' sample size divided by the number of samples) for all parameters in
#' the model (that do not seem to be constant). This function provides
#' no guarantees that the desired relative effective sample size
#' (rESS) will actually be achieved - it is best treated as a rough
#' guide for this.
#'
#' The approach used involves a pilot run of the sampler, followed by
#' a computation of the acf (autocorrelation function) for each
#' component. The acf is used only up to (and excluding) the point
#' used where it becomes negative for the first time. This part of the
#' acf is then used to approximate the rESS and to determine the
#' amount of thinning needed. The reported result is the thinning
#' needed to achieve the rESS for all components (the matrix as well
#' as the parameter theta). The initial pilot run may not be
#' sufficient and further pilot runs may have to be started.
#'
#' @inheritParams sample_HierarchicalModel
#' @param relESStarget Target for the relative effective sample size,
#' must be in (0,1). Default 0.3.
#' @param maxthin Upper bound on thinning to consider. Default 10000.
#'
#' @examples
#' n <- 10
#' m <- Model.Indep.p.lambda(Model.p.BetaPrior(n),
#'                           Model.lambda.GammaPrior(n,scale=1e-1))
#' x <- genL(m)
#' l <- rowSums(x$L)
#' a <- colSums(x$L)
#' \dontrun{
#' choosethin(l,a,model=m)
#' choosethin(l,a,model=m,relESStarget=0.7)
#' }
#'@export
choosethin <- function(l, a, L_fixed=NA,
                         model,
                         relESStarget=0.3, burnin = 100,
                         matrpertheta=length(l)^2,
                         silent=FALSE,
                         maxthin=10000) {
    if (relESStarget<=0 || relESStarget>=1) stop("relESStarget needs to be in (0,1).")
    if (!silent)
        cat("Determining thinning needed to achieve relative ESS of ",relESStarget,".\n",sep="")
    ## adjust target to build in a safety margin
    relESStarget <- min(0.1,(1-relESStarget)/2)+relESStarget
    if (!silent)
        cat("To include a safety margin, aiming for relative ESS of ",relESStarget,".\n",sep="")

    nsamples <- 1e3
    thin <- 1
    while(thin<maxthin){
        if (!silent) cat("Pilot run with thin=",thin," started.\n",sep="")
        res <- sample_HierarchicalModel(l=l,a=a,L_fixed=L_fixed,model=model,burnin=burnin,matrpertheta=matrpertheta,silent=silent,thin=thin,nsamples=nsamples)
        allvariables <- sapply(res$L,function(x) c(x))
        allvariables <- rbind(allvariables,res$theta)
        allacf <- sapply(1:(dim(allvariables)[1]), function(i)acf(allvariables[i,],plot=FALSE,lag.max=100)$acf)
        thinfac <- apply(allacf,2,function(x) {
            if (all(is.na(x))) return(0);
            if (all(x>0)) return(Inf)
            negind <- which(x<0) ##find first negative autocorrelation
            if (length(negind)==0)
                Inf
            else{
                x <- x[1:(min(negind)-1)]
                if (length(x)==1) 1
                else {
                    rESS <- sapply(1:(length(x)),function(i){
                        if (i==length(x)) 1
                        else 1/(1+2*sum(x[1+(1:floor((length(x)-1)/i))*i]))})
                    if (all(rESS<relESStarget))
                        Inf
                    else
                        min(which(rESS>=relESStarget))
                }
            }
        })
        thinres <- max(thinfac*thin)
        if (thinres<=maxthin) {
            if (!silent) cat("Recommend thin=",thinres,".\n",sep="")
            return(thinres)
        } else {
            thin <- thin*10;
            if (!silent) {
                cat("Need to increase thinning in pilot sample.\n")
                names <- c(outer(1:length(l),1:length(a),FUN=function(x,y) paste("L[",x,",",y,"]",sep="")))
                if (!is.null(res$theta))
                    names <- c(names,paste("theta[",1:(dim(res$theta)[1]),"]",sep=""))
                cat("Variables for which determining required thinning failed are:\n",
                    names[thinfac*thin>maxthin],"\n")
                 cat("Using now thin=",thin,".\n",sep="")
            }
        }
    }
    if (!silent) warning("Upper bound for thinnig reached.")
    return(maxthin)
}

#' Generate Liabilities Matrix from Prior
#'
#' Generates a libabilities matrix using a the prior distribution from
#' a given model for p and lambda.
#'
#' @param model a model for p and lambda.
#' @return A list consisting of a liabilities matrix and the parameter
#' vector theta.
#'
#' @examples
#' n <- 5
#' m <- Model.Indep.p.lambda(Model.p.BetaPrior(n),
#'                           Model.lambda.GammaPrior(n,scale=1e-1))
#' genL(m)
#'
#' @export
genL <- function(model){
    theta <- model$rtheta()
    u <- model$matr(theta)
    n <- dim(u$p)[1]
    res <- matrix((runif(n*n)<=u$p)*rexp(n*n,rate=u$lambda),nrow=n)
    list(L=res,theta=theta)
}


getfixed <- function(l,a,printresults=TRUE){

    ##Ignore columns and rows with sum 0
    NArow <- l==0
    if (sum(NArow)>0){
        cat("Row sum ",paste(which(NArow),collapse=" "),"are zero.\n")
    }
    NAcol <- a==0
    if (sum(NArow)>0){
        cat("Column sum ",paste(which(NArow),collapse=" "),"are zero.\n")
    }

    ## this is not a full check of the condition, but it helps a bit
    ## check for individual matches
    wl <- l!=0&l%in%a
    if (length(wl)>0){
        for (i in which(wl)){
            wa <- which(a==l[i])
            if (length(wa)==1){
                cat("Sum of Row ",i," = Sum of Column ",wa,". ",
                    "\n",sep="")
                NArow[i] <- TRUE
                NAcol[wa] <- TRUE
            }else{
                cat("Sum of Row ",i," matches sums of columns ",wa,". ",
                    "Rows and Columns NOT ignored.\n",sep="")
            }
        }
    }
    ## check for individual matches against sum of two elements
    { ##one row element
        s <- outer(a,a,"+")
        diag(s)=0;
        s[NAcol,] <- NA
        s[,NAcol] <- NA
        for (i in 1:length(l)){
            if (!NArow[i]&&l[i]>0){
                w <- l[i]==s
                if (any(w,na.rm=TRUE)){
                    wa <- which(rowSums(w,na.rm=TRUE)>0)
                    if (length(wa)==2){
                        cat("Row sum ",i," equal to sum of columns ",paste(wa,collapse=" "),". ",
                            "\n",sep="")
                        NArow[i] <- TRUE
                        NAcol[wa] <- TRUE
                    }else{
                        cat("Row sum ",i,"equal to some subsets of  the columns ",paste(wa,collapse=" "),". ",
                            "Matrix not necessarily deterministics.\n",sep="")
                    }


                }
            }
        }
    }
    { ##one column element
        s <- outer(l,l,"+")
        diag(s)=0;
        s[NArow,] <- NA
        s[,NArow] <- NA
        for (i in 1:length(l)){
            if (!NAcol[i]&&a[i]>0){
                w <- a[i]==s
                if (any(w,na.rm=TRUE)){
                    wl <- which(rowSums(w,na.rm=TRUE)>0)
                    if (length(wl)==2){
                        cat("Column sum ",i," equal to sum of rows ",paste(wl,collapse=" "),".\n",sep="")
                        NAcol[i] <- TRUE
                        NArow[wl] <- TRUE
                    }else{
                        cat("Column sum  ",i,"equal to two of the row sums: ",paste(wl,collapse=" "),".",
                            "Matrix not necessarily deterministics.\n",sep="")
                    }


                }
            }
        }
    }
    list(row=NArow,col=NAcol)
}

#' Outputs Effective Sample Size Diagonis for MCMC run
#'
#' Computes the Effective Sample Size using the method
#' \code{effectiveSize} in of the package \code{coda}.
#'
#' Currently only works with L where the diagonal is 0. The function
#' ignores the diagonal and tries to determine from the row and column
#' sums which parts of the matrix are 0.
#'
#' @param res output from \code{\link{sample_HierarchicalModel}}.
#'
#' @export
diagnose <- function(res){
    if (!all(sapply(res$L,diag)==0))
        stop("Model does not have an enforced 0 on the diagonal. This function is not designed for this.")

    fixed <- getfixed(l=rowSums(res$L[[1]]),a=colSums(res$L[[1]]))
    n <- dim(res$L[[1]])[1]
    if (dim(res$L[[1]])[2]!=n) stop("Works only for quadratic matrices")

    NArow <- fixed$row
    NAcol <- fixed$col


    rmNA <- function(L) {
        diag(L) <- NA;
        L[NArow,] <- NA;
        L[,NAcol] <- NA;
        na.omit(c(L))
    }
    if (length(rmNA(res$L[[1]]))==0){
        stop("Entire matrix deterministic. No remaining elements to sample over! \n")
    }
    allL <- sapply(res$L, rmNA)
    cat("Analysis does not consider ",prod(dim(res$L[[1]]))-dim(allL)[1]," entries of matrix \n",
        "that are deterministic (diagonal elements, row/column sum=0 or forced result).\n",sep="")
    ## check how many matrix entries have moved
    notmoved <- rowSums((allL[,1]!=allL))==0
    if (sum(notmoved)==0)
        cat("All remaining elements of the liabilities matrix have moved during sample run.\n")
    else{
        cat("Number of matrix entries not moving: ",sum(notmoved),"\n")
        cat("Of those always equal to 0: ",sum(allL[notmoved,1]==0),"\n")
        rownotmoved <- rmNA(matrix(1:n,nrow=n,ncol=n))[notmoved]
        colnotmoved <- rmNA(t(matrix(1:n,nrow=n,ncol=n)))[notmoved]
        o <- order(rownotmoved,colnotmoved)
        rownotmoved <- rownotmoved[o]
        colnotmoved <- colnotmoved[o]
        lastr <- 0
        for (i in 1:min(100,length(rownotmoved))){
            if (rownotmoved[i]!=lastr){
                lastr <- rownotmoved[i]
                cat("\n Row ",lastr,"; Columns:")
            }
            cat(" ",colnotmoved[i])
        }
        cat("\n")
        if (length(rownotmoved)>100) cat("Output truncated to 100 not-moving values\n")
    }
    cat("ESS in matrix:\n")
    print(summary(coda::effectiveSize(t(allL))))
    if ("theta" %in% names(res)){
        cat("ESS in theta:\n")
        print(summary(coda::effectiveSize(t(res$theta))))
    }else
        cat("No hyperparameter theta\n")
}

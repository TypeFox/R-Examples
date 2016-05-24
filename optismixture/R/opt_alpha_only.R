test1 <- function(){
    get.initial.alpha <- function(eps, J){
        if(length(eps) != J)
            stop("wrong eps input")
        k <- (1 - sum(eps))/(J + sqrt(J))
        alpha <- eps + k
        return(alpha)
    }
    J <- 2
    y <- 1
    z <- matrix(c(2,3),1,2)
    x <- matrix(0, 1,1)
    eps <- rep(0.01/J, J)
    reltol <- 10^-7
    relerr <- 10^-7
    e <- penoptpersp.alpha.only(y, z, eps, reltol ,relerr, a0 = get.initial.alpha(rep(0.1/J, J), J))
}

#' penalized optimization of the constrained linearized perspective function
#' @export
#' @export
#' @param y length \eqn{n} vector
#' @param z \eqn{n \times J} matrix
#' @param a0 length \eqn{J} vector
#' @param eps length \eqn{J} vector, default to be rep(0.1/J, J)
#' @param reltol  relative tolerence for Newton step, between 0 to 1, default to be \eqn{10^{-3}}. For each inner loop, we optimize \eqn{f_0 + \rho \times \mathrm{pen} } for a fixed \eqn{\rho}, we stop when the Newton decrement \eqn{f(x) - inf_y \hat{f}(y) \leq f(x)* \mathrm{reltol}}, where \eqn{\hat{f}} is the second-order approximation of \eqn{f} at \eqn{x}
#' @param relerr relerr stop when within (1+\emph{relerr}) of minimum variance, default to be \eqn{10^{-3}}, between 0 to 1.
#' @param rho0 initial value for \eqn{\rho}, default to be 1
#' @param maxin maximum number of inner iterations
#' @param maxout maximum number of outer iterations
#' @return a list of \describe{
#'   \item{y}{input y}
#'   \item{z}{input z}
#'   \item{alpha}{optimized alpha}
#'   \item{rho}{value of rho}
#'   \item{f}{value of the objective function}
#'   \item{rhopen}{value of rho*pen when returned}
#'   \item{outer}{number of outer loops}
#'   \item{relerr}{relative error}
#'   \item{alphasum}{sum of optimized alpha}
#' }
#' @details To minimize \eqn{\sum_i  \frac{y_i^2}{z_i^T\alpha}} over \eqn{\alpha}
#'   subject to \eqn{\alpha_j > \epsilon_j} for \eqn{j = 1, \cdots, J} and \eqn{\sum_{j=1}^J \alpha_j < 1},
#'
#' Instead we minimize \eqn{ \sum_i  \frac{y_i^2}{z_i^T\alpha} + \rho \times \mathrm{pen}} for a decreasing sequence of \eqn{\rho}
#'
#' where \eqn{ \mathrm{pen} = -( \sum_{j = 1}^J( \log(\alpha_j-\epsilon_j) ) + \log(1-\sum_{j = 1}^J \alpha_j) )}
#'
#' starting values are \eqn{\alpha = a0} and can be missing.
#'
#' The optimization stops when within (1+\emph{relerr}) of minimum variance.

penoptpersp.alpha.only = function( y,z, a0, eps = NULL, reltol = NULL, relerr = NULL,  rho0 = NULL, maxin= NULL, maxout = NULL){
    J = ncol(z)
    if(sum(y^2) == 0){
        print("All y are zero. Cannot optimize for alpha")
        alpha <- a0
        return(list(y=y, z=z, alpha=alpha,f=0, alphasum = sum(alpha)))
    }
                                        # Verify that inputs are ok and feasible

    if(is.null(eps))
        eps <- rep(0.1/J, J)
    else if( length(eps)==1 )
        eps = rep(eps,J)
    else if(length(eps) != J)
        stop("Wrong dimension of eps")
    if( min(eps)<0 )stop("Negative epsilon")
    if(is.null(reltol)) reltol <- 10^-3
    else if(reltol < 0 | reltol > 1) reltol <- 10^-3
    if(is.null(relerr)) relerr <- 10^-3
    else if(relerr > 1 | relerr < 0) relerr <- 10^-3
    if(is.null(rho0)) rho0 <- 1
    if(is.null(maxin))    maxin <- 20
    if(is.null(maxout))    maxout <- 30
    if( min(eps)<0 )stop("Negative epsilon")
    if( !is.null(a0) ){
        if( any(a0<eps) )stop("Infeasible start point, any(a0<eps)")
        if( sum(a0)>1 )stop("Infeasible start point, sum(a0)>1")
    }
    if(any(apply(z, 1, function(x) sum(x^2)) == 0)) stop("Some rows of z is zero")


    delta = (1 -sum(eps))/(J+1)
    if( delta < 0 )stop("No feasible alpha")
    ## if( K==0 )warning("No regression portion") DLT

    ## svdrange = range(svd(x)$d) DLT
    ## print(svdrange)


                                        # main workhorse function
    ## fgH = function(alpha,beta,rho,do=3){ DLT
    fgH = function(alpha,rho,do=3){
                                        # do = 3    function, gradient, Hessian
                                        # do = 2    function, gradient
                                        # do = 1    function
                                        # do = 0    feasibility

                                        # The penalty
        if( any(alpha<eps) )
            pen = Inf
        else if( sum(alpha)>1  )
            pen = Inf
        else
            pen = -( sum( log(alpha-eps) ) + log(1-sum(alpha) ) )

        if( do <= 0 )return( pen < Inf )

                                        # The objective function
        ## res = y - x%*%beta  # n-vector of residuals DLT
        res = y  # n-vector of residuals
        zal = z%*%alpha     # n-vector of mixture probabilities
        ebz = res/zal

        f0   = sum( res^2/zal )
        f   = f0 + rho * pen

        if( is.nan(f) ){
            print("Hit a NaN")
            print(pen)
            print(alpha)
            print(sum(alpha))
        }

        if( do <=1 )return(list(f=f, f0 = f0))

                                        # The gradient

        gpen  = -(1/(alpha-eps) - 1/(1-sum(alpha)))
        ## gbeta = -2 * t(x) %*% ebz DLT
        galph = -t(z) %*% ebz^2
        ## g0 <- c(gbeta, galph) DLT
        ## gphi <- c(rep(0, K), rho*gpen)
        ## g     = c( gbeta, galph + rho*gpen )
        g0 <- c(galph)
        gphi <- rho*gpen
        g     = c(galph + rho*gpen)

        if( do <=2 )return(list(f=f,g=g,f0 = f0, g0 = g0, gphi = gphi))

                                        # The Hessian
        ## DLT
        ## xbyrootq = x
        ## for( j in 1:ncol(x) )
        ##     xbyrootq[,j] = xbyrootq[,j]/sqrt(zal)
        ## Hbb = 2*t(xbyrootq) %*% xbyrootq
        ## DLT

        zebyqrootq = z
        for( j in 1:ncol(z) )
            zebyqrootq[,j] = zebyqrootq[,j] * res /zal^(3/2)
        Haa = 2 * t(zebyqrootq) %*% zebyqrootq

        ## Hba = 2 * t(xbyrootq) %*% zebyqrootq DLT


        Hpen = diag( (alpha-eps)^-2 ) + 1/(1-sum(alpha))^2

        ## H = rbind( cbind(    Hbb, Hba ),
        ##     cbind( t(Hba), Haa + rho*Hpen ) ) DLT
        H = Haa + rho*Hpen

        return(list(f=f,g=g,H=H, f0 = f0, g0 = g0, gphi = gphi))
    }

    ## linesearch = function(alpha,beta,rho,direction,oldf,oldg,LSalpha=.15,LSbeta=0.45){
    linesearch = function(alpha,rho,direction,oldf,oldg,LSalpha=.15,LSbeta=0.45){
                                        # backtracking line search, Boyd and Vandenberghe p 464, Alg 9.2
                                        # line search parameters:
                                        #    LSbeta  typically between .1 (crude search) and .8
                                        #    LSalpha typically between .01 and 0.30

                                        #print("in linesearch")
                                        #print(oldf)
                                        #print(oldg)
                                        #print(direction)
        tval = 1
        while(1){
            ## newbeta  = beta  + tval * direction[1:K] DLT
            ## newalpha = alpha + tval * direction[K+(1:J)] DLT
            newalpha = alpha + tval * direction
            ## newfg = fgH(newalpha,newbeta,rho,do=1) DLT
            newfg = fgH(newalpha,rho,do=1)
                                        #  print("newfg");print(newfg)
            ## if( fgH(newalpha,newbeta,rho,do=0) && DLT
            if( fgH(newalpha,rho,do=0) &&
               (newfg$f <= oldf + LSalpha * tval * direction%*%oldg) )
                break
            tval = tval * LSbeta
                                        #  print(tval)
                                        #  print("newalpha");print(newalpha)
                                        #  print("newbeta");print(newbeta)
        }
        tval
    }

    svdsolve = function(A,b,epsrel=1e-9,epsabs=1e-100){
                                        # Solve Av = b via SVD
        A.svd <- svd(A)
        d <- A.svd$d
        index1 <- which(d>=(epsrel*d[1] +epsabs))
        index2 <- which(d<(epsrel*d[1] +epsabs))
        d.trun.inv <- d
        d.trun.inv[index1] <- 1/d[index1]
        d.trun.inv[index2] <- 0
        A.inv <- A.svd$v%*%diag(d.trun.inv)%*%t(A.svd$u)
        return(A.inv%*%b)
    }
    ## DLT
    ## preconditionsolve = function(A, b, reltol=.Machine$double.eps){
    ##     J <- (dim(A)[1] + 1)/2
    ##     p <- sqrt(median(abs(A[1:(J-1), 1:(J-1)]))/median(abs(A[J:(2*J-1), J:(2*J - 1)])))
    ##     P.vec <- c(rep(1, J - 1), rep(p, J))
    ##     P <- diag(P.vec)
    ##     A.pc <- P%*%A%*%P
    ##     if(kappa(A.pc) < kappa(A)){
    ##         x <- try(as.vector(P %*% svdsolve(A.pc, P%*%b)), silent = TRUE)
    ##     }else{
    ##         x <- try(as.vector(svdsolve(A, b)), silent = TRUE)
    ##     }
    ##     return(x)
    ## }

    ## dampednewton = function(alpha,beta,rho,reltol){ DLT
    dampednewton = function(alpha,rho,reltol){
                                        # From Boyd and Vandenberghe p 487, Alg 9.5
        done = FALSE
        inct <- 0
        oldf <- -Inf
        while( 1 ){
            ## vals  = fgH(alpha,beta,rho,do=3) DLT
            vals  = fgH(alpha,rho,do=3)
            if(identical(oldf,vals$f)){
                solstatus <- "exact"
                break
            }
            oldf <- vals$f
                                        #  print(vals)
            ## print(svd(vals$H)$d)
                                        #  newtonstep = -solve(vals$H,vals$g,reltol=1e-50) # aggressively small reltol here
            newtonstep = try(-svdsolve(vals$H,vals$g), silent = TRUE)
                                        #  print("newtonstep"); print(newtonstep)
            if(class(newtonstep) == "try-error"){
                if(isTRUE(all.equal(vals$g, rep(0, length(vals$g))))){
                    solstatus <- "exact"
                    break
                }else{
                    print("returning inexact solution")
                    solstatus <- "inexact"
                    break
                    ## H <- vals$H
                    ## save(x, y, z, H, alpha, beta, rho, eps, file = "../errorworkspace/largeconditionerrorH.RData")
                    ## print("largeconditionerrorH.RData saved")
                    ## stop("condition number still too large after preconditioning")
                }
            }
            newtonstep <- drop(newtonstep)
            decrement  = - (vals$g %*% newtonstep)[1,1]/2
            ## print(paste("vals$f", vals$f))
            ## print(paste("vals$f0", vals$f0))
            ## print(paste("decrement", decrement))
            reldecrement = decrement/vals$f
            if( reldecrement <= reltol ){
                solstatus <- "exact"
                break
            }
            ## tval  = linesearch(alpha,beta,rho,newtonstep,vals$f,vals$g) DLT
            tval  = linesearch(alpha,rho,newtonstep,vals$f,vals$g)
                                        #  print("tval from linesearch");print(tval)
            ## beta  = beta  + tval * newtonstep[1:K]
            ## alpha = alpha + tval * newtonstep[K+(1:J)] DLT
            alpha = alpha + tval * newtonstep
                                        #  print("alpha");print(alpha)
                                        #  print("beta");print(beta)
            inct = inct + 1
            if(inct >= maxin){
                print("Reaching maximum inner iterations")
                solstatus <- "semiexact"
                print(paste("decrement", decrement))
                print(paste("vals$f", vals$f))
                break
            }
        }
        ## print(paste("inct", inct))
        ## print(paste("sum(alpha)", sum(alpha)))
        ## list(alpha=alpha,beta=beta, solstatus = solstatus, inner = inner) DLT
        list(alpha=alpha,solstatus = solstatus, inct = inct)
    }

    if( is.null(a0))
        alpha = eps + delta
    else
        alpha = a0
    ## if( missing(b0) ) DLT
    ##     beta = rep(0,K)
    ## else
    ##     beta = b0

                                        # Use Boyd and Vandenberghe p 569, Alg 11.1, Barrier method


    ## print("starting alpha")
    ## print(alpha)
    ## initfgH <- fgH(alpha, beta, rho=1, do=2) DLT
    initfgH <- fgH(alpha, rho=1, do=2)
    rho <- rho0
    mu  = 10    # penalty increase factor
    outer.ct <- 0
    inct.vec <- c()
    while(1){
        ## thedn  = dampednewton(alpha,beta,rho,reltol=reltol) DLT
        thedn  = dampednewton(alpha,rho,reltol=reltol)
        inct.vec <- c(inct.vec, thedn$inct)
        if(thedn$solstatus == "inexact"){
            warning("cannot achieve required relerr")
            ## thevar = fgH(alpha,beta,rho=0,do=1)$f DLT
            thevar = fgH(alpha,rho=0,do=1)$f
            dual.opt <- thevar - (J+1)*rho*mu
            if(dual.opt > 0){
                relerr <- (J+1)*rho*mu/dual.opt
                print("relative error of the last available exact solution")
                print(relerr)
            }else{
                abserr <- (J+1)*rho*mu
                print("absolute error of the last available exact solution")
                print(abserr)
                print("thevar at last available exact solution, i.e. the value of the unpenalized objective function at the solution")
                print(thevar)
            }
            break
        }

        alpha  = thedn$alpha
        ## beta   = thedn$beta DLT
        ## thevar = fgH(alpha,beta,rho=0,do=1)$f DLT
        thevar = fgH(alpha,rho=0,do=1)$f
                                        #  print(thevar)
        dual.opt <- thevar - (J+1)*rho
        if( (J+1)*rho < dual.opt*relerr ){ #Boyd p 242 relerr
            ## print("relative error reached relerr")
            relerr <- (J+1)*rho/thevar
            ## print(relerr)
            break
        }
        ## f      = fgH(alpha,beta,rho=0,do=1)$f DLT
        ## rhopen = fgH(alpha,beta,rho=rho,do=1)$f-f DLT
        f      = fgH(alpha,rho=0,do=1)$f
        rhopen = fgH(alpha,rho=rho,do=1)$f-f
        rho = rho/mu
        outer.ct <- outer.ct + 1
        if(outer.ct >= maxout){
            print("Reaching maximum outer iterations")
            break
        }
    }
    ## print("alpha")
    ## print(alpha)
    ## print("sum(alpha)")
    ## print(sum(alpha))
    ## print("number of outer iterations")
    ## print(outer.ct)
    ## print("number of inner iterations")
    ## print(inct.vec)
    ## print("total inner iterations")
    ## print(inner)

    ## f      = fgH(alpha,beta,rho=0,do=1)$f DLT
    ## rhopen = fgH(alpha,beta,rho=rho,do=1)$f-f DLT
    f      = fgH(alpha,rho=0,do=1)$f
    rhopen = fgH(alpha,rho=rho,do=1)$f-f
    ## print(paste("outer loops", log(rho0/rho,mu)))
    ## list(x=x, y=y, z=z, alpha=alpha,beta=beta,rho=rho,f=f,rhopen=rhopen,outer=log(rho0/rho,mu), relerr = relerr, alphasum = sum(alpha)) DLT
    list(y=y, z=z, alpha=alpha,rho=rho,f=f,rhopen=rhopen,outer.count =outer.ct, inner.counts = inct.vec, relerr = relerr, alphasum = sum(alpha))
                                        #list(allalpha=allalpha,allbeta=allbeta,fgH=fgH(alpha,beta,rho,do=3))
}

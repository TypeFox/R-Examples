#' Internal function. Calculate the initial alpha vector for the optimization of \emph{alpha} with a lower bound constraint
#' @param eps a vector of lower bound values for vector \emph{alpha}
#' @param J number of mixture components
#' @return The initial alpha vector for optimization
get.initial.alpha <- function(eps, J){
    if(length(eps) != J)
        stop("input eps is not of length J")
    k <- (1 - sum(eps))/(J + 1)
    alpha <- eps + k
    return(alpha)
}

#' Test the compatibility of user defined functions \emph{fname, rpname, rqname, dpname, dqname} with \emph{mixture.param}
#' @export
#' @inheritParams do.mixture.sample
#' @return Stop with error or print the success message.
compatible.test <- function(fname, rpname, dpname, rqname, dqname, mixture.param){
    a <- try({
        for(j in 1:(mixture.param$J - 1)){
            y1 <- eval(parse(text = rqname))(1, j, mixture.param)
            y2 <- eval(parse(text = rqname))(2, j, mixture.param)
            a1 <- eval(parse(text = fname))(y1, j, mixture.param)
            a2 <- eval(parse(text = fname))(y2, j, mixture.param)
            d1 <- eval(parse(text = dqname))(y1,j, mixture.param)
            d2 <- eval(parse(text = dqname))(y2,j, mixture.param)
        }
        x1 <- eval(parse(text = rpname))(1, mixture.param)
        x2 <- eval(parse(text = rpname))(2, mixture.param)
        e1 <- eval(parse(text = fname))(x1, j, mixture.param)
        e2 <- eval(parse(text = fname))(x2, j, mixture.param)
        b1 <- eval(parse(text = dpname))(x1,mixture.param)
        b2 <- eval(parse(text = dpname))(x2, mixture.param)
    })
    if(class(a) == "try-error")
        stop("Compatible test of functions and mixture.param did not pass. Please check definition of function or mixture.param.")
    if(length(a1) != 1 | length(a2) != 2
       | length(e1) != 1 |length(e2) !=2){
        stop("Returned value of f has wrong length")
    }
    if(length(d1) != 1 | length(d2) !=2){
        stop("Returned value of dq has wrong length")
    }
    if(length(b1) != 1 | length(b2) != 2){
        stop("Returned value of dp has wrong length")
    }
    print("Definition of f, rp, dp, rq, dq passed test")
}



#' Internal function. With stratified samples, calculate the variance of the estimate from importance sampling without control variates
#' @param Y vector of stratified samples of length \eqn{n}. i.e. \eqn{Y_1 = Y[1:nvec[1]]} are sampled from \eqn{q_1}, \eqn{Y_i = Y[(nvec[i-1]+1):nvec[i]]} are sample from \eqn{q_i}.
#' @param nvec the vector of number of samples from each mixture component. It sums up to \eqn{n}.
#' @return the variance estimate of \eqn{\hat{\mu} = 1/n \sum_{i=1}^n Y[i]}
#' @details Suppose we sample Y from a mixture \eqn{q_{\alpha} = \alpha_1*q_1 + ... + \alpha_J*q_J}. To estimate \eqn{\mathrm{mean}(Y)}, fixing the number of samples from each mixture component and getting a stratified sample would reduce the variance of the estimate. The formula for \eqn{\mathrm{Var}(\hat{\mu})} with stratified samples is \deqn{\mathrm{Var}(\hat{\mu}) = 1/n \times \sum_{j=1}^J \alpha_j \mathrm{Var}(Y_j)}
get.var <- function(Y, nvec){
    if(length(Y) != sum(nvec)){
        print(length(Y) - sum(nvec))
        print(paste("nvec", nvec))
        stop("length(Y) != sum(nvec)")
    }
    if(any(nvec == 1)){
        warning(paste0("Cannot estimate sd due to small sample size(only one sample) in some of the mixture components in batch 2. Consider providing a batch.size[2] at least twice as large as the current batch.size[2]"))
    }
    nvec.cum <- c(0, cumsum(nvec))
    variance <- 0
    for(i in 1:length(nvec)){
        variance <- variance + nvec[i]/sum(nvec)*var(Y[(nvec.cum[i]+1):nvec.cum[i+1]])
    }
    variance <- variance/length(Y)
    return(variance)
}

#' Internal function. Get the row index in the stacked sample matrices for the \eqn{b^{th}} batch
#' @param batch.size  vector of batch sizes
#' @param b  the index of the batch of interest
#' @return the row index or sample matrices of batch b for all b's in b.vec
#' @details In the computation, samples from each batches are stacked be rows to form x.mat,
#' as well as fx.mat, px.mat, qx.mat. To extract the samples from batch b and their corresponding
#' fx, px, qx values, the indices for batch b will be needed. They can be obtained by calling
#' get.index.b(batch.size, b)
get.index.b <- function(batch.size, b){
    index.b <- (sum(batch.size[0:(b-1)]) + 1):(sum(batch.size[0:b]))
    return(index.b)
}


#' Internal function. convert mixture proportions to mixture sample size with a fixed total sample size
#' @param n total number of samples
#' @param alpha the vector of mixture proportions
#' @return the vector of sample sizes for each mixture component
alpha2N <- function(n, alpha){
    alpha <- alpha/sum(alpha)
    alphacum <- cumsum(alpha)
    J <- length(alpha)
    nvec <- round(n * alphacum)
    nvec <- c(nvec[1], nvec[2:J] - nvec[1:(J-1)])
    return(nvec)
}


#' Internal function. sample from the mixture distribution \eqn{q_{\alpha}}
#' @param seed seed for sampling
#' @param b batch index for the samples
#' @param n total number of samples
#' @param J number of mixture components, including the defensive one
#' @param mixture.param mixture.param = list(p, J, ...), where \eqn{p} is the dimension of the sample, and \eqn{J} is the number of mixture components, including the defensive one. mixture.param should be compatible with user defined functions \emph{f(n, j, mixture.param), rp(n, mixture.param), rq(n, j, mixture.param), dp(xmat, mixture.param), dq(xmat, j, mixture.param)}
#' @param alpha vector of mixture proportions
#' @param fname name of user defined function \emph{fname(xmat, j, mixture.param)}. \emph{xmat} is an \eqn{n \times p} matrix of \eqn{n} samples with \eqn{p} dimensions. \emph{fname} returns a vector of function values for each row in \emph{xmat}. \emph{fname} is defined for \eqn{j = 1, \cdots, J}. \eqn{j = 1, \cdots, J - 1} corresponds to different proposal mixture components, and \eqn{j = J} corresponds to the defensive mixture component.
#' @param rpname name of user definded function \emph{rpname(n, mixture.param)}. It generates \eqn{n} random samples from target distribution \emph{pname}. Parameters can be specified in \emph{mixture.param}. \emph{rpname} returns an \eqn{n \times p} matrix.
#' @param rqname name of user defined function \emph{rqname(n, j, mixture.param)}. It generate \eqn{n} random samples from the \eqn{j^{th}} mixture component of proposal mixture distribution. \emph{rqname} returns an \eqn{n \times p} matrix. \emph{rqname} is defined for \eqn{j = 1, \cdots, J - 1}.
#' @param dpname name of user defined function \emph{dpname(xmat, mixture.param)}. It returns the density of \emph{xmat} from the target distribution \eqn{pname} as a vector of length \emph{nrow(xmat)}. Note that only the ratio between \emph{dpname} and \emph{dqname} matters. So \emph{dpname} and \emph{dqname} can return values of \eqn{C \times}\emph{dpname} and \eqn{C \times}\emph{dqname} respectively.
#' @param dqname name of user defined function \emph{dqname(xmat, j, mixture.param)}. It returns the density of \emph{xmat} from the proposal distribution \eqn{q_j} as a vector of length \emph{nrow(xmat)}. \emph{dqname} is defined for \eqn{j = 1, \cdots, J - 1}. Note that only the ratio between \emph{dpname} and \emph{dqname} matters. So \emph{dpname} and \emph{dqname} can  return values of \eqn{C \times}\emph{dpname} and \eqn{C \times}\emph{dqname} respectively.
#' @return a list of \describe{
#'   \item{x}{the matrix(size \eqn{n \times p}) of samples from the mixture distribution \eqn{q_{\alpha}}}
#'   \item{fx}{the vector of \emph{fname} evaluated at sample matrix \emph{x}}
#'   \item{qx}{the matrix of densities under each mixture component, i.e. \emph{qx[,j]} is the density under \eqn{q_j}. The defensive component is  the \eqn{J^{th}} column}
#'   \item{px}{the vector of densities under proposal distribution \emph{pname}}
#'   \item{alpha}{the rounded \emph{alpha} used for sampling, i.e. \eqn{alpha2N(n, alpha)/n}}
#' }
do.mixture.sample <- function(seed, b, n, J, mixture.param, alpha, fname, rpname, rqname, dpname, dqname){
    set.seed(seed)
    rpname <- as.character(rpname)
    rqname <- as.character(rqname)
    dpname <- as.character(dpname)
    dqname <- as.character(dqname)
    fname <- as.character(fname)

    nvec <- alpha2N(n, alpha)
    if(any(nvec ==0))
        stop(paste0("Some mixture components have 0 samples. ", "batch.size[", b, "]", " should be at least ", ceiling(1/min(alpha))))
    alpha <- nvec/sum(nvec)
    xmat <- matrix(NA, sum(nvec), mixture.param$p)
    fx <- matrix(NA, sum(nvec), 1)
    for(j in 1:J){
        index.j <- get.index.b(nvec, j)
        if(length(index.j) > 0){
            if(j < J){
                xmat[index.j, ] <- eval(parse(text = rqname))(nvec[j], j, mixture.param)
            }else if(j == J)
                 xmat[index.j, ] <- eval(parse(text = rpname))(nvec[j], mixture.param)
        }
        fx[index.j,] <- eval(parse(text = fname))(matrix(xmat[index.j, ], length(index.j), mixture.param$p), j, mixture.param)
    }
    qx <- matrix(NA, n, J)
    for(j in 1:(J-1)){
        qx[,j] <- eval(parse(text = dqname))(xmat, j, mixture.param)
    }
    px <- eval(parse(text = dpname))(xmat, mixture.param)
    qx[,J] <- px
    return(list(x = xmat, fx = fx, qx = qx, px = px, alpha = alpha))
}

#' Two stage estimation, a pilot estimate of mixing alpha and a following importance sampling, with or without control variates
#' @export
#' @param batch.size length two vector of batch sizes
#' @param eps the lower bound for optimizing \emph{alpha}. if \emph{eps} is of length 1, it is expanded to \emph{rep(eps, J)}, default to be \emph{rep(0.1/J, J)}
#' @param cv TRUE indicates optimizing \emph{alpha} and \emph{beta} at the same time,
#' and estimate with the formula \deqn{\hat{\mu}_{\alpha_{**},\beta} = \frac{1}{n_2}\sum\limits_{i = 1}^{n_2}\frac{f(x_{i})p(x_{i}) - {\beta}^{\mathsf{T}}\bigl(\boldsymbol{q}(x_{i}) - p(x_{i})\boldsymbol{1}\bigr)}{q_{\alpha_{**}}(x_{i})}.}
#' FALSE indicates optimizing \emph{alpha} only, and estimate with the formula \deqn{\hat{\mu}_{\alpha_{*}}= \frac{1}{n_2}\sum\limits_{i = 1}^{n_2}\frac{f(x_{i})p(x_{i})}{q_{\alpha_{*}}(x_{i})}.}
#' @param opt.info logical value indicating whether to save the returned value of the optimization procedure. See \code{\link{penoptpersp}} and \code{\link{penoptpersp.alpha.only}} for the returned value.
#' @param opt.param a list of paramters for the  damped Newton method with backtracking line search \describe{
#'   \item{reltol}{relative tolerence in dampednewton step, default to be \eqn{10^{-2}}}
#'   \item{relerr}{Newton step stop when within (1+\emph{relerr}) of minimum variance, default to be \eqn{10^{-3}}}
#'   \item{rho0}{initial value for \emph{rho}, default to be 1}
#' } Only need to supply part of the list to change the default value.
#' @inheritParams do.mixture.sample
#' @return  a list of \describe{
#'   \item{mu.hat}{the estimate for mu}
#'   \item{sd.hat}{estimated sd of mu.hat}
#'   \item{alpha.opt}{the estimated optimal alpha}
#'   \item{beta.opt}{if cv = TRUE, the estimated optimal beta}
#'   \item{opt.info}{if opt.info = TRUE, also return the list(x=x, y=y, z=z, alpha=alpha,beta=beta,rho=rho,f=f}
#'   \item{}{rhopen=rhopen, outer=log(rho0/rho,mu), relerr = relerr, alphasum = sum(alpha))}
#'   \item{}{from the optimization after batch 1.}
#' }
#' @details Estimate \eqn{E_p f} with a two step importance sample procedure. See He & Owen(2014) for details.
#' @references Boyd, S. and Vandeberghe, L. (2004). \emph{Convex optimization.} Cambridge University Press, Cambridge
#' @examples
#' library(optismixture)
#' seed <- 1
#' p <- 5
#' rho <- 1/2
#' gamma <- 2.4
#' sigma.dvar <- function(rho, p){
#'    sigma <- matrix(0, p, p)
#'    for(i in 1:(p-1)){
#'        for(j in (i+1):p){
#'            sigma[i,j] <- rho^(abs(i-j))
#'        }
#'    }
#'    sigma <- sigma + t(sigma)
#'    diag(sigma) <- 1
#'    return(sigma)
#'}
#' sigma <- sigma.dvar(rho, p)
#' batch.size <- c(10^4, 1002)
#' j.vec <- 2^-(seq(1,5,1))
#' eps.status <- 1
#' eps.safe <- 0.1
#' ## initialization and construct derived parameters
#' mu <- rep(0, p)
#' x0 <- matrix(1, 1, p)
#' x0.mat <- rbind(rep(1,p), rep(-1, p))
#' j.mat <- data.frame(centerid = rep(1:dim(x0.mat)[1], each = length(j.vec)),
#'                     variance = rep(j.vec, 2))
#' J <- dim(j.mat)[1] + 1
#' eps <- rep(0.1/J, J)
#' mixture.param <- list(x0 = x0, x0.mat = x0.mat, p = p,
#' sigma = sigma, gamma = gamma, j.mat = j.mat, J = J)
#' f <- function(x, j, mixture.param){
#'     f1 <- function(x, mixture.param){
#'         x0 <- mixture.param$x0
#'         gamma <- mixture.param$gamma
#'         return(sum((x - x0)^2)^(-gamma/2))
#'     }
#'     return(apply(x, 1, f1, mixture.param))
#' }
#' dq <- function(x, j, mixture.param){
#'     centerid <- mixture.param$j.mat[j, 1]
#'     j.param <- mixture.param$j.mat[j, 2]
#'     return(mvtnorm::dmvnorm(x, mixture.param$x0.mat[centerid,], j.param*diag(mixture.param$p)))
#' }
#' dp <- function(x, mixture.param){
#'     return(mvtnorm::dmvnorm(x, rep(0, mixture.param$p), mixture.param$sigma))
#' }
#' rq <- function(n, j, mixture.param){
#'     centerid <- mixture.param$j.mat[j, 1]
#'     j.param <- mixture.param$j.mat[j,2]
#'     return(mvtnorm::rmvnorm(n, mixture.param$x0.mat[centerid, ], j.param*diag(mixture.param$p)))
#' }
#' rp <- function(n, mixture.param){
#'     mu <- rep(0, mixture.param$p)
#'     sigma <- mixture.param$sigma
#'     return(mvtnorm::rmvnorm(n, mu, sigma))
#' }
#' a <- batch.estimation(seed, batch.size, mixture.param, eps, cv = FALSE,
#' fname = "f", rpname = "rp", rqname = "rq", dpname = "dp", dqname = "dq")
batch.estimation <- function(seed, batch.size, mixture.param, eps = rep(0.1/mixture.param$J, mixture.param$J), fname = "f", rpname = "rp", rqname = "rq", dpname = "dp", dqname = "dq", cv = TRUE, opt.info = FALSE, opt.param = list(reltol = 1e-3, relerr = 1e-3, rho0 = 1, maxin = 20, maxout = 30)){
    if(length(batch.size) != 2){
        stop("batch.size should be length 2")
    }
    mixture.param.names <- c("J", "p")
    for(name in mixture.param.names){
        if(length(mixture.param) < 2 || !exists(name, where = mixture.param))
            stop(do.call(paste, as.list(c("mixture.param should contain", mixture.param.names))))
    }

    J <- mixture.param$J
    K <- J - 1                          #number of control variates
    B <- length(batch.size)
    N <- sum(batch.size)
    sample.density.list <- list()
    sample.density.list$x.mat <- matrix(data = NA, N, mixture.param$p)
    sample.density.list$fx.mat <- matrix(data = NA, N, 1)
    sample.density.list$qx.mat <- matrix(data = NA, N, J)
    sample.density.list$px.mat <- matrix(data = NA, N, 1)

    alpha.mat <- matrix(data = NA, B, J)
    beta.mat <- matrix(data = NA, B, K)
    alpha.init <- get.initial.alpha(eps, J)
    alpha.opt <- alpha.init
    alpha.opt.norm <- alpha.opt/sum(alpha.opt)
    beta.init <- rep(0, K)
    beta.opt <- beta.init
    opt.list <- list()                #optimize alpha and beta simultaneously

    Y.mat <- matrix(NA, N, 1)
    X.mat <- matrix(NA, N, K)

    opt.input <- list()
    opt.input$Y <- matrix(data = NA, N, 1)
    opt.input$Z <- matrix(data = NA, N, J)
    opt.input$X <- matrix(data = NA, N, K)

    for(b in 1:B){
        index.b <- get.index.b(batch.size, b)
        mixture.sample <- do.mixture.sample(seed, b, batch.size[b], J, mixture.param, alpha.opt.norm, fname, rpname, rqname, dpname, dqname)
        alpha.mat[b,] <- mixture.sample$alpha
        weights <- c(rep(0, b - 1), 1)

        sample.density.list$x.mat[index.b,] <- mixture.sample$x
        sample.density.list$fx.mat[index.b,] <- mixture.sample$fx
        sample.density.list$qx.mat[index.b,] <- mixture.sample$qx
        sample.density.list$px.mat[index.b,] <- mixture.sample$px

        opt.input$Y[index.b,] <- sample.density.list$fx.mat[index.b,] * sample.density.list$px.mat[index.b,]
        opt.input$X[index.b,] <- sample.density.list$qx.mat[index.b,1:K] - sample.density.list$px.mat[index.b,]
        opt.input$Z[index.b,] <- drop(sample.density.list$qx.mat[index.b,]%*%alpha.mat[b,]) * sample.density.list$qx.mat[index.b,] * batch.size[b]
        if(b < B){
            if(cv){
                ## tmp <- list(x = as.matrix(opt.input$X[get.index.b(batch.size, b),]), y = opt.input$Y[get.index.b(batch.size, b),], z = as.matrix(opt.input$Z[get.index.b(batch.size, b),]), a0 = alpha.opt, b0 = beta.opt, eps = eps,  reltol = opt.param$reltol, relerr = opt.param$relerr, rho0 = opt.param$rho0)
                ## save(tmp, file = "../errordata/error.RData")
                ## stop("saved")
                opt.list[[b]] <- penoptpersp(x = as.matrix(opt.input$X[get.index.b(batch.size, b),]), y = opt.input$Y[get.index.b(batch.size, b),], z = as.matrix(opt.input$Z[get.index.b(batch.size, b),]), a0 = alpha.opt, b0 = beta.opt, eps = eps,  reltol = opt.param$reltol, relerr = opt.param$relerr, rho0 = opt.param$rho0, maxin = opt.param$maxin, maxout = opt.param$maxout)
                ## saveRDS(opt.list[[b]], file = "../errorworkspace/optlist.rds")
                ## stop("saved")

                                        #use alpha.opt beta.opt from previous iterations as starting points
                alpha.opt <- opt.list[[b]]$alpha
                alpha.opt.norm <- alpha.opt/sum(alpha.opt)
            }else if(!cv){
                ## tmp <- list(y = opt.input$Y[get.index.b(batch.size, b),], z = as.matrix(opt.input$Z[get.index.b(batch.size, b),]), a0 = alpha.opt, b0 = beta.opt, eps = eps,  reltol = opt.param$reltol, relerr = opt.param$relerr, rho0 = opt.param$rho0)
                ## save(tmp, file = "../errordata/error.RData")
                ## stop("saved")
                opt.list[[b]] <- penoptpersp.alpha.only(y = opt.input$Y[get.index.b(batch.size, b),], z = as.matrix(opt.input$Z[get.index.b(batch.size, b),]), a0 = alpha.opt, eps = eps, reltol = opt.param$reltol,relerr = opt.param$relerr, rho0 = opt.param$rho0, maxin = opt.param$maxin, maxout = opt.param$maxout)

                                        #use alpha.opt beta.opt from previous iterations as starting points
                alpha.opt <- opt.list[[b]]$alpha
                alpha.opt.norm <- alpha.opt/sum(alpha.opt)
                beta.opt <- opt.list[[b]]$beta
            }
            if(identical(alpha.opt.norm, alpha.init/sum(alpha.init))){
                warning("alpha may have not converged yet")
            }
        }
        qalpha <- drop(mixture.sample$qx %*% mixture.sample$alpha)
        Y.mat[index.b,] <- mixture.sample$fx * mixture.sample$px/ qalpha
        X.mat[index.b,] <-  (mixture.sample$qx[, -dim(mixture.sample$qx)[2]] - drop(mixture.sample$px))/qalpha
    }
    Y <- Y.mat[get.index.b(batch.size, B),]
    X <- X.mat[get.index.b(batch.size, B),]
    if(!cv){
        ## No control variates
        mu.hat <- mean(Y)
        nvec <- round(alpha.mat[b,]*batch.size[b])
        sd.hat <- sqrt(get.var(Y, nvec))
    }
    if(cv){
        ## control variates via normal regression
        MLR <- lm(Y ~ X)
        beta.lm.coef <- MLR$coef
        beta.lm.coef[is.na(beta.lm.coef)] <- 0
        mu.hat <- beta.lm.coef[[1]]
        sd.hat <- summary(MLR)$coef[1,2]
        beta.opt <- beta.lm.coef
    }
    result <- list(mu.hat = mu.hat, sd.hat = sd.hat, alpha.opt = alpha.opt.norm)
    if(cv){
        result[["beta.opt"]] <- beta.opt
    }
    if(opt.info){
        result[["opt.info"]] = opt.list[[B-1]]
    }
    return(result)
}

#' For a given mixture weight alpha, use importance sample with or withour control variates for estimation
#' @export
#' @param N sample size
#' @param alpha the mixture weight, sum up to 1
#' @inheritParams batch.estimation
mixture.is.estimation <- function(seed, N, mixture.param, alpha, fname = "f", rpname = "rp", rqname = "rq", dpname = "dp", dqname = "dq", cv = TRUE){
    mixture.param.names <- c("J", "p")
    for(name in mixture.param.names){
        if(length(mixture.param) < 2 || !exists(name, where = mixture.param))
            stop(do.call(paste, as.list(c("mixture.param should contain", mixture.param.names))))
    }
    N <- as.integer(N)
    if(length(N) != 1 | N <= 0){
        stop("N should be a positive integer")
    }
    if(any(alpha < 0)){stop("alpha need to be positive")}

    J <- mixture.param$J
    K <- J - 1                          #number of control variates
    sample.density.list <- list()
    sample.density.list$x.mat <- matrix(data = NA, N, mixture.param$p)
    sample.density.list$fx.mat <- matrix(data = NA, N, 1)
    sample.density.list$qx.mat <- matrix(data = NA, N, J)
    sample.density.list$px.mat <- matrix(data = NA, N, 1)
    alpha.opt.norm <- alpha/sum(alpha)
    Y.mat <- matrix(NA, N, 1)
    X.mat <- matrix(NA, N, K)
    mixture.sample <- do.mixture.sample(seed, 1, N, J, mixture.param, alpha.opt.norm, fname, rpname, rqname, dpname, dqname)
    alpha.mat <- mixture.sample$alpha
    weights <- 1
    sample.density.list$x.mat <- mixture.sample$x
    sample.density.list$fx.mat <- mixture.sample$fx
    sample.density.list$qx.mat <- mixture.sample$qx
    sample.density.list$px.mat <- mixture.sample$px
    qalpha <- drop(mixture.sample$qx %*% mixture.sample$alpha)
    Y <- mixture.sample$fx * mixture.sample$px/ qalpha
    X <-  (mixture.sample$qx[, -dim(mixture.sample$qx)[2]] - drop(mixture.sample$px))/qalpha
    if(!cv){
        ## No control variates
        mu.hat <- mean(Y)
        nvec <- round(alpha.mat*N)
        sd.hat <- sqrt(get.var(Y, nvec))
    }
    if(cv){
        ## control variates via normal regression
        MLR <- lm(Y ~ X)
        beta.lm.coef <- MLR$coef
        beta.lm.coef[is.na(beta.lm.coef)] <- 0
        mu.hat <- beta.lm.coef[[1]]
        sd.hat <- summary(MLR)$coef[1,2]
        beta.opt <- beta.lm.coef
    }

    result <- list(mu.hat = mu.hat, sd.hat = sd.hat, alpha = alpha.mat)
    if(cv){
        result[["beta"]] <- beta
    }
    return(result)
}

#' Do plain monte carlo with target density
#' @export
#' @param plainmc.N number of samples
#' @inheritParams do.mixture.sample
#' @return  a list of \describe{
#'   \item{plainmc.N}{number of samples for the plain monte carlo}
#'   \item{mu.hat}{estimated \eqn{E_p f} from plain monte carlos samples}
#'   \item{sd.hat}{estimated sd for \emph{mu.hat}}
#' }
do.plain.mc <- function(plainmc.N, mixture.param, fname = "f", rpname = "rp"){
    x <- eval(parse(text = rpname))(plainmc.N, mixture.param)
    fx <- eval(parse(text = fname))(x, mixture.param$J, mixture.param)
    mu.hat <- mean(fx)
    sd.hat <- sd(fx)
    return(list(plainmc.N = plainmc.N, mu.hat = mu.hat, sd.hat = sd.hat))
}

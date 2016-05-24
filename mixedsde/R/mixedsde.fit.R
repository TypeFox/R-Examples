#' Estimation Of The Random Effects In Mixed Stochastic Differential Equations
#' 
#' @description Estimation of the random effects \eqn{(\alpha_j, \beta_j)} and of their density, parametrically or nonparametrically in the mixed SDE
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma a(X_j(t)) dW_j(t)}.
#' @param times vector of observation times
#' @param X matrix of the M trajectories (each row is a trajectory with as much columns as observations)
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross)
#' @param random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects
#' @param fixed fixed effect in the drift: value of the fixed effect when there is only one random effect and it is not estimated, 0 otherwise
#' @param estim.method estimation method: 'paramML' for a Gaussian parametric estimation by maximum likelihood, 'paramBayes' for a Gaussian parametric Bayesian estimation or 'nonparam' for a non-parametric estimation
#' @param gridf if nonparametric estimation: grid of values on which the density is estimated, a matrix with two rows if two random effects; NULL by default and
#' then grid is chosen as a function of the estimated values of the random effects. For the plots this grid is used.
#' @param estim.fix default 0, 1 if random = 1 or 2, method = 'paramML' and an estimator of the fixed parameter is needed (to lead the nonparametric estimation after for example)
#' @param prior if method = 'paramBayes', list of prior parameters: mean and variance of the Gaussian prior on the mean mu, shape and scale of the inverse Gamma prior for the variances omega, shape and scale of the inverse Gamma prior for sigma  
#' @param nMCMC if method = 'paramBayes', number of iterations of the MCMC algorithm
#' @return
#'
#' \item{index}{is the vector of subscript in \eqn{1,...,M} where the estimation of \eqn{phi} has been done,  most of the time \eqn{index= 1:M}}
#' \item{estimphi}{matrix of estimators of \eqn{\phi=\alpha, or \beta, or (\alpha,\beta)} from the efficient statitics (see \code{\link{UV}}), matrix of two lines if random =c(1,2), numerical type otherwise}
#' \item{estim.fixed}{if estim.fix is TRUE and random = 1 or 2, estimator of \eqn{\phi=\alpha, or \beta} from the efficient statitics (see \code{\link{UV}}), 0 otherwise}
#' \item{gridf}{grid of values on which the estimated is done for the nonparametric method, otherwise, grid used for the plots, matrix form}
#' \item{estimf}{estimator of the density of \eqn{\phi} from a kernel estimator from package: stats, function: density, or package: MASS, function: kde2D. Matrix form: one line if one random effect or square matrix otherwise}
#' If there is a truncation threshold in the estimator
#' \item{cutoff}{the binary vector of cutoff, FALSE otherwise}
#' \item{estimphi.trunc}{troncated estimator of \eqn{\phi}, vector or matrix of 0 if we do not use truncation, matrix of two lines if random =c(1,2), numerical type otherwise}
#' \item{estimf.trunc}{troncated estimator of the density of \eqn{\phi}, vector or matrix of 0 if we do not use truncation, matrix if random =c(1,2), numerical type otherwise}
#' For the parametric maximum likelihood estimation 
#' \item{mu}{estimator of the mean of the random effects normal density, 0 if we do nonparametric estimation}
#' \item{omega}{estimator of the standard deviation of the random effects normal density, 0 if we do nonparametric estimation}
#' \item{bic}{BIC criterium, 0 if we do nonparametric estimation}
#' \item{aic}{AIC criterium, 0 if we do nonparametric estimation}
#' \item{model}{initial choice}
#' \item{random}{initial choice}
#' \item{fixed}{initial choice}
#' \item{times}{initial choice}
#' \item{X}{initial choice}
#' 
#' For the parametric Bayesian estimation 
#' \item{alpha}{posterior samples (Markov chain) of \eqn{\alpha}}
#' \item{beta}{posterior samples (Markov chain) of \eqn{\beta}}
#' \item{mu}{posterior samples (Markov chain) of \eqn{\mu}}
#' \item{omega}{posterior samples (Markov chain) of \eqn{\Omega}}
#' \item{sigma2}{posterior samples (Markov chain) of \eqn{\sigma^2}}
#' \item{model}{initial choice}
#' \item{random}{initial choice}
#' \item{burnIn}{proposal for burn-in period}
#' \item{thinning}{proposal for thinning rate}
#' \item{prior}{initial choice or calculated by the first 10\% of series}
#' \item{times}{initial choice}
#' \item{X}{initial choice}
#' \item{ind.4.prior}{in the case of calculation of prior parameters: the indices of used series}
#' @details
#' Estimation of the random effects density from M independent trajectories of the SDE (the Brownian motions \eqn{W_j} are independent), with linear drift. Two diffusions are implemented, with one or two random effects:
#' \subsection{Ornstein-Uhlenbeck model (OU)}{
#' 
#' If random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma dW_j(t)  } 
#' 
#' If random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma dW_j(t)  }
#' 
#' If random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma dW_j(t)  } 
#' }
#' \subsection{Cox-Ingersoll-Ross model (CIR)}{
#' 
#' If random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma \sqrt{X_(t)} dWj_(t)  } 
#' 
#' If random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma \sqrt{X_j(t)} dW_j(t)  } 
#' 
#' If random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma \sqrt{X_j(t)}  dW_j(t)  } 
#'}
#' The nonparametric method estimates the density of the random effects with a kernel estimator (one-dimensional or two-dimensional density).
#' The parametric method estimates the mean and standard deviation of the Gaussian distribution of the random effects. 

#' @examples
#'# Frequentist estimation
#' # Two random effects
#' model = 'CIR'; T <- 10
#' delta <- 0.1; M <- 100 # delta <- 0.001 and M <- 200 would yield good results
#' N <- floor(T/delta); sigma <- 0.01 ;
#' random <- c(1,2); density.phi <- 'gammainvgamma2'; param<- c(1.8, 0.8, 8, 0.05);  
#' simu <- mixedsde.sim(M=M, T=T, N=N, model=model,random=random, density.phi=density.phi, 
#'                param=param, sigma=sigma, invariant = 1)
#' X <- simu$X ; phi <- simu$phi; times <- simu$times
#' estim.method<- 'nonparam'
#' estim <- mixedsde.fit(times=times, X=X, model=model, random=random, estim.method= 'nonparam') 
#' #To stock the results of the function, use method \code{out}
#' #which put the outputs of the function on a list according to the frequentist or 
#' # Bayesian approach.
#' outputsNP <- out(estim)
#' 
#' \dontrun{
#' plot(estim)}
#' # It represents the bidimensional density, the histogram of the first estimated random 
#' # effect \eqn{\alpha} with the  marginal of \eqn{\hat{f}} from the first coordonate which 
#' # estimates  the density of \eqn{\alpha}. And the same for the second random effect 
#' # \eqn{\beta}. More, it plots a qq-plot for the sample of estimator of the random effects 
#' # compared with the quantiles of a Gaussian sample with the same mean and standard deviation.
#' 
#' summary(estim)
#' print(estim)
#' # Validation 
#' # If numj is fixed by the user: this function simulates Mrep =100 (by default) new 
#' # trajectories with the value of the estimated random effect. Then it plots on the 
#' # left graph the Mrep new trajectories \eqn{(Xnumj^{k}(t1), ... Xnumj^{k}(tN)),
#' # k= 1, ... Mrep} with in red the true trajectory \eqn{(Xnumj(t1), ... Xnumj(tN))}. 
#' #The right graph is a qq-plot of the quantiles of samples 
#' # \eqn{(Xnumj^{1}(ti), ... Xnumj^{Mrep}(ti))}
#' # for each time \eqn{ti} compared with the uniform quantiles. The outputs of the function  
#' # are: a matrix \code{Xnew} dimension Mrepx N+1, vector of quantiles \code{quantiles} length 
#' # N and the number of the trajectory for the plot \code{plotnumj= numj} 
#' # If numj is not precised by the user, then, this function simulates Mrep =100 (by default) 
#' # new trajectories for each estimated random effect. Then left graph is a plot of the Mrep 
#' # new trajectories \eqn{(Xj^{k}(t1), ... Xj^{k}(tN)), k= 1, ... Mrep}
#' #for a randomly chosen number j with in red the true trajectory \eqn{(Xj(t1), ... Xj(tN))}. 
#' #The right graph is a qq-plot of the quantiles of samples \eqn{(Xj^{1}(ti), ... Xj^{Mrep}(ti))},
#' # for the same j and for each time \eqn{ti}. The outputs of the function are: a list of 
#' # matrices \code{Xnew} length M, matrix of quantiles \code{quantiles} dimension MxN
#' # and the number of the trajectory for the plot \code{plotnumj} 
#' 
#' validation <- valid(estim,  numj=floor(runif(1,1,M)))
#' 
#' # Parametric estimation
#' estim.method<-'paramML'
#' estim_param <- mixedsde.fit(times= times, X= X, model= model, random= random, 
#'            estim.method = 'paramML') 
#' outputsP <- out(estim_param)
#' 
#' #plot(estim_param)
#' summary(estim_param)
#' 
#' # Prediction for the frequentist approach
#' # This function uses the estimation of the density function to simulate a 
#' # new sample of random effects according to this density. If \code{plot.pred =1} (default)
#' # is plots on the top the predictive random effects versus the estimated random effects
#' # from the data. On the bottom, the left graph is the true trajectories, on the right
#' #the predictive trajectories and the empiric prediciton intervals at level 
#' # \code{level=0.05} (defaut). The function return on a list the prediction of phi 
#' # \code{phipred}, the prediction of X \code{Xpred}, and the indexes of the 
#' # corresponding true trajectories \code{indexpred} 
#' 
#' # Not run
#' \dontrun{
#' test1 <- pred(estim,  invariant  = 1)
#' test2 <- pred(estim_param, invariant  = 1)
#' }
#' # More graph
#' fhat <- outputsNP$estimf  
#' fhat_trunc <- outputsNP$estimf.trunc 
#' fhat_param <- outputsP$estimf
#' 
#' gridf <- outputsNP$gridf; gridf1 <- gridf[1,]; gridf2 <- gridf[2,]
#' 
#' marg1 <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat,1,sum) 
#' marg1_trunc <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_trunc,1,sum) 
#' marg2 <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat,2,sum)
#' marg2_trunc <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_trunc,2,sum)
#' 
#' marg1_param <- ((max(gridf2)-min(gridf2))/length(gridf2))*apply(fhat_param,1,sum) 
#' marg2_param <- ((max(gridf1)-min(gridf1))/length(gridf1))*apply(fhat_param,2,sum)
#' f1 <-  (gridf1^(param[1]-1))*exp(-gridf1/param[2])/((param[2])^param[1]*gamma(param[1]))
#' f2 <-  (gridf2^(-param[3]-1)) * exp(-(1/param[4])*(1/gridf2)) * 
#'  ((1/param[4])^param[3])*(1/gamma(param[3]))
#' par(mfrow=c(1,2))
#' plot(gridf1,f1,type='l', lwd=1,  xlab='', ylab='')
#' lines(gridf1,marg1_trunc,col='blue', lwd=2)
#' lines(gridf1,marg1,col='blue', lwd=2, lty=2)
#' lines(gridf1,marg1_param,col='red', lwd=2, lty=2)
#' plot(gridf2,f2,type='l', lwd=1, xlab='', ylab='')
#' lines(gridf2,marg2_trunc,col='blue', lwd=2)
#' lines(gridf2,marg2,col='blue', lwd=2, lty=2)
#' lines(gridf2,marg2_param,col='red', lwd=2, lty=2)
#' 
#' cutoff <- outputsNP$cutoff
#' phihat <- outputsNP$estimphi 
#' phihat_trunc <- outputsNP$estimphi.trunc
#' par(mfrow=c(1,2))
#' plot.ts(phi[1,], phihat[1,], xlim=c(0, 15), ylim=c(0,15), pch=18); abline(0,1)
#' points(phi[1,]*(1-cutoff), phihat[1,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20),pch=18, col='red'); 
#' abline(0,1)
#' plot.ts(phi[2,], phihat[2,], xlim=c(0, 15), ylim=c(0,15),pch=18); abline(0,1)
#' points(phi[2,]*(1-cutoff), phihat[2,]*(1-cutoff), xlim=c(0, 20), ylim=c(0,20),pch=18, col='red'); 
#' abline(0,1)
#' 
#' # one random effect: 
#' \dontrun{
#' model <-'OU'
#' random <- 1
#' M <- 80; T <- 100  ; N <- 2000 
#' sigma <- 0.1 ; X0 <- 0
#' density.phi <- 'normal'
#' fixed <- 2 ; param <- c(1, 0.2)
#'#-------------------
#'#- simulation
#' simu <- mixedsde.sim(M,T=T,N=N,model=model,random=random, fixed=fixed,density.phi=density.phi, 
#'                param=param, sigma=sigma, X0=X0)
#' X <- simu$X
#' phi <- simu$phi
#' times <-simu$times
#' plot(times, X[10,], type='l')
#' 
#' #- parametric estimation
#' estim.method<-'paramML'
#' estim_param <- mixedsde.fit(times, X=X, model=model, random=random, estim.fix= 1, 
#'                estim.method=estim.method)
#' outputsP <- out(estim_param)
#' estim.fixed <- outputsP$estim.fixed #to compare with fixed
#' #or 
#' estim_param2 <- mixedsde.fit(times, X=X, model=model, random=random, fixed = fixed, 
#'              estim.method=estim.method)
#' outputsP2 <- out(estim_param2)
#' #- nonparametric estimation
#' estim.method <- 'nonparam'
#' estim <- mixedsde.fit(times, X, model=model, random=random, fixed = fixed, 
#'            estim.method=estim.method) 
#' outputsNP <- out(estim)
#'
#' plot(estim)
#' print(estim)
#' summary(estim)
#'
#' plot(estim_param)
#' print(estim_param)
#' summary(estim_param)
#'
#' valid1 <- valid(estim,  numj=floor(runif(1,1,M)))
#' test1 <- pred(estim )
#' test2 <- pred(estim_param)
#'}
#'
#' # Parametric Bayesian estimation 
#' # one random effect
#' random <- 1; sigma <- 0.1; fixed <- 5; param <- c(3, 0.5)
#' sim <- mixedsde.sim(M = 20, T = 1, N = 50, model = 'OU', random = random, fixed = fixed,
#'        density.phi = 'normal',param= param, sigma= sigma, X0 = 0, op.plot = 1)
#' 
#' # here: only 100 iterations for example - should be much more!
#' prior <- list( m = c(param[1], fixed), v = c(param[1], fixed), alpha.omega = 11,
#'             beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
#' estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random, 
#'            estim.method = 'paramBayes', prior = prior, nMCMC = 100) 
#' 
#' validation <- valid(estim_Bayes, numj = 10)
#' plot(estim_Bayes)
#' outputBayes <- out(estim_Bayes)
#' summary(outputBayes)
#' (results_Bayes <- summary(estim_Bayes))
#' plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi)
#' print(estim_Bayes)
#' \dontrun{
#' pred.result <- pred(estim_Bayes)
#' summary(out(pred.result))
#' plot(pred.result)
#' 
#' pred.result.trajectories <- pred(estim_Bayes, trajectories = TRUE)
#' }
#' # second example
#' \dontrun{
#' random <- 2; sigma <- 0.2; fixed <- 5; param <- c(3, 0.5)
#' sim <- mixedsde.sim(M = 20, T = 1, N = 100, model = 'CIR', random = random, 
#'         fixed = fixed, density.phi = 'normal',param = param, sigma = sigma, X0 = 0.1, op.plot = 1)
#' 
#' prior <- list(m = c(fixed, param[1]), v = c(fixed, param[1]), alpha.omega = 11, 
#'          beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
#'
#' estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'CIR', random = random, 
#'                  estim.method = 'paramBayes', prior = prior, nMCMC = 1000) 
#' 
#' pred.result <- pred(estim_Bayes)
#' }
#' 
#' 
#' 
#' # invariant case
#' \dontrun{
#' random <- 1; sigma <- 0.1; fixed <- 5; param <- c(3, 0.5)
#' sim <- mixedsde.sim(M = 50, T = 5, N = 100, model = 'OU', random = random, fixed = fixed, 
#'            density.phi = 'normal',param = param, sigma = sigma, invariant = 1, op.plot = 1)
#' 
#' prior <- list(m = c(param[1], fixed), v = c(param[1], 1e-05), alpha.omega = 11, 
#'        beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
#' estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random, 
#'        estim.method = 'paramBayes', prior = prior, nMCMC = 100) 
#' plot(estim_Bayes)
#'
#' pred.result <- pred(estim_Bayes, invariant = 1)
#' pred.result.traj <- pred(estim_Bayes, invariant = 1, trajectories = TRUE)
#'}

#' 
#' @keywords estimation
#' @references For the parametric estimation see:
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}
#' 
#' Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model.  S. Hermann, K. Ickstadt and C. Mueller, \emph{appearing in: Applied Stochastic Models in Business and Industry 2016}.
#'
#' For the nonparametric estimation see:
#' 
#' Nonparametric estimation for stochastic differential equations with random effects, F. Comte, V. Genon-Catalot and A. Samson, \emph{Stochastic Processes and Their Applications 2013}, Vol 7, \bold{2522--2551}
#' 
#' Estimation for stochastic differential equations with mixed effects, V. Genon-Catalot and C. Laredo 2014 \emph{e-print: hal-00807258 }
#' 
#' Bidimensional random effect estimation in mixed stochastic differential model, C. Dion and V. Genon-Catalot,  \emph{Stochastic Inference for Stochastic Processes 2015, Springer Netherlands}, \bold{1--28}

mixedsde.fit <- function(times, X, model = c("OU", "CIR"), random, fixed = 0, estim.fix = 0, estim.method = c("nonparam", "paramML", "paramBayes"), 
    gridf = NULL, prior, nMCMC = NULL) {
    model <- match.arg(model)
    estim.method <- match.arg(estim.method)
    
    if (is.matrix(X)) {
        if (nrow(X) == length(times)) {
            X <- t(X)
        } else {
            if (ncol(X) != length(times)) {
                stop("length of times has to be equal to the columns of X")
            }
        }
    }
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    delta <- round(diff(times), 10)  #diff(times)[1]
    Tend <- times[length(times)]
    
    if (estim.method == "paramBayes") {
        if (model == "CIR") {
            if (any(X < 0)) {
                Xold <- X
                indices <- sapply(1:M, function(i) any(X[i, ] < 0))
                X <- X[!indices, ]
                message("attention: series ", indices, " are skipped for estimation because of negative values")
            }
        }
        
        if (missing(prior)) {
            ind.4.prior <- 1:max(3, ceiling(M/10))
            X.4.prior <- X[ind.4.prior, ]
            estimUV <- UV(X.4.prior, model, random = c(1, 2), fixed = 0, times)  # fixed is only used for random == 1
            U <- estimUV$U
            V <- estimUV$V
            deter <- lapply(V, det)
            index <- which((deter != Inf) & (deter != 0))
            Mindex <- length(index)
            V <- V[index]
            U <- U[, index]
            
            
            if (length(V) == 0) {
                l.prior <- length(ind.4.prior)
                while (length(V) == 0 && max(ind.4.prior) + l.prior - 1 < M) {
                  ind.4.prior <- ind.4.prior + l.prior
                  X.4.prior <- X[ind.4.prior, ]
                  estimUV <- UV(X.4.prior, model, random = c(1, 2), fixed = 0, times)  # fixed is only used for random == 1
                  U <- estimUV$U
                  V <- estimUV$V
                  deter <- lapply(V, det)
                  index <- which((deter != Inf) & (deter != 0))
                  Mindex <- length(index)
                  V <- V[index]
                  U <- U[, index]
                }
            }
            if (length(V) == 0) {
                message("please specify prior parameters")
                prior <- list(m = c(1, 1), v = c(10, 10), alpha.omega = rep(3, length(random)), beta.omega = rep(10, length(random)) * 
                  2, alpha.sigma = 3, beta.sigma = 1 * 2)
                message("parameters are set to: ", unlist(prior))
            } else {
                
                A <- matrix(0, 2, Mindex)
                for (j in 1:Mindex) {
                  A[, j] <- (1/det(V[[j]])) * matrix(c(V[[j]][2, 2], -V[[j]][1, 2], -V[[j]][1, 2], V[[j]][1, 1]), 2, 2) %*% U[, j]
                }
                estimphi <- A
                
                mu <- apply(estimphi, 1, mean)
                Omega <- apply(estimphi, 1, var)
                if (model == "OU") 
                  var.fun <- function(x) 1
                if (model == "CIR") 
                  var.fun <- function(x) x
                
                sigma2 <- mean(sapply(1:Mindex, function(i) mean(diff(X.4.prior[index[i], 2:K])^2 * (1/delta[1:(K - 2)]) * var.fun(1/X.4.prior[index[i], 
                  3:K]))))
                prior <- list(m = mu, v = abs(mu), alpha.omega = rep(3, length(random)), beta.omega = Omega[random] * 2, alpha.sigma = 3, 
                  beta.sigma = sigma2 * 2)
                message("attention: series ", ind.4.prior, " are used for prior parameter calculation")
            }
        } else {
            ind.4.prior <- M + 1
        }
        
        res <- BayesianNormal(times, X[-ind.4.prior, ], model, prior, start = list(mu = prior$m, sigma = prior$beta.sigma/(prior$alpha.sigma - 
            1)), random, nMCMC)
        he <- diagnostic(res, random)
        return(new(Class = "Bayes.fit", prior = prior, alpha = as.matrix(res$alpha), beta = as.matrix(res$beta), random = random, mu = as.matrix(res$mu), 
            omega = as.matrix(res$omega), sigma2 = res$sigma2, burnIn = he$burnIn, thinning = he$thinning, model = model, times = times, 
            X = X, ind.4.prior = ind.4.prior))
        
    } else {
        
        if (sum(random) > 2) {
          
          if (missing(fixed)==0){
            message('the parameter fixed is not used because random = c(1,2)')
                  } 
          
          # estimation of sigma^2
          if (model == "OU") {
            
            meanU <- rep(0, M)
            for (i in 1:M) {
              meanU[i] <- mean((diff(X[i, 2:K])^2) * (1/delta[1:(K - 2)]))
            }
            sigma2 <- mean(meanU)
            Mindex <- M
            index <- 1 : Mindex
          }
          if (model == "CIR") {
            
            index <- which(rowSums(X <= 0) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
              #message("All the trajectories have non positive values the model CIR cannot be used")
              warning("All the trajectories have non positive values the model CIR cannot be used", 
                      call. = FALSE)
              estimf <- 0
              estimphi <- 0
              bic <- 0 
              aic <- 0
              gridf <- 0
              mu <- 0 
              omega <- 0
              cutoff <- 0
              sigma2 <- 0 
              estimf.trunc <- 0
              estimphi.trunc <- 0
              estim.fixed <- 0                   
            }
            if (Mindex > 0) {
              meanU <- rep(0, Mindex)
              for (i in 1:Mindex) {
                meanU[i] <- mean(diff(X[index[i], 2:K])^2 * (1/delta[1:(K - 2)]) * (1/X[index[i], 3:K]))
              }
              sigma2 <- mean(meanU)
            }
          }
            # -- computation of the sufficient statistics
            U <- matrix(0, 2, Mindex)
            V <- as.list(1:Mindex)
            b <- as.list(1:Mindex)
            
            estimUV <- UV(X[index,], model, random, fixed, times)
            U <- estimUV$U
            V <- estimUV$V
            
            deter <- lapply(V, det)
            index2 <-  which((deter != Inf) & (deter != 0)) # indexes in 1:Mindex 
            Mindex2 <- length(index2)
            V <- V[index2]
            U <- U[, index2]
            
            # estimation of the random effects phi
            A <- matrix(0, Mindex2, 2)
            for (j in 1:Mindex2) {
                A[j, ] <- (1/det(V[[j]])) * matrix(c(V[[j]][2, 2], -V[[j]][1, 2], -V[[j]][1, 2], V[[j]][1, 1]), 2, 2) %*% U[, j]
            }
            estimphi <- t(A)
            eigenvalues <- eigenvaluesV(V)
            
            if (is.null(gridf) == 1) {
                gridf <- matrix(0, 2, 500)
                gridf[1, ] <- seq(min(estimphi[1, ]) * 0.8, max(estimphi[1, ]) * 1.2, length = 500)
                gridf[2, ] <- seq(min(estimphi[2, ]) * 0.8, max(estimphi[2, ]) * 1.2, length = 500)
            }
            
            if (is.null(gridf) == 0) {
                gridf <- gridf
            }
            
            estim.fixed <- 0
            
            if (estim.method == "nonparam") {
                # troncation of the phi estimators
                
                kap <- 0.125
                cutoff <- apply(eigenvalues, 1, min) * (1/sigma2) > kap * sqrt(Tend)
                estimphi.trunc <- estimphi * matrix(c(cutoff, cutoff), 2, dim(estimphi)[2], byrow = TRUE)
                
                # estimation of the density
                
                estimf <- kde2d(estimphi[1, ], estimphi[2, ], n = length(gridf[1, ]), lims = c(min(gridf[1, ]), max(gridf[1, ]), min(gridf[2, 
                  ]), max(gridf[2, ])))$z
                
                if (sum(cutoff) >= 0.25 * Mindex2) {
                  estimf.trunc <- kde2d(estimphi.trunc[1, ], estimphi.trunc[2, ], n = length(gridf[1, ]), lims = c(min(gridf[1, ]), max(gridf[1, 
                    ]), min(gridf[2, ]), max(gridf[2, ])))$z
                }
                if (sum(cutoff) < 0.25 * Mindex2) {
                  message("More than 75 percents of the estimated values of the random effect have been put to zero")

                  estimf.trunc <- matrix(0, 2, length(gridf[1, ]))
                }
                
                bic <- 0
                aic <- 0
                mu <- 0
                omega <- 0
            }
            
            if (estim.method == "paramML") {
                # estimphi has Mindex colomns maximization of the likelihood
                
                Vsigma2 <- as.list(1:length(V))
                for (i in 1:length(V)) {
                  Vsigma2[[i]] <- V[[i]] * (1/sigma2)
                }
                Usigma2 <- U/sigma2
                
                res <- EstParamNormal(Usigma2, Vsigma2, K = K, random = random, estim.fix = 0)
                
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                omega <- res$omega
                # computation of the densities
                estimf1 <- dnorm(gridf[1, ], mean = mu[1], sd = abs(omega[1]))
                estimf2 <- dnorm(gridf[2, ], mean = mu[2], sd = abs(omega[2]))
                
                estimf <- estimf1 %*% t(estimf2)
                
                estimf.trunc <- estimf
                estimphi.trunc <- estimphi
                
                cutoff <- FALSE
            }
           
            index <- index[index2]

        }

#######################################################################################################################
        
        if (length(random) == 1) {
            
          if (model == "OU") {
            
            meanU <- rep(0, M)
            for (i in 1:M) {
              
              meanU[i] <- mean((diff(X[i, 2:K])^2) * (1/delta[1:(K - 2)]))
            }
            sigma2 <- mean(meanU)
            Mindex <- M
            index <- 1 : Mindex
          }
          
          if (model == "CIR") {
            index <- which(apply(X <= 0, 1, sum) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
              #message("All the trajectories have non positive values the model CIR cannot be used ")
              warning("All the trajectories have non positive values the model CIR cannot be used", 
                      call. = FALSE)
              
              estimf <- 0
              estimphi <- 0
              bic <- 0 
              aic <- 0
              gridf <- 0
              mu <- 0 
              omega <- 0
              cutoff <- 0
              sigma2<- 0 
              estimf.trunc <- 0
              estimphi.trunc <- 0
              estim.fixed <- 0                   
            }
            if (Mindex != 0){
              meanU <- rep(0, Mindex)
              for (i in 1:Mindex) {
                meanU[i] <- mean(diff(X[index[i], 2:K])^2 * (1/delta[1:(K - 2)]) * (1/X[index[i], 3:K]))
             }
              sigma2 <- mean(meanU)
            }
          }
          
          if (estim.method == "nonparam") {     
              
                if (estim.fix == 1) {
                  #message("wrong argument estim.fix with method nonparam, fixed as to be specify and estim.fix = 0")
                  warning("wrong argument estim.fix with method nonparam, fixed as to be specify and estim.fix = 0", 
                          call. = FALSE)
                  
                }
                
                if (estim.fix == 0) {
                  
                  U <- rep(0, Mindex)
                  V <- rep(0, Mindex)
                  estimUV <- UV(X[index,], model, random, fixed = fixed, times)
                  U <- estimUV$U
                  V <- estimUV$V
                  
                  # estimation of the random effect phi
                  index2 <- which((V != Inf) & (V != 0))
                  Mindex2 <- length(index2)
                  V <- V[index2]
                  U <- U[index2]
                  A <- U/V
                  estimphi <- A
                
                  
                  if (is.null(gridf) == 1) {
                    gridf <- seq(min(estimphi) * 0.8, max(estimphi) * 1.2, length = 500)
                  }
                  
                  if (is.null(gridf) == 0) {
                    gridf <- gridf
                  }
                  # estimation of the density
                  test <- density(estimphi, from = min(gridf), to = max(gridf), bw = "ucv", n = length(gridf))
                  
                  if (test$bw < 0.1) {
                    estimf <- density(estimphi, from = min(gridf), to = max(gridf), n = length(gridf))$y
                  }
                  if (test$bw >= 0.1) {
                    estimf <- test$y
                  }
                  
                  if (random == 2 & fixed == 0 & model == "OU") {
                    # troncation of the phi estimators
                    
                    kap <- 0.2
                    cutoff <- V * (1/sigma2) > kap * sqrt(Tend)
                    estimphi.trunc <- A * cutoff
                    
                    if (sum(cutoff) < 0.25 * Mindex2) {
                      message("More than 75 percents of the estimated values of the random effect have been put to zero")
                    }
                    
                    test2 <- density(estimphi.trunc, from = min(gridf), bw = "ucv", to = max(gridf), n = length(gridf))
                    
                    if (test2$bw < 0.1) {
                      estimf.trunc <- density(estimphi.trunc, from = min(gridf), to = max(gridf), n = length(gridf))$y
                      estimf.trunc <- matrix(estimf.trunc, 1, length(estimf.trunc), byrow = TRUE)
                      estimphi.trunc <- as.matrix(estimphi.trunc, 1, length(estimphi.trunc), byrow = TRUE)
                    }
                    if (test2$bw >= 0.1) {
                      estimf.trunc <- test2$y
                      estimf.trunc <- matrix(estimf.trunc, 1, length(estimf.trunc), byrow = TRUE)
                      estimphi.trunc <- matrix(estimphi.trunc, 1, length(estimphi.trunc), byrow = TRUE)
                    }
                    
                  } else {
                    
                    estimf.trunc <- matrix(estimf, 1, length(estimf), byrow = TRUE)
                    estimphi.trunc <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
                    cutoff <- FALSE
                  }
                  gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                  estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
                  estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
                  
                  bic <- 0
                  aic <- 0
                  mu <- 0
                  omega <- 0
                  estim.fixed <- 0
                }
              index <- index[index2]
            }
          if (estim.method == "paramML") {
            
                U12 <- matrix(0, 2, Mindex)
                V12 <- as.list(1:Mindex)
                b12 <- as.list(1:Mindex)
                
                estimUV12 <- UV(X[index, ], model, random = c(1, 2), fixed = 0, times)
                U12 <- estimUV12$U
                V12 <- estimUV12$V
                
                deter12 <- lapply(V12, det)
                index12 <- which((deter12 != Inf) & (deter12 != 0))
                Mindex12 <- length(index12)
                V12 <- V12[index12]
                U12 <- U12[, index12]
                # estimation of sigma2
#                 if (model == "OU") {
#                   
#                   meanU <- rep(0, Mindex12)
#                   for (i in 1:Mindex12) {
#                     
#                     meanU[i] <- mean((diff(X[index12[i], 2:K])^2) * (1/delta[1:(K - 2)]))
#                   }
#                   sigma2 <- mean(meanU)
#                 }
#                 if (model == "CIR") {
#                   index12 <- intersect(which(apply(X <= 0, 1, sum) == 0), index12)
#                   Mindex12 <- length(index12)
#                   meanU <- rep(0, Mindex12)
#                   for (i in 1:Mindex12) {
#                     meanU[i] <- mean(diff(X[index12[i], 2:K])^2 * (1/delta[1:(K - 2)]) * (1/X[index12[i], 3:K]))
#                   }
#                   sigma2 <- mean(meanU)
#                 }
                
                Vsigma212 <- as.list(1:length(V12))
                Usigma212 <- U12 * (1/sigma2)
                for (i in 1:length(V12)) {
                  Vsigma212[[i]] <- V12[[i]] * (1/sigma2)
                }
                
                
                if (estim.fix == 1) {
                  
                  res <- EstParamNormal(U = Usigma212, V = Vsigma212, K = K, random = random, estim.fix = 1)
                  bic <- res$BIChere
                  aic <- res$AIChere
                  mu <- res$mu[random]
                  omega <- res$omega
                  estim.fixed <- res$mu[(random == 1) + 1]
                  
                  U <- rep(0, Mindex)
                  V <- rep(0, Mindex)
                  estimUV <- UV(X[index,], model, random, fixed = estim.fixed, times)
                  V <- estimUV$V
                  U <- estimUV$U
                  # estimation of the random effect phi
                  index2 <- which((V != Inf) & (V != 0))
                  
                  V <- V[index2]
                  U <- U[index2]
                  A <- U/V
                  estimphi <- A
                  index <- index[index2]
                }
                
                if (estim.fix == 0) {
                  if (fixed == 0) {
                    print("be careful: estim.fix=0 and fixed=0 thus your fixed effect is 0 and it is not estimated")
                  }
                  
                  res <- EstParamNormal(U = Usigma212, V = Vsigma212, K = K, random, estim.fix = 0, fixed = fixed)
                  bic <- res$BIChere
                  aic <- res$AIChere
                  mu <- res$mu
                  omega <- res$omega
                  estim.fixed <- 0
                  
                  U <- rep(0, Mindex)
                  V <- rep(0, Mindex)
                  estimUV <- UV(X[index,], model, random, fixed = fixed, times)
                  V <- estimUV$V
                  U <- estimUV$U
                  # estimation of the random effect phi
                  index2 <- which((V != Inf) & (V != 0))
        
                  V <- V[index2]
                  U <- U[index2]
                  A <- U/V
                  estimphi <- A
                  index <- index[index2]
                }
                
                if (is.null(gridf) == 1) {
                  gridf <- seq(min(estimphi) * 0.8, max(estimphi) * 1.2, length = 500)
                }
                
                if (is.null(gridf) == 0) {
                  gridf <- gridf
                }
                
                # estimation of the density
                estimf <- matrix(dnorm(gridf, mean = mu, sd = omega), 1, length(gridf), byrow = TRUE)
                gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
                estimphi.trunc <- estimphi
                estimf.trunc <- estimf
                cutoff <- FALSE
                
            }
        }
    }
    return(new(Class = "Freq.fit", model = model, random = random, fixed = fixed, gridf = gridf, mu = mu, omega = omega, cutoff = cutoff, 
        sigma2 = sigma2, estimf.trunc = estimf.trunc, estimphi.trunc = estimphi.trunc, estimf = estimf, estimphi = estimphi, estim.fixed = estim.fixed, 
        estim.fix = estim.fix, index = index, bic = bic, aic = aic, times = times, X = X))
    
}



#' S4 class for the frequentist estimation results  
#'  
#' @slot model character 'OU' or 'CIR'
#' @slot random numeric 1, 2, or c(1,2)
#' @slot fixed numeric value of the fixed effect if there is one 
#' @slot gridf matrix of values on which the estimated is done
#' @slot mu numeric MLE estimator for parametric approach
#' @slot omega numeric  MLE estimator for parametric approach
#' @slot cutoff value of the cutoff if there is one
#' @slot sigma2 numeric estimated value of \eqn{\sigma^2}
#' @slot estimf.trunc matrix estimator of the density of \eqn{\phi} for the truncated estimateur of the random effects
#' @slot estimphi.trunc matrix truncated estimator of the random effects
#' @slot index index of the used trajectories
#' @slot estimphi matrix of the estimator of the random effects
#' @slot estimf estimator of the density of \eqn{\phi}
#' @slot estim.fixed estimator of the fixed parameter if option estim.fix = 1
#' @slot estim.fix 1 if the user asked for the estimation of fixed parameter 
#' @slot bic numeric bic 
#' @slot aic numeric aic
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable

setClass(Class = "Freq.fit", representation = representation(model = "character", random = "numeric", fixed = "numeric", gridf = "matrix", 
    mu = "numeric", omega = "numeric", cutoff = "logical", sigma2 = "numeric", estimf.trunc = "matrix", estimphi.trunc = "matrix", index = "numeric", 
    estimphi = "matrix", estimf = "matrix", estim.fixed = "numeric", estim.fix = "numeric", bic = "numeric", aic = "numeric", times = "numeric", X = "matrix"))



#' S4 class for the Bayesian estimation results
#' @slot sigma2 vector of posterior samples for \eqn{\sigma^2}
#' @slot mu matrix of posterior samples for \eqn{\mu}
#' @slot omega matrix of posterior samples for \eqn{\omega}
#' @slot alpha matrix of posterior samples for \eqn{\alpha}
#' @slot beta matrix of posterior samples for \eqn{\beta}
#' @slot random 1, 2 or c(1,2)
#' @slot burnIn proposal for the burn-in phase
#' @slot thinning proposal for the thinning rate
#' @slot model 'OU' or 'CIR'
#' @slot prior list of prior values, input variable or calculated by the first 10\% of series
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable
#' @slot ind.4.prior indices of series used for the prior parameter calculation, if prior knowledge is availabe it is set to M+1
#' 
setClass(Class = "Bayes.fit", representation = representation(sigma2 = "numeric", mu = "matrix", omega = "matrix", alpha = "matrix", beta = "matrix", 
    random = "numeric", burnIn = "numeric", thinning = "numeric", model = "character", prior = "list", times = "numeric", X = "matrix", 
    ind.4.prior = "numeric"))

#' S4 class for the Bayesian prediction results
#' @slot phi.pred matrix of predictive samples for the random effect
#' @slot Xpred matrix of predictive samples for observations
#' @slot coverage.rate amount of covering prediction intervals
#' @slot qu.u upper prediction interval bound
#' @slot qu.l lower prediction interval bound
#' @slot estim list of Bayes.fit object entries, storage of input variable
setClass(Class = "Bayes.pred", representation = representation(phi.pred = "matrix", Xpred = "matrix", coverage.rate = "numeric", qu.u = "numeric", 
    qu.l = "numeric", estim = "list"))


########################################################### OUTPUTS
#' Transfers the class object to a list
#' 
#' @description Method for the S4 classes
#' @param x Freq.fit, Bayes.fit or Bayes.pred class
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 

out <- function(x){
    sN <- slotNames(x)
    res <- lapply(sN, function(name) slot(x, name))
    names(res) <- sN
    res
}




############################################################### SUMMARY

#' Short summary of the results of class object Freq.fit
#' @description Method for the S4 class Freq.fit
#' @param object Freq.fit class
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#'
setMethod("summary", "Freq.fit", function(object) {
    if (object@bic != 0) {
        
        if (dim(object@gridf)[1] == 1) {
            
            if (object@estim.fix == 1) {
                print(matrix(c("BIC", object@bic, "AIC", object@aic), 2, 2, byrow = TRUE))
                print(matrix(c("kurtosis", kurtosis(object@estimphi[1, ]), "skewness", skewness(object@estimphi[1, ])), 2, 2, byrow = TRUE))
                print(matrix(c("sigma", sqrt(object@sigma2), "estim.fixed", object@estim.fixed, "empiric mean", mean(object@estimphi[1, 
                  ]), "MLE mean", object@mu, "empiric sd", sd(object@estimphi[1, ]), "MLE sd", object@omega), 6, 2, byrow = TRUE))
            }
            if (object@estim.fix == 0) {
                print(matrix(c("BIC", object@bic, "AIC", object@aic), 2, 2, byrow = TRUE))
                print(matrix(c("kurtosis", kurtosis(object@estimphi[1, ]), "skewness", skewness(object@estimphi[1, ])), 2, 2, byrow = TRUE))
                print(matrix(c("sigma", sqrt(object@sigma2), "empiric mean", mean(object@estimphi[1, ]), "MLE mean", object@mu, "empiric sd", 
                  sd(object@estimphi[1, ]), "MLE sd", object@omega), 5, 2, byrow = TRUE))
            }
            
            
        }
        
        if (dim(object@gridf)[1] == 2) {
            print(matrix(c("BIC", object@bic, "AIC", object@aic), 2, 2, byrow = TRUE))
            print(matrix(c("sigma", sqrt(object@sigma2)), 1, 2, byrow = TRUE))
            print(matrix(c("empiric mean 1", mean(object@estimphi[1, ]), "MLE mean 1", object@mu[1], "empiric sd 1", sd(object@estimphi[1, 
                ]), "MLE sd 1", object@omega[1], "kurtosis 1", kurtosis(object@estimphi[1, ]), "skewness 1", skewness(object@estimphi[1, 
                ])), 6, 2, byrow = TRUE))
            print(matrix(c("empiric mean 2", mean(object@estimphi[2, ]), "MLE mean 2", object@mu[2], "empiric sd 2", sd(object@estimphi[2, 
                ]), "MLE sd 2", object@omega[2], "kurtosis 2", kurtosis(object@estimphi[2, ]), "skewness 2", skewness(object@estimphi[2, 
                ])), 6, 2, byrow = TRUE))
            
        }
    }
    
    if (object@bic == 0) {
        
        if (dim(object@gridf)[1] == 1) {
            
            if (sum(object@cutoff) != 0) {
                print(matrix(c("kurtosis", kurtosis(object@estimphi[1, ]), "skewness", skewness(object@estimphi[1])), 2, 2, byrow = TRUE))
                print(matrix(c("sigma", sqrt(object@sigma2), "number of truncated values", length(object@cutoff) - sum(object@cutoff), 
                  "empiric mean", mean(object@estimphi[1, ]), "empiric sd", sd(object@estimphi[1, ])), 4, 2, byrow = TRUE))
            }
            if (sum(object@cutoff) == 0) {
                print(matrix(c("kurtosis", kurtosis(object@estimphi[1, ]), "skewness", skewness(object@estimphi[1, ])), 2, 2, byrow = TRUE))
                print(matrix(c("sigma", sqrt(object@sigma2), "empiric mean", mean(object@estimphi), "empiric sd", sd(object@estimphi)), 
                  3, 2, byrow = TRUE))
            }
        }
        
        if (dim(object@gridf)[1] == 2) {
            if (sum(object@cutoff) != 0) {
                print(matrix(c("sigma", sqrt(object@sigma2), "number of truncated values", length(object@cutoff) - sum(object@cutoff)), 
                  2, 2, byrow = TRUE))

                cat("\nRandom effects:\n")
                reff1 <- matrix(c(mean(object@estimphi[1, ]), sd(object@estimphi[1, ]), kurtosis(object@estimphi[1, ]),
                    skewness(object@estimphi[1, ])), 4, 1, byrow = TRUE)
                rownames(reff1) <- c("empiric mean 1", "empiric sd 1", "kurtosis 1", "skewness 1")
                print(reff1, quote = FALSE, right = TRUE)


                print(matrix(c("empiric mean 2", mean(object@estimphi[2, ]), "empiric sd 2", sd(object@estimphi[2, ]), "kurtosis 2", kurtosis(object@estimphi[2, 
                  ]), "skewness 2", skewness(object@estimphi[2, ])), 4, 2, byrow = TRUE), quote = FALSE, right = TRUE)
                return(invisible(list(reff1)))
            }
            if (sum(object@cutoff) == 0) {
                print(matrix(c("sigma", sqrt(object@sigma2)), 1, 2, byrow = TRUE))
                print(matrix(c("empiric mean 1", mean(object@estimphi[1, ]), "empiric sd 1", sd(object@estimphi[1, ]), "kurtosis 1", kurtosis(object@estimphi[1, 
                  ]), "skewness 1", skewness(object@estimphi[1, ])), 4, 2, byrow = TRUE))
                print(matrix(c("empiric mean 2", mean(object@estimphi[2, ]), "empiric sd 2", sd(object@estimphi[2, ]), "kurtosis 2", kurtosis(object@estimphi[1, 
                  ]), "skewness 2", skewness(object@estimphi[2, ])), 4, 2, byrow = TRUE))
            }
        }
    }
})



#' Short summary of the results of class object Bayes.fit
#' @description Method for the S4 class Bayes.fit
#' @param object Bayes.fit class
#' @param level default is 0.05
#' @param burnIn optional
#' @param thinning optional
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#'
setMethod("summary", signature = "Bayes.fit", definition = function(object, level = 0.05, burnIn, thinning) {
    if (missing(burnIn)) 
        burnIn <- object@burnIn
    if (missing(thinning)) 
        thinning <- object@thinning
    ind.samples <- seq(burnIn + 1, length(object@sigma2), by = thinning)
    if (length(object@random) == 2) {
        out <- list(sigma2.mean = mean(object@sigma2[ind.samples]), sigma2.cred_int = quantile(object@sigma2[ind.samples], c(level/2, 
            1 - level/2)), mu.mean = apply(object@mu[ind.samples, ], 2, mean), mu.cred_int = apply(object@mu[ind.samples, ], 2, quantile, 
            c(level/2, 1 - level/2)), omega.mean = apply(object@omega[ind.samples, ], 2, mean), omega.cred_int = apply(object@omega[ind.samples, 
            ], 2, quantile, c(level/2, 1 - level/2)), alpha.mean = apply(object@alpha[ind.samples, ], 2, mean), alpha.cred_int = apply(object@alpha[ind.samples, 
            ], 2, quantile, c(level/2, 1 - level/2)), beta.mean = apply(object@beta[ind.samples, ], 2, mean), beta.cred_int = apply(object@beta[ind.samples, 
            ], 2, quantile, c(level/2, 1 - level/2)))
    } else {
        if (object@random == 1) {
            out <- list(sigma2.mean = mean(object@sigma2[ind.samples]), sigma2.cred_int = quantile(object@sigma2[ind.samples], c(level/2, 
                1 - level/2)), mu.mean = mean(object@mu[ind.samples]), mu.cred_int = quantile(object@mu[ind.samples], c(level/2, 1 - level/2)), 
                omega.mean = mean(object@omega[ind.samples]), omega.cred_int = quantile(object@omega[ind.samples], c(level/2, 1 - level/2)), 
                alpha.mean = apply(object@alpha[ind.samples, ], 2, mean), alpha.cred_int = apply(object@alpha[ind.samples, ], 2, quantile, 
                  c(level/2, 1 - level/2)), beta.mean = mean(object@beta[ind.samples]), beta.cred_int = quantile(object@beta[ind.samples], 
                  c(level/2, 1 - level/2)))
        } else {
            out <- list(sigma2.mean = mean(object@sigma2[ind.samples]), sigma2.cred_int = quantile(object@sigma2[ind.samples], c(level/2, 
                1 - level/2)), mu.mean = mean(object@mu[ind.samples]), mu.cred_int = quantile(object@mu[ind.samples], c(level/2, 1 - level/2)), 
                omega.mean = mean(object@omega[ind.samples]), omega.cred_int = quantile(object@omega[ind.samples], c(level/2, 1 - level/2)), 
                alpha.mean = mean(object@alpha[ind.samples]), alpha.cred_int = quantile(object@alpha[ind.samples], c(level/2, 1 - level/2)), 
                beta.mean = apply(object@beta[ind.samples, ], 2, mean), beta.cred_int = apply(object@beta[ind.samples, ], 2, quantile, 
                  c(level/2, 1 - level/2)))
        }
    }

    return(out)
})


########################################################### PRINT

#' Description of print
#' @description Method for the S4 class Freq.fit
#' @param x Freq.fit class
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#'
setMethod("print", "Freq.fit", function(x) {
    if (x@bic != 0) {
        
        print(c("number of used trajectories", length(x@index)))
    }
    
    if (x@bic == 0) {
        
        if (sum(x@cutoff) != 0) {
            print(matrix(c("number of used trajectories", length(x@index), "number of truncated values", length(x@cutoff) - sum(x@cutoff)), 
                2, 2, byrow = TRUE))
        }
        if (sum(x@cutoff) == 0) {
            print(c("number of used trajectories", length(x@index)))
        }
    }
})


#' Print of acceptance rates of the MH steps
#' @description Method for the S4 class Bayes.fit
#' @param x Bayes.fit class
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#'
setMethod("print", "Bayes.fit", function(x) {
    if (length(x@random) == 2) {
      alpha <- sapply(1:ncol(x@alpha), function(i) length(unique(x@alpha[, i])))
      beta <- sapply(1:ncol(x@alpha), function(i) length(unique(x@beta[, i])))
      cat("\nacceptance rate for phi:")
      print(summary(cbind(alpha, beta)))
        
    } else {
        if (x@random == 1) {
          alpha <- apply(x@alpha, 2, function(vec) length(unique(vec))/length(vec))
          beta <- length(unique(x@beta))/length(x@beta)  
          cat("\nacceptance rates for random effect:\n")
          print(summary(alpha))
          cat(c("\nacceptance rate for fixed effect:", beta))
        }
        if (x@random == 2) {
          alpha <- length(unique(x@alpha))/length(x@alpha)
          beta <- apply(x@beta, 2, function(vec) length(unique(vec))/length(vec))
          cat("\nacceptance rates for random effect:\n")
          print(summary(beta))
          cat(c("\nacceptance rate for fixed effect:", alpha))
        }
    }
    if (x@model == "CIR") {
      sigma2 <- length(unique(x@sigma2))/length(x@sigma2)
      cat(c("\nacceptance rate for sigma:", sigma2))
    }else{
      sigma2 <- 1
    }
  
  invisible(list(alpha = alpha, beta = beta, sigma2 = sigma2))
    
})

########################################################### PLOT
#' Plot method for the frequentist estimation class object
#' 
#' @description Plot method for the S4 class Freq.fit
#' @param x Freq.fit class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 

setMethod(f = "plot", signature = "Freq.fit", definition = function(x, newwindow = FALSE, ...) {
    if (newwindow) {
        x11(width = 14)
    }
    
    
    if (dim(x@gridf)[1] == 1) {
        
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
            cex.axis = 0.7)
        hist(x@estimphi, main = "Density of the random effect", freq = FALSE, xlab = "", ylab = "", xlim = c(min(x@estimphi) * 0.8, max(x@estimphi) * 
            1.2), ylim = c(0, max(x@estimf) * 1.5), breaks = 12)
        lines(x@gridf, x@estimf, col = "red")
        
        if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(x@gridf)) {
            lines(x@gridf, x@estimf.trunc, col = "red", lty = 2)
        }
        
        if (x@bic == 0) {
            qqplot(x@estimphi, rnorm(length(x@estimphi), mean(x@estimphi), sd(x@estimphi)), xlab = "Normal Quantiles", ylab = "Sample Quantiles", 
                pch = 18, xlim = c(min(x@estimphi) * 0.8, max(x@estimphi) * 1.2), ylim = c(min(x@estimphi) * 0.8, max(x@estimphi) * 1.2), 
                cex.lab = 1.1)
            abline(0, 1)
        }
        if (x@bic != 0) {
            qqplot(x@estimphi, rnorm(length(x@estimphi), x@mu, x@omega), xlab = "Normal Quantiles", ylab = "Sample Quantiles", pch = 18, 
                xlim = c(min(x@estimphi) * 0.8, max(x@estimphi) * 1.2), ylim = c(min(x@estimphi) * 0.8, max(x@estimphi) * 1.2), cex.lab = 1.1)
            abline(0, 1)
        }
        
    }
    
    if (dim(x@gridf)[1] == 2) {
        
        if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(x@gridf[1, ])^2) {
            op <- par(mfrow = c(2, 4), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            persp3D(x@gridf[1, ], x@gridf[2, ], x@estimf, main = "Estimator", theta = 45, phi = 25, expand = 0.75, colkey = FALSE, bty = "b2")
            persp3D(x@gridf[1, ], x@gridf[2, ], x@estimf.trunc, main = "Truncated Estimator", theta = 45, phi = 25, expand = 0.75, colkey = FALSE, 
                bty = "b2")
            
            gridf1 <- x@gridf[1, ]
            gridf2 <- x@gridf[2, ]
            marg1 <- ((max(gridf2) - min(gridf2))/length(gridf2)) * rowSums(x@estimf)
            marg2 <- ((max(gridf1) - min(gridf1))/length(gridf1)) * colSums(x@estimf)
            
            hist(x@estimphi[1, ], main = "", freq = FALSE, xlab = "", ylab = "", xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, 
                ]) * 1.2), ylim = c(0, max(marg1) * 1.5), breaks = 12)
            lines(gridf1, marg1, col = "red")
            if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(gridf1)^2) {
                marg1_trunc <- ((max(gridf2) - min(gridf2))/length(gridf2)) * rowSums(x@estimf.trunc)
                lines(gridf1, marg1_trunc, col = "red", lty = 2)
            }
            
            hist(x@estimphi[2, ], main = "", freq = FALSE, xlab = "", ylab = "", xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, 
                ]) * 1.2), ylim = c(0, max(marg2) * 1.5), breaks = 12)
            lines(gridf2, marg2, col = "red")
            if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(gridf1)^2) {
                marg2_trunc <- ((max(gridf1) - min(gridf1))/length(gridf1)) * colSums(x@estimf.trunc)
                lines(gridf2, marg2_trunc, col = "red", lty = 2)
            }
            
            
            if (x@bic == 0) {
                
                qqplot(x@estimphi[1, ], rnorm(length(x@estimphi[1, ]), mean(x@estimphi[1, ]), sd(x@estimphi[1, ])), xlab = "Normal Quantiles", 
                  ylab = "Sample Quantiles", pch = 18, xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, ]) * 1.2), ylim = c(min(x@estimphi[1, 
                    ]) * 0.8, max(x@estimphi[1, ]) * 1.2), cex.lab = 1.1)
                abline(0, 1)
                qqplot(x@estimphi[2, ], rnorm(length(x@estimphi[2, ]), mean(x@estimphi[2, ]), sd(x@estimphi[2, ])), xlab = "Normal Quantiles", 
                  ylab = "Sample Quantiles", pch = 18, xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, ]) * 1.2), ylim = c(min(x@estimphi[2, 
                    ]) * 0.8, max(x@estimphi[2, ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
            }
            
            if (x@bic != 0) {
                
                qqplot(x@estimphi[1, ], rnorm(length(x@estimphi[1, ]), x@mu[1], x@omega[1]), xlab = "Normal Quantiles", ylab = "Sample Quantiles", 
                  pch = 18, xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, ]) * 1.2), ylim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, 
                    ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
                
                qqplot(x@estimphi[2, ], rnorm(length(x@estimphi[2, ]), x@mu[2], x@omega[2]), xlab = "Normal Quantiles", ylab = "Sample Quantiles", 
                  pch = 18, xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, ]) * 1.2), ylim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, 
                    ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
            }
        } else {
            
            op <- par(mfrow = c(2, 3), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            persp3D(x@gridf[1, ], x@gridf[2, ], x@estimf, main = "Estimator", theta = 45, phi = 25, expand = 0.75, colkey = FALSE, bty = "b2")
            gridf1 <- x@gridf[1, ]
            gridf2 <- x@gridf[2, ]
            marg1 <- ((max(gridf2) - min(gridf2))/length(gridf2)) * rowSums(x@estimf)
            marg2 <- ((max(gridf1) - min(gridf1))/length(gridf1)) * colSums(x@estimf)
            
            hist(x@estimphi[1, ], main = "", freq = FALSE, xlab = "", ylab = "", xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, 
                ]) * 1.2), ylim = c(0, max(marg1) * 1.5), breaks = 12)
            lines(gridf1, marg1, col = "red")
            if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(gridf1)^2) {
                marg1_trunc <- ((max(gridf2) - min(gridf2))/length(gridf2)) * rowSums(x@estimf.trunc)
                lines(gridf1, marg1_trunc, col = "red", lty = 2)
            }
            
            hist(x@estimphi[2, ], main = "", freq = FALSE, xlab = "", ylab = "", xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, 
                ]) * 1.2), ylim = c(0, max(marg2) * 1.5), breaks = 12)
            lines(gridf2, marg2, col = "red")
            if (sum(x@estimf.trunc) != 0 & sum(x@estimf.trunc == x@estimf) != length(gridf1)^2) {
                marg2_trunc <- ((max(gridf1) - min(gridf1))/length(gridf1)) * colSums(x@estimf.trunc)
                lines(gridf2, marg2_trunc, col = "red", lty = 2)
            }
            
            
            if (x@bic == 0) {
                
                qqplot(x@estimphi[1, ], rnorm(length(x@estimphi[1, ]), mean(x@estimphi[1, ]), sd(x@estimphi[1, ])), xlab = "Normal Quantiles", 
                  ylab = "Sample Quantiles", pch = 18, xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, ]) * 1.2), ylim = c(min(x@estimphi[1, 
                    ]) * 0.8, max(x@estimphi[1, ]) * 1.2), cex.lab = 1.1)
                abline(0, 1)
                qqplot(x@estimphi[2, ], rnorm(length(x@estimphi[2, ]), mean(x@estimphi[2, ]), sd(x@estimphi[2, ])), xlab = "Normal Quantiles", 
                  ylab = "Sample Quantiles", pch = 18, xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, ]) * 1.2), ylim = c(min(x@estimphi[2, 
                    ]) * 0.8, max(x@estimphi[2, ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
            }
            
            if (x@bic != 0) {
                
                qqplot(x@estimphi[1, ], rnorm(length(x@estimphi[1, ]), x@mu[1], x@omega[1]), xlab = "Normal Quantiles", ylab = "Sample Quantiles", 
                  pch = 18, xlim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, ]) * 1.2), ylim = c(min(x@estimphi[1, ]) * 0.8, max(x@estimphi[1, 
                    ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
                
                qqplot(x@estimphi[2, ], rnorm(length(x@estimphi[2, ]), x@mu[2], x@omega[2]), xlab = "Normal Quantiles", ylab = "Sample Quantiles", 
                  pch = 18, xlim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, ]) * 1.2), ylim = c(min(x@estimphi[2, ]) * 0.8, max(x@estimphi[2, 
                    ])) * 1.2, cex.lab = 1.1)
                abline(0, 1)
            }
        }
        
    }
    
    
})

######## 
#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x Bayes.fit class
#' @param plot.priorMean logical(1), if TRUE, prior means are added to the plots
#' @param reduced logical(1), if TRUE, the chains are reduced with the burn-in and thin rate
#' @param style one out of 'chains', 'acf', 'density' or 'cred.int'
#' @param level alpha for the credibility intervals, only for style 'cred.int', default = 0.05
#' @param true.phi only for style 'cred.int', for the case of known true values, e.g. for simulation
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
#' @examples
#' random <- c(1,2); sigma <- 0.1; param <- c(3, 0.5, 5, 0.2)
#' 
#' sim <- mixedsde.sim(M = 20, T = 1, N = 50, model = 'OU', random = random, 
#'        density.phi = 'normalnormal', param = param, sigma = sigma, X0 = 0, op.plot = 1)
#' 
#' # here: only 100 iterations for example - should be much more!
#' prior <- list(m = param[c(1,3)], v = param[c(1,3)], alpha.omega = c(11,11), 
#'            beta.omega = param[c(2,4)]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
#' estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random = random,
#'                 estim.method = 'paramBayes', prior = prior, nMCMC = 100) 
#' plot(estim_Bayes)
#' plot(estim_Bayes, style = 'cred.int', true.phi = sim$phi)
#' plot(estim_Bayes, style = 'acf')
#' plot(estim_Bayes, style = 'density')
#' 

setMethod(f = "plot", signature = "Bayes.fit", definition = function(x, plot.priorMean = FALSE, reduced = FALSE, style = c("chains", "acf", 
    "density", "cred.int"), level = 0.05, true.phi, newwindow = FALSE, ...) {
    if (newwindow) {
        x11(width = 10)
    }
    original.settings <- par(no.readonly = TRUE)
    style <- match.arg(style)
    ind <- 1:length(x@sigma2)
    if (reduced) 
        ind <- seq(x@burnIn + 1, length(x@sigma2), by = x@thinning)
    
    if (length(x@random) == 2) {
        op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.9, cex.lab = 0.8, 
            cex.axis = 0.8)
        
        if (style == "chains") {
            if (plot.priorMean) {
                
                # layout(matrix(c(1,2,3,4,5,6), ncol=2, byrow=TRUE), heights=c(4, 1, 4))
                
                plot(x@mu[ind, 1], ylim = range(c(x@mu[ind, 1], x@prior$m[1])), main = "Markov Chain", ylab = expression(mu[1]), type = "l")
                abline(h = x@prior$m[1], col = 2)  #if parametrisation of the package
                
                plot(x@mu[ind, 2], ylim = range(c(x@mu[ind, 2], x@prior$m[2])), main = "Markov Chain", ylab = expression(mu[2]), type = "l")
                abline(h = x@prior$m[2], col = 2)
                
                plot(x@omega[ind, 1], ylim = range(c(x@prior$beta.omega[1]/(x@prior$alpha.omega[1] - 1), x@omega[ind, 1])), main = "Markov Chain", 
                  ylab = expression(omega[1]^2), type = "l")
                abline(h = x@prior$beta.omega[1]/(x@prior$alpha.omega[1] - 1), col = 2)
                
                plot(x@omega[ind, 2], ylim = range(c(x@prior$beta.omega[2]/(x@prior$alpha.omega[2] - 1), x@omega[ind, 2])), main = "Markov Chain", 
                  ylab = expression(omega[2]^2), type = "l")
                abline(h = x@prior$beta.omega[2]/(x@prior$alpha.omega[2] - 1), col = 2)
                
                plot(x@sigma2[ind], ylim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = "Markov Chain", 
                  ylab = expression(sigma^2), type = "l")
                abline(h = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                
                par(mai = c(0, 0, 0, 0))
                plot.new()
                legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                
                
            } else {
                plot(x@mu[ind, 1], main = "Markov Chain", ylab = expression(mu[1]), type = "l")
                plot(x@mu[ind, 2], main = "Markov Chain", ylab = expression(mu[2]), type = "l")
                plot(x@omega[ind, 1], main = "Markov Chain", ylab = expression(omega[1]^2), type = "l")
                plot(x@omega[ind, 2], main = "Markov Chain", ylab = expression(omega[2]^2), type = "l")
                plot(x@sigma2[ind], main = "Markov Chain", ylab = expression(sigma^2), type = "l")
                
            }
            
        }
        if (style == "acf") {
            he <- acf(x@mu[ind, 1], plot = FALSE)
            plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ mu[1]), xlab = "lag", ylab = "acf", ylim = c(0, 
                1))
            he <- acf(x@mu[ind, 2], plot = FALSE)
            plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ mu[2]), xlab = "lag", ylab = "acf", ylim = c(0, 
                1))
            he <- acf(x@omega[ind, 1], plot = FALSE)
            plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ omega[1]^2), xlab = "lag", ylab = "acf", ylim = c(0, 
                1))
            he <- acf(x@omega[ind, 2], plot = FALSE)
            plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ omega[2]^2), xlab = "lag", ylab = "acf", ylim = c(0, 
                1))
            he <- acf(x@sigma2[ind], plot = FALSE)
            plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ sigma^2), xlab = "lag", ylab = "acf", ylim = c(0, 
                1))
            
        }
        if (style == "density") {
            if (plot.priorMean) {
                
                plot(density(x@mu[ind, 1]), xlim = range(c(x@mu[ind, 1], x@prior$m[1])), main = expression("Posterior of " ~ mu[1]), xlab = expression(mu[1]))
                abline(v = x@prior$m[1], col = 2)
                
                plot(density(x@mu[ind, 2]), xlim = range(c(x@mu[ind, 2], x@prior$m[2])), main = expression("Posterior of " ~ mu[2]), xlab = expression(mu[2]))
                abline(v = x@prior$m[2], col = 2)
                
                plot(density(x@omega[ind, 1]), xlim = range(c(x@omega[ind, 1], x@prior$beta.omega[1]/(x@prior$alpha.omega[1] - 1))), main = expression("Posterior of " ~ 
                  omega[1]^2), xlab = expression(omega[1]^2))
                abline(v = x@prior$beta.omega[1]/(x@prior$alpha.omega[1] - 1), col = 2)
                
                plot(density(x@omega[ind, 2]), xlim = range(c(x@omega[ind, 2], x@prior$beta.omega[2]/(x@prior$alpha.omega[2] - 1))), main = expression("Posterior of " ~ 
                  omega[2]^2), xlab = expression(omega[2]^2))
                abline(v = x@prior$beta.omega[2]/(x@prior$alpha.omega[2] - 1), col = 2)
                
                plot(density(x@sigma2[ind]), xlim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = expression("Posterior of " ~ 
                  sigma^2), xlab = expression(sigma^2))
                abline(v = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                
                par(mai = c(0, 0, 0, 0))
                plot.new()
                legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                
                
            } else {
                plot(density(x@mu[ind, 1]), main = expression("Posterior of " ~ mu[1]), xlab = expression(mu[1]))
                
                plot(density(x@mu[ind, 2]), main = expression("Posterior of " ~ mu[2]), xlab = expression(mu[2]))
                
                plot(density(x@omega[ind, 1]), main = expression("Posterior of " ~ omega[1]^2), xlab = expression(omega[1]^2))
                
                plot(density(x@omega[ind, 2]), main = expression("Posterior of " ~ omega[2]^2), xlab = expression(omega[2]^2))
                
                plot(density(x@sigma2[ind]), main = expression("Posterior of " ~ sigma^2), xlab = expression(sigma^2))
                
            }
            
            
            
        }
        if (style == "cred.int") {
            layout(matrix(c(1, 2, 3)), heights = c(4, 4, 1))
            
            par(mai = rep(0.5, 4))
            
            alpha.mean <- apply(x@alpha[ind, ], 2, mean)
            alpha.cred_int <- apply(x@alpha[ind, ], 2, quantile, c(level/2, 1 - level/2))
            
            
            if (!missing(true.phi)) {
                ra.plot <- range(c(range(alpha.cred_int), true.phi[1, -x@ind.4.prior]))
            } else {
                ra.plot <- range(alpha.cred_int)
            }
            plot(alpha.mean, ylim = ra.plot, ylab = expression(alpha), main = "credibility intervals", pch = 20)
            segments(1:length(alpha.mean), alpha.cred_int[1, ], 1:length(alpha.mean), alpha.cred_int[2, ])
            segments(1:length(alpha.mean) - 0.2, alpha.cred_int[1, ], 1:length(alpha.mean) + 0.2, alpha.cred_int[1, ])
            segments(1:length(alpha.mean) - 0.2, alpha.cred_int[2, ], 1:length(alpha.mean) + 0.2, alpha.cred_int[2, ])
            
            if (!missing(true.phi)) 
                points(true.phi[1, -x@ind.4.prior], col = 2, pch = 20)
            
            beta.mean <- apply(x@beta[ind, ], 2, mean)
            beta.cred_int <- apply(x@beta[ind, ], 2, quantile, c(level/2, 1 - level/2))
            
            if (!missing(true.phi)) {
                ra.plot <- range(c(range(beta.cred_int), true.phi[2, -x@ind.4.prior]))
            } else {
                ra.plot <- range(beta.cred_int)
            }
            plot(beta.mean, ylim = ra.plot, ylab = expression(beta), main = "credibility intervals", pch = 20)
            segments(1:length(beta.mean), beta.cred_int[1, ], 1:length(beta.mean), beta.cred_int[2, ])
            segments(1:length(beta.mean) - 0.2, beta.cred_int[1, ], 1:length(beta.mean) + 0.2, beta.cred_int[1, ])
            segments(1:length(beta.mean) - 0.2, beta.cred_int[2, ], 1:length(beta.mean) + 0.2, beta.cred_int[2, ])
            
            if (!missing(true.phi)) 
                points(true.phi[2, -x@ind.4.prior], col = 2, pch = 20)
            
            
            par(mai = c(0, 0, 0, 0))
            plot.new()
            if (!missing(true.phi)) {
                legend(x = "center", legend = c("credibility interval", "posterior mean", "true values"), col = c(1, 1, 2), lty = c(1, 
                  -1, -1), pch = c(-1, 20, 20), cex = 0.8, horiz = TRUE)
                
            } else {
                legend(x = "center", legend = c("credibility interval", "posterior mean"), lty = c(1, -1), pch = c(-1, 20), horiz = TRUE)
            }
            
        }
    } else {
        if (x@random == 1) {
            op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            if (style == "chains") {
                if (plot.priorMean) {
                  
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  
                  par(mai = rep(0.5, 4))
                  plot(x@mu[ind], ylim = range(c(x@mu[ind], x@prior$m[1])), main = "Markov Chain", ylab = expression(mu[1]), type = "l")
                  abline(h = x@prior$m[1], col = 2)  #if parametrisation of the package
                  
                  plot(x@beta[ind], ylim = range(c(x@beta[ind], x@prior$m[2])), main = "Markov Chain", ylab = expression(beta), type = "l")
                  abline(h = x@prior$m[2], col = 2)
                  
                  plot(x@omega[ind], ylim = range(c(x@prior$beta.omega/(x@prior$alpha.omega - 1), x@omega[ind])), main = "Markov Chain", 
                    ylab = expression(omega[1]^2), type = "l")
                  abline(h = x@prior$beta.omega/(x@prior$alpha.omega - 1), col = 2)
                  
                  plot(x@sigma2[ind], ylim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = "Markov Chain", 
                    ylab = expression(sigma^2), type = "l")
                  abline(h = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                  
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                  
                } else {
                  plot(x@mu[ind], main = "Markov Chain", ylab = expression(mu[1]), type = "l")
                  plot(x@beta[ind], main = "Markov Chain", ylab = expression(beta), type = "l")
                  plot(x@omega[ind], main = "Markov Chain", ylab = expression(omega[1]^2), type = "l")
                  plot(x@sigma2[ind], main = "Markov Chain", ylab = expression(sigma^2), type = "l")
                  
                }
                
            }
            if (style == "acf") {
                he <- acf(x@mu[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ mu[1]), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                he <- acf(x@beta[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ beta), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                he <- acf(x@omega[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ omega[1]^2), xlab = "lag", ylab = "acf", 
                  ylim = c(0, 1))
                he <- acf(x@sigma2[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ sigma^2), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                
            }
            if (style == "density") {
                if (plot.priorMean) {
                  
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  
                  par(mai = rep(0.5, 4))
                  
                  plot(density(x@mu[ind]), xlim = range(c(x@mu[ind], x@prior$m[1])), main = expression("Posterior of " ~ mu[1]), xlab = expression(mu[1]))
                  abline(v = x@prior$m[1], col = 2)
                  
                  plot(density(x@beta[ind]), xlim = range(c(x@beta[ind], x@prior$m[2])), main = expression("Posterior of " ~ beta), xlab = expression(beta))
                  abline(v = x@prior$m[2], col = 2)
                  
                  plot(density(x@omega[ind]), xlim = range(c(x@omega[ind], x@prior$beta.omega/(x@prior$alpha.omega - 1))), main = expression("Posterior of " ~ 
                    omega[1]^2), xlab = expression(omega[1]^2))
                  abline(v = x@prior$beta.omega/(x@prior$alpha.omega - 1), col = 2)
                  
                  plot(density(x@sigma2[ind]), xlim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = expression("Posterior of " ~ 
                    sigma^2), xlab = expression(sigma^2))
                  abline(v = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                  
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                  
                  
                } else {
                  plot(density(x@mu[ind]), main = expression("Posterior of " ~ mu[1]), xlab = expression(mu[1]))
                  
                  plot(density(x@beta[ind]), main = expression("Posterior of " ~ beta), xlab = expression(beta))
                  
                  plot(density(x@omega[ind]), main = expression("Posterior of " ~ omega[1]^2), xlab = expression(omega[1]^2))
                  
                  plot(density(x@sigma2[ind]), main = expression("Posterior of " ~ sigma^2), xlab = expression(sigma^2))
                  
                }
                
                
                
            }
            if (style == "cred.int") {
                layout(matrix(c(1, 2)), heights = c(5, 1))
                
                # par(mai=rep(0.5, 4))
                
                alpha.mean <- apply(x@alpha[ind, ], 2, mean)
                alpha.cred_int <- apply(x@alpha[ind, ], 2, quantile, c(level/2, 1 - level/2))
                
                if (!missing(true.phi)) {
                  ra.plot <- range(c(range(alpha.cred_int), true.phi[-x@ind.4.prior]))
                } else {
                  ra.plot <- range(alpha.cred_int)
                }
                plot(alpha.mean, ylim = ra.plot, ylab = expression(alpha), main = "credibility intervals", pch = 20)
                segments(1:length(alpha.mean), alpha.cred_int[1, ], 1:length(alpha.mean), alpha.cred_int[2, ])
                segments(1:length(alpha.mean) - 0.2, alpha.cred_int[1, ], 1:length(alpha.mean) + 0.2, alpha.cred_int[1, ])
                segments(1:length(alpha.mean) - 0.2, alpha.cred_int[2, ], 1:length(alpha.mean) + 0.2, alpha.cred_int[2, ])
                
                if (!missing(true.phi)) 
                  points(true.phi[-x@ind.4.prior], col = 2, pch = 20)
                
                par(mai = c(0, 0, 0, 0))
                plot.new()
                if (!missing(true.phi)) {
                  legend(x = "center", legend = c("credibility interval", "posterior mean", "true values"), col = c(1, 1, 2), lty = c(1, 
                    -1, -1), pch = c(-1, 20, 20), cex = 0.8, horiz = TRUE)
                  
                } else {
                  legend(x = "center", legend = c("credibility interval", "posterior mean"), lty = c(1, -1), pch = c(-1, 20), horiz = TRUE)
                }
                
            }
            
        } else {
            op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            if (style == "chains") {
                if (plot.priorMean) {
                  
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  
                  par(mai = rep(0.5, 4))
                  
                  plot(x@alpha[ind], ylim = range(c(x@alpha[ind], x@prior$m[1])), main = "Markov Chain", ylab = expression(alpha), type = "l")
                  abline(h = x@prior$m[1], col = 2)
                  
                  plot(x@mu[ind], ylim = range(c(x@mu[ind], x@prior$m[2])), main = "Markov Chain", ylab = expression(mu[2]), type = "l")
                  abline(h = x@prior$m[2], col = 2)
                  
                  plot(x@omega[ind], ylim = range(c(x@prior$beta.omega/(x@prior$alpha.omega - 1), x@omega[ind])), main = "Markov Chain", 
                    ylab = expression(omega[2]^2), type = "l")
                  abline(h = x@prior$beta.omega/(x@prior$alpha.omega - 1), col = 2)
                  
                  plot(x@sigma2[ind], ylim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = "Markov Chain", 
                    ylab = expression(sigma^2), type = "l")
                  abline(h = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                  
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                  
                } else {
                  plot(x@alpha[ind], main = "Markov Chain", ylab = expression(alpha), type = "l")
                  plot(x@mu[ind], main = "Markov Chain", ylab = expression(mu[2]), type = "l")
                  plot(x@omega[ind], main = "Markov Chain", ylab = expression(omega[2]^2), type = "l")
                  plot(x@sigma2[ind], main = "Markov Chain", ylab = expression(sigma^2), type = "l")
                }
            }
            if (style == "acf") {
                he <- acf(x@alpha[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ alpha), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                he <- acf(x@mu[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ mu[2]), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                he <- acf(x@omega[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ omega[2]^2), xlab = "lag", ylab = "acf", 
                  ylim = c(0, 1))
                he <- acf(x@sigma2[ind], plot = FALSE)
                plot(he$lag, he$acf, type = "h", main = expression("Acf of simulations for " ~ sigma^2), xlab = "lag", ylab = "acf", ylim = c(0, 
                  1))
                
            }
            if (style == "density") {
                if (plot.priorMean) {
                  
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  
                  par(mai = rep(0.5, 4))
                  
                  plot(density(x@alpha[ind]), xlim = range(c(x@alpha[ind], x@prior$m[1])), main = expression("Posterior of " ~ alpha), 
                    xlab = expression(alpha))
                  abline(v = x@prior$m[1], col = 2)
                  
                  plot(density(x@mu[ind]), xlim = range(c(x@mu[ind], x@prior$m[2])), main = expression("Posterior of " ~ mu[2]), xlab = expression(mu[2]))
                  abline(v = x@prior$m[2], col = 2)
                  
                  plot(density(x@omega[ind]), xlim = range(c(x@omega[ind], x@prior$beta.omega/(x@prior$alpha.omega - 1))), main = expression("Posterior of " ~ 
                    omega[2]^2), xlab = expression(omega[2]^2))
                  abline(v = x@prior$beta.omega/(x@prior$alpha.omega - 1), col = 2)
                  
                  plot(density(x@sigma2[ind]), xlim = range(c(x@sigma2[ind], x@prior$beta.sigma/(x@prior$alpha.sigma - 1))), main = expression("Posterior of " ~ 
                    sigma^2), xlab = expression(sigma^2))
                  abline(v = x@prior$beta.sigma/(x@prior$alpha.sigma - 1), col = 2)
                  
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "prior mean", col = 2, lty = 1)
                  
                } else {
                  plot(density(x@alpha[ind]), main = expression("Posterior of " ~ alpha), xlab = expression(alpha))
                  
                  plot(density(x@mu[ind]), main = expression("Posterior of " ~ mu[2]), xlab = expression(mu[2]))
                  
                  plot(density(x@omega[ind]), main = expression("Posterior of " ~ omega[2]^2), xlab = expression(omega[2]^2))
                  
                  plot(density(x@sigma2[ind]), main = expression("Posterior of " ~ sigma^2), xlab = expression(sigma^2))
                  
                }
                
                
                
            }
            if (style == "cred.int") {
                layout(matrix(c(1, 2)), heights = c(5, 1))
                
                # par(mai=rep(0.5, 4))
                beta.mean <- apply(x@beta[ind, ], 2, mean)
                beta.cred_int <- apply(x@beta[ind, ], 2, quantile, c(level/2, 1 - level/2))
                
                
                if (!missing(true.phi)) {
                  ra.plot <- range(c(range(beta.cred_int), true.phi[-x@ind.4.prior]))
                } else {
                  ra.plot <- range(beta.cred_int)
                }
                plot(beta.mean, ylim = ra.plot, ylab = expression(beta), main = "credibility intervals", pch = 20)
                segments(1:length(beta.mean), beta.cred_int[1, ], 1:length(beta.mean), beta.cred_int[2, ])
                segments(1:length(beta.mean) - 0.2, beta.cred_int[1, ], 1:length(beta.mean) + 0.2, beta.cred_int[1, ])
                segments(1:length(beta.mean) - 0.2, beta.cred_int[2, ], 1:length(beta.mean) + 0.2, beta.cred_int[2, ])
                
                if (!missing(true.phi)) 
                  points(true.phi[-x@ind.4.prior], col = 2, pch = 20)
                
                par(mai = c(0, 0, 0, 0))
                plot.new()
                if (!missing(true.phi)) {
                  legend(x = "center", legend = c("credibility interval", "posterior mean", "true values"), col = c(1, 1, 2), lty = c(1, 
                    -1, -1), pch = c(-1, 20, 20), cex = 0.8, horiz = TRUE)
                  
                } else {
                  legend(x = "center", legend = c("credibility interval", "posterior mean"), lty = c(1, -1), pch = c(-1, 20), horiz = TRUE)
                }
                
            }
        }
        
    }
    # set plot settings back
    par(original.settings)
})


######## 
#' Plot method for the Bayesian prediction class object
#' 
#' @description Plot method for the S4 class Bayes.pred
#' @param x Bayes.fit class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.legend logical(1)
#' @param ylim optional
#' @param xlab optional, default 'times'
#' @param ylab optional, default 'X'
#' @param col color for the prediction intervals, default 3
#' @param lwd linewidth for the prediction intervals, default 2 
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setMethod(f = "plot", signature = "Bayes.pred", definition = function(x, newwindow = FALSE, plot.legend = TRUE, ylim, xlab = "times", 
    ylab = "X", col = 3, lwd = 2, ...) {
    if (newwindow) {
        x11(width = 10)
    }
    original.settings <- par(no.readonly = TRUE)
    estim <- x@estim
    
    if (length(x@qu.u) == 0) {
        op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
            cex.axis = 0.7)
        plot(estim$times, x@Xpred[1, ], type = "l", ylim = range(c(range(x@Xpred), range(estim$X))), xlab = xlab, ylab = ylab, col = col, 
            lwd = lwd)
        for (i in 2:nrow(x@Xpred[i, ])) lines(estim$times, x@Xpred[1, ], col = col, lwd = lwd)
        for (i in 1:nrow(estim$X)) lines(estim$times, estim$X[i, ])
        
        if (plot.legend) 
            legend("bottomright", c("data", "drawn trajectories"), lty = 1, col = c(1, col), lwd = c(1, lwd), box.lty = 0, inset = 0.01)
        
    } else {
        qu.l <- x@qu.l
        qu.u <- x@qu.u
        cr <- x@coverage.rate
        
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
            cex.axis = 0.7)
        if (missing(ylim)) 
            ylim <- range(c(min(qu.l), max(qu.u), range(estim$X)))
        
        plot(estim$times[-1], qu.l, type = "l", ylim = ylim, xlab = xlab, ylab = ylab, col = col, lwd = lwd, ...)
        lines(estim$times[-1], qu.u, col = col, lwd = lwd, ...)
        
        for (i in 1:nrow(estim$X)) lines(estim$times, estim$X[i, ])
        lines(estim$times[-1], qu.l, col = col, lwd = lwd, ...)
        lines(estim$times[-1], qu.u, col = col, lwd = lwd, ...)
        
        if (plot.legend) 
            legend("bottomright", c("data", "prediction intervals"), lty = 1, col = c(1, col), lwd = c(1, lwd), cex = 0.7, box.lty = 0, 
                inset = 0.01)
        
        plot(estim$times[-1], cr, ylim = c(min(cr) * 0.9, max(c(cr), 1)), type = "l", xlab = xlab, ylab = "coverage rates")
        abline(h = 0.95, col = 2, lty = 2)
        if (plot.legend) 
            legend("bottomright", "95%", lty = 2, col = 2, cex = 0.7, box.lty = 0, inset = 0.01)
    }
    
    # set plot settings back
    par(original.settings)
})

#' Comparing plot method
#' 
#' @description Method for classes
#' @param x Bayes.fit or Bayes.pred class
#' @param y Bayes.fit or Bayes.pred class
#' @param z Bayes.fit or Bayes.pred class (optional)
#' @param ... other parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setGeneric("plot2compare", function(x, y, z, ...) {
    standardGeneric("plot2compare")
})
######## 
#' Comparing plot method plot2compare for three Bayesian estimation class objects
#' 
#' @description Comparison of the posterior densities for up to three S4 class Bayes.fit objects
#' @param x Bayes.fit class
#' @param y Bayes.fit class
#' @param z Bayes.fit class (optional)
#' @param names character vector of names for x, y and z
#' @param true.values list of parameters to compare with the estimations, if available
#' @param reduced logical(1), if TRUE, the chains are reduced with the burn-in and thin rate
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' 
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setMethod(f = "plot2compare", signature = "Bayes.fit", definition = function(x, y, z, names, true.values, reduced = TRUE, newwindow = FALSE) {
    if (newwindow) {
        x11(width = 10)
    }
    if (!missing(y)) {
        list.classes <- list(x, y)
    }
    if (!missing(z)) {
        list.classes <- list(x, y, z)
    }
    original.settings <- par(no.readonly = TRUE)
    
    l.li <- length(list.classes)
    
    
    if (reduced) {
        ind <- lapply(1:l.li, function(i) seq(list.classes[[i]]@burnIn + 1, length(list.classes[[i]]@sigma2), by = list.classes[[i]]@thinning))
    } else {
        ind <- lapply(1:l.li, function(i) 1:length(list.classes[[i]]@sigma2))
    }
    
    random <- sapply(1:l.li, function(i) list.classes[[i]]@random)
    if (is.list(random)) 
        he <- FALSE
    if (is.matrix(random)) {
        he <- TRUE
        random <- c(1, 2)
    }
    if (is.numeric(random)) {
        if (length(unique(random)) == 1) {
            he <- TRUE
            random <- unique(random)
        } else {
            he <- FALSE
        }
    }
    if (!he) {
        print("comparison of parameters not possible, only variance is plotted")
        
        op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.9, cex.lab = 0.8, 
            cex.axis = 0.8)
        ra.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$x))
        ra.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$y))
        plot(density(list.classes[[1]]@sigma2[ind[[1]]]), xlim = ra.x, ylim = ra.y, main = expression("Posterior of " ~ sigma^2), xlab = expression(sigma^2))
        for (i in 2:l.li) lines(density(list.classes[[i]]@sigma2[ind[[i]]]), col = i)
        if (!missing(names)) 
            legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
    }
    if (he) {
        if (length(random) == 2) {
            
            op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.9, cex.lab = 0.8, 
                cex.axis = 0.8)
            
            if (!missing(true.values)) {
                ra.mu1.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 1])$x)), true.values$mu[1]))
                ra.mu2.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 2])$x)), true.values$mu[2]))
                ra.omega1.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 1])$x)), true.values$omega[1]))
                ra.omega2.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 2])$x)), true.values$omega[2]))
                ra.sigma2.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$x)), true.values$sigma2))
            } else {
                ra.mu1.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 1])$x))
                ra.mu2.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 2])$x))
                ra.omega1.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 1])$x))
                ra.omega2.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 2])$x))
                ra.sigma2.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$x))
            }
            ra.mu1.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 1])$y))
            ra.mu2.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]], 2])$y))
            ra.omega1.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 1])$y))
            ra.omega2.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]], 2])$y))
            ra.sigma2.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$y))
            
            
            plot(density(list.classes[[1]]@mu[ind[[1]], 1]), xlim = ra.mu1.x, ylim = ra.mu1.y, main = expression("Posterior of " ~ mu[1]), 
                xlab = expression(mu[1]))
            for (i in 2:l.li) lines(density(list.classes[[i]]@mu[ind[[i]], 1]), col = i)
            if (!missing(true.values)) 
                abline(v = true.values$mu[1], lty = 2)
            if (!missing(names)) 
                legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
            
            plot(density(list.classes[[1]]@mu[ind[[1]], 2]), xlim = ra.mu2.x, ylim = ra.mu2.y, main = expression("Posterior of " ~ mu[2]), 
                xlab = expression(mu[2]))
            for (i in 2:l.li) lines(density(list.classes[[i]]@mu[ind[[i]], 2]), col = i)
            if (!missing(true.values)) 
                abline(v = true.values$mu[2], lty = 2)
            if (!missing(names)) 
                legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
            
            plot(density(list.classes[[1]]@omega[ind[[1]], 1]), xlim = ra.omega1.x, ylim = ra.omega1.y, main = expression("Posterior of " ~ 
                omega[1]), xlab = expression(omega[1]))
            for (i in 2:l.li) lines(density(list.classes[[i]]@omega[ind[[i]], 1]), col = i)
            if (!missing(true.values)) 
                abline(v = true.values$omega[1], lty = 2)
            if (!missing(names)) 
                legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
            
            plot(density(list.classes[[1]]@omega[ind[[1]], 2]), xlim = ra.omega2.x, ylim = ra.omega2.y, main = expression("Posterior of " ~ 
                omega[2]), xlab = expression(omega[2]))
            for (i in 2:l.li) lines(density(list.classes[[i]]@omega[ind[[i]], 2]), col = i)
            if (!missing(true.values)) 
                abline(v = true.values$omega[2], lty = 2)
            if (!missing(names)) 
                legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
            
            plot(density(list.classes[[1]]@sigma2[ind[[1]]]), xlim = ra.sigma2.x, ylim = ra.sigma2.y, main = expression("Posterior of " ~ 
                sigma^2), xlab = expression(sigma^2))
            for (i in 2:l.li) lines(density(list.classes[[i]]@sigma2[ind[[i]]]), col = i)
            if (!missing(true.values)) 
                abline(v = true.values$sigma2, lty = 2)
            if (!missing(names)) 
                legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
            
            if (!missing(true.values)) {
                par(mai = c(0, 0, 0, 0))
                plot.new()
                legend(x = "center", legend = "true value", lty = 2)
            }
            
        } else {
            if (!missing(true.values)) {
                ra.mu.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]]])$x)), true.values$mu))
                ra.omega.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]]])$x)), true.values$omega))
                ra.sigma2.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$x)), true.values$sigma2))
            } else {
                ra.mu.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]]])$x))
                ra.omega.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]]])$x))
                ra.sigma2.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$x))
            }
            ra.mu.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@mu[ind[[i]]])$y))
            ra.omega.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@omega[ind[[i]]])$y))
            ra.sigma2.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@sigma2[ind[[i]]])$y))
            
            if (random == 1) {
                op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                  cex.axis = 0.7)
                
                if (!missing(true.values)) {
                  ra.beta.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@beta[ind[[i]]])$x)), true.values$beta))
                } else {
                  ra.beta.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@beta[ind[[i]]])$x))
                }
                ra.beta.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@beta[ind[[i]]])$y))
                
                if (!missing(true.values)) {
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  par(mai = rep(0.5, 4))
                }
                
                plot(density(list.classes[[1]]@mu[ind[[1]]]), xlim = ra.mu.x, ylim = ra.mu.y, main = expression("Posterior of " ~ mu[1]), 
                  xlab = expression(mu[1]))
                for (i in 2:l.li) lines(density(list.classes[[i]]@mu[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$mu, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@beta[ind[[1]]]), xlim = ra.beta.x, ylim = ra.beta.y, main = expression("Posterior of " ~ 
                  beta), xlab = expression(beta))
                for (i in 2:l.li) lines(density(list.classes[[i]]@beta[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$beta, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@omega[ind[[1]]]), xlim = ra.omega.x, ylim = ra.omega.y, main = expression("Posterior of " ~ 
                  omega[1]), xlab = expression(omega[1]))
                for (i in 2:l.li) lines(density(list.classes[[i]]@omega[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$omega, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@sigma2[ind[[1]]]), xlim = ra.sigma2.x, ylim = ra.sigma2.y, main = expression("Posterior of " ~ 
                  sigma^2), xlab = expression(sigma^2))
                for (i in 2:l.li) lines(density(list.classes[[i]]@sigma2[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$sigma2, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                if (!missing(true.values)) {
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "true value", lty = 2)
                }
                
            } else {
                op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                  cex.axis = 0.7)
                
                if (!missing(true.values)) {
                  ra.alpha.x <- range(c(range(sapply(1:l.li, function(i) density(list.classes[[i]]@alpha[ind[[i]]])$x)), true.values$alpha))
                } else {
                  ra.alpha.x <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@alpha[ind[[i]]])$x))
                }
                ra.alpha.y <- range(sapply(1:l.li, function(i) density(list.classes[[i]]@alpha[ind[[i]]])$y))
                
                if (!missing(true.values)) {
                  layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 2, byrow = TRUE), heights = c(4, 4, 1))
                  par(mai = rep(0.5, 4))
                }
                
                plot(density(list.classes[[1]]@alpha[ind[[1]]]), xlim = ra.alpha.x, ylim = ra.alpha.y, main = expression("Posterior of " ~ 
                  alpha), xlab = expression(alpha))
                for (i in 2:l.li) lines(density(list.classes[[i]]@alpha[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$alpha, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@mu[ind[[1]]]), xlim = ra.mu.x, ylim = ra.mu.y, main = expression("Posterior of " ~ mu[2]), 
                  xlab = expression(mu[2]))
                for (i in 2:l.li) lines(density(list.classes[[i]]@mu[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$mu, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@omega[ind[[1]]]), xlim = ra.omega.x, ylim = ra.omega.y, main = expression("Posterior of " ~ 
                  omega[2]), xlab = expression(omega[2]))
                for (i in 2:l.li) lines(density(list.classes[[i]]@omega[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$omega, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                plot(density(list.classes[[1]]@sigma2[ind[[1]]]), xlim = ra.sigma2.x, ylim = ra.sigma2.y, main = expression("Posterior of " ~ 
                  sigma^2), xlab = expression(sigma^2))
                for (i in 2:l.li) lines(density(list.classes[[i]]@sigma2[ind[[i]]]), col = i)
                if (!missing(true.values)) 
                  abline(v = true.values$sigma2, lty = 2)
                if (!missing(names)) 
                  legend("topleft", names, col = 1:l.li, lty = 1, cex = 0.7, box.lty = 0, inset = 0.001)
                
                if (!missing(true.values)) {
                  par(mai = c(0, 0, 0, 0))
                  plot.new()
                  legend(x = "center", legend = "true value", lty = 2)
                  
                }
            }
        }
    }
    
    # set plot settings back
    par(original.settings)
})

######## 
#' Comparing plot method plot2compare for three Bayesian prediction class objects
#' 
#' @description Comparison of the results for up to three S4 class Bayes.pred objects
#' @param x Bayes.pred class
#' @param y Bayes.pred class
#' @param z Bayes.pred class (optional)
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.legend logical(1), if TRUE, a legend is added
#' @param names character vector with names for the three objects appearing in the legend
#' @param ylim optional
#' @param xlab optional, default 'times'
#' @param ylab optional, default 'X'
#' @param ... optional plot parameters
#' 
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' @examples 
#' random <- 1; sigma <- 0.1; fixed <- 5; param <- c(3, 0.5)
#' sim <- mixedsde.sim(M = 20, T = 1, N = 50, model = 'OU', random = random, fixed = fixed,
#'        density.phi = 'normal',param= param, sigma= sigma, X0 = 0, op.plot = 1)
#' 
#' # here: only 100 iterations for example - should be much more!
#' estim_Bayes_withoutprior <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU',
#'              random, estim.method = 'paramBayes',  nMCMC = 100)
#' prior <- list( m = c(param[1], fixed), v = c(param[1], fixed), alpha.omega = 11,
#'             beta.omega = param[2]^2*10, alpha.sigma = 10, beta.sigma = sigma^2*9)
#' estim_Bayes <- mixedsde.fit(times = sim$times, X = sim$X, model = 'OU', random, 
#'            estim.method = 'paramBayes', prior = prior, nMCMC = 100) 
#' plot2compare(estim_Bayes, estim_Bayes_withoutprior, names = c('with prior', 'without prior'))

setMethod(f = "plot2compare", signature = "Bayes.pred", definition = function(x, y, z, newwindow = FALSE, plot.legend = TRUE, names, ylim, 
    xlab = "times", ylab = "X", ...) {
    if (newwindow) {
        x11(width = 10)
    }
    original.settings <- par(no.readonly = TRUE)
    
    if (!missing(y)) {
        list.classes <- list(x, y)
    }
    if (!missing(z)) {
        list.classes <- list(x, y, z)
    }
    
    l.li <- length(list.classes)
    
    X <- lapply(1:l.li, function(i) list.classes[[i]]@estim$X)
    times <- lapply(1:l.li, function(i) list.classes[[i]]@estim$times)
    
    if (any(sapply(1:l.li, function(i) length(list.classes[[i]]@qu.u)) == 0)) {
        op <- par(mfrow = c(1, l.li), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
            cex.axis = 0.7)
        if (missing(names)) {
            for (i in 1:l.li) {
                plot(times[[i]], list.classes[[i]]@Xpred[1, ], type = "l", ylim = range(list.classes[[i]]@Xpred), xlab = xlab, ylab = ylab, 
                  col = 2, lwd = 2, ...)
                for (j in 1:nrow(list.classes[[i]]@Xpred)) lines(times[[i]], list.classes[[i]]@Xpred[j, ], lwd = 2, col = 2)
                for (j in 1:nrow(X[[i]])) lines(times[[i]], X[[i]][j, ])
                if (plot.legend) 
                  legend("bottomright", c("data", "drawn trajectories"), lty = 1, col = c(1, 2), lwd = c(1, 2), cex = 0.7, box.lty = 0, 
                    inset = 0.01)
            }
        } else {
            for (i in 1:l.li) {
                plot(times[[i]], list.classes[[i]]@Xpred[1, ], type = "l", ylim = range(list.classes[[i]]@Xpred), main = names[i], xlab = xlab, 
                  ylab = ylab, col = 2, lwd = 2, ...)
                for (j in 1:nrow(list.classes[[i]]@Xpred)) lines(times[[i]], list.classes[[i]]@Xpred[j, ], lwd = 2, col = 2)
                for (j in 1:nrow(X[[i]])) lines(times[[i]], X[[i]][j, ])
                if (plot.legend) 
                  legend("bottomright", c("data", "drawn trajectories"), lty = 1, col = c(1, 2), lwd = c(1, 2), cex = 0.7, box.lty = 0, 
                    inset = 0.01)
            }
        }
    } else {
        qu.l <- lapply(1:l.li, function(i) list.classes[[i]]@qu.l)
        qu.u <- lapply(1:l.li, function(i) list.classes[[i]]@qu.u)
        cr <- lapply(1:l.li, function(i) list.classes[[i]]@coverage.rate)
        
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
            cex.axis = 0.7)
        if (missing(ylim)) 
            ylim <- c(min(unlist(qu.l)), max(unlist(qu.u)))
        
        plot(times[[1]][-1], qu.l[[1]], type = "l", ylim = ylim, xlab = xlab, ylab = ylab, lwd = 2, ...)
        lines(times[[1]][-1], qu.u[[1]], lwd = 2, ...)
        
        for (j in 1:l.li) for (i in 1:nrow(X[[1]])) lines(times[[j]], X[[j]][i, ], col = "grey")
        for (i in 2:l.li) lines(times[[i]][-1], qu.l[[i]], col = i, lwd = 2, ...)
        for (i in 2:l.li) lines(times[[i]][-1], qu.u[[i]], col = i, lwd = 2, ...)
        
        if (plot.legend) {
            if (!missing(names)) {
                legend("bottomright", c("data", names), lty = 1, col = c("grey", 1:l.li), lwd = c(1, rep(2, l.li)), cex = 0.7, box.lty = 0, 
                  inset = 0.01)
            } else {
                legend("bottomright", c("data", "prediction intervals"), lty = 1, col = c("grey", 1:l.li), lwd = c(1, rep(2, l.li)), cex = 0.7, 
                  box.lty = 0, inset = 0.01)
            }
        }
        
        plot(times[[1]][-1], cr[[1]], ylim = c(min(unlist(cr)) * 0.9, max(unlist(cr), 1)), lwd = 2, type = "l", xlab = xlab, ylab = "coverage rates")
        for (i in 2:l.li) lines(times[[i]][-1], cr[[i]], col = i, lwd = 2)
        abline(h = 0.95, lty = 2)
        if (plot.legend) {
            if (!missing(names)) {
                legend("bottomright", c("95%", names), lty = c(2, rep(1, l.li)), col = c(1, 1:l.li), lwd = c(1, rep(2, l.li)), cex = 0.7, 
                  box.lty = 0, inset = 0.01)
            } else {
                legend("bottomright", "95%", lty = 2, cex = 0.7, box.lty = 0, inset = 0.01)
            }
        }
        
    }
    # set plot settings back
    par(original.settings)
})

########################################################### VALIDATION

#' Validation of the chosen model.
#' 
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Freq.fit or Bayes.fit class
#' @param ... other optional parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setGeneric("valid", function(x, ...) {
    standardGeneric("valid")
})


#' Validation of the chosen model.
#' 
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Freq.fit class
#' @param Mrep number of trajectories to be drawn
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.valid logical(1), if TRUE, the results are depicted grafically
#' @param numj optional number of series to be validated
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 

setMethod(f = "valid", signature = "Freq.fit", definition = function(x, Mrep = 100, newwindow = FALSE, plot.valid = TRUE, 
    numj, ...) {
    
    if (newwindow) {
        x11(width = 10)
    }
    
    ## local sde.sim to sink undesired output away into a tempfile
    con <- file(tempfile(), open="w")
    on.exit(close(con))
    sde.sim <- function(...){
      sink(con)
      res <- sde::sde.sim(...)
      sink(NULL)
      res
    }
    
    Xtrue <- x@X
    times <- round(x@times, 10)
    Tend <- max(times)
    del <- round(min(diff(times)), 10)
    timessimu <- round(seq(del, Tend, by = del), 10)
    
    M <- dim(Xtrue)[1]
    sig <- sqrt(x@sigma2)
    
    if (missing(numj)) {
        
        if (dim(x@gridf)[1] == 2) {
            
            phihat <- x@estimphi
            
            Xnew <- as.list(1:M)  # for each phihat_j a new sample size Mrep
            
            if (x@model == "OU") {
                for (j in 1:M) {
                  Xnew[[j]] <- matrix(0, Mrep, N + 1)
                  Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "EA", theta = c(phihat[, 
                    j], sig), model = "OU", M = Mrep))
                }
            }
            if (x@model == "CIR") {
                for (j in 1:M) {
                  Xnew[[j]] <- matrix(0, Mrep, N + 1)
                  Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(phihat[, 
                    j], sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
                }
            }
            
            vecttimes <- intersect(round(timessimu, 10), round(times, 10))
            N <- length(vecttimes)
            
            q <- matrix(0, M, N)
            for (j in 1:M) {
                for (i in 1:N) {
                  q[j, i] <- sum(Xtrue[j, i] > Xnew[[j]][, which(timessimu == vecttimes[i])[1]])/Mrep
                }
            }
            if (plot.valid == 1) {
                op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                  cex.axis = 0.7)
                plotnumj <- floor(runif(1, 1, M + 1))
                
                plot(c(0, timessimu), Xnew[[plotnumj]][1, ], type = "l", ylim = c(min(Xnew[[plotnumj]]) * 0.8, max(Xnew[[plotnumj]]) * 
                  1.2), xlab = "", ylab = "")
                for (k in 1:Mrep) {
                  lines(c(0, timessimu), Xnew[[plotnumj]][k, ])
                }
                lines(times, Xtrue[plotnumj, ], col = "red", lwd = 2)
                
                plot(1:N/N, sort(q[plotnumj, ]), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
                abline(0, 1)
            }
        }
        if (dim(x@gridf)[1] == 1) {
            
            phihat <- x@estimphi
            if (x@estim.fix == 1) {
                paramfixed <- x@estim.fixed
            }
            if (x@estim.fix == 0) {
                paramfixed <- x@fixed
            }
            
            Xnew <- as.list(1:M)  # for each phihat_j a new sample size Mrep
            
            if (sum(x@random) == 1) {
                if (x@model == "OU") {
                  
                  for (j in 1:M) {
                    Xnew[[j]] <- matrix(0, Mrep, length(timessimu))
                    Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "EA", theta = c(phihat[j], 
                      paramfixed, sig), model = "OU", M = Mrep))
                  }
                }
                if (x@model == "CIR") {
                  
                  for (j in 1:M) {
                    Xnew[[j]] <- matrix(0, Mrep, length(timessimu))
                    Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(phihat[j], 
                      paramfixed, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
                  }
                }
            }
            
            if (sum(x@random) == 2) {
                
                if (x@model == "OU") {
                  
                  for (j in 1:M) {
                    Xnew[[j]] <- matrix(0, Mrep, length(timessimu))
                    Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "EA", theta = c(paramfixed, 
                      phihat[j], sig), model = "OU", M = Mrep))
                  }
                }
                if (x@model == "CIR") {
                  
                  for (j in 1:M) {
                    Xnew[[j]] <- matrix(0, Mrep, length(timessimu))
                    Xnew[[j]] <- t(sde.sim(T = Tend, X0 = Xtrue[j, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(paramfixed, 
                      phihat[j], sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
                    
                  }
                }
            }
            
            vecttimes <- intersect(round(timessimu, 10), round(times, 10))
            N <- length(vecttimes)
            
            q <- matrix(0, M, N)
            for (j in 1:M) {
                for (i in 1:N) {
                  q[j, i] <- sum(Xtrue[j, i] > Xnew[[j]][, which(timessimu == vecttimes[i])[1]])/Mrep
                }
            }
            
            if (plot.valid == 1) {
                op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                  cex.axis = 0.7)
                
                
                plotnumj <- floor(M * runif(1, 0, 1) + 1)
                plot(c(0, timessimu), Xnew[[plotnumj]][1, ], type = "l", ylim = c(min(Xnew[[plotnumj]]) * 0.8, max(Xnew[[plotnumj]]) * 
                  1.2), xlab = "", ylab = "")
                for (k in 1:Mrep) {
                  lines(c(0, timessimu), Xnew[[plotnumj]][k, ])
                }
                lines(times, Xtrue[plotnumj, ], col = "red", lwd = 2)
                
                plot(1:N/N, sort(q[plotnumj, ]), xlab = "", ylab = "")
                abline(0, 1)
            }
        }
        
    }
    
    
    
    if (!missing(numj)) {
        
        
        if (dim(x@gridf)[1] == 2) {
            
            phihat <- x@estimphi[, numj]
            
            if (x@model == "OU") {
                Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "EA", theta = c(phihat, 
                  sig), model = "OU", M = Mrep))
            }
            if (x@model == "CIR") {
                Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(phihat, 
                  sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
            }
        }
        
        if (dim(x@gridf)[1] == 1) {
            phihat <- x@estimphi[numj]
            
            if (x@estim.fix == 1) {
                paramfixed <- x@estim.fixed
            }
            if (x@estim.fix == 0) {
                paramfixed <- x@fixed
            }
            if (sum(x@random) == 1) {
                if (x@model == "OU") {
                  Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "EA", theta = c(phihat, 
                    paramfixed, sig), model = "OU", M = Mrep))
                }
                
                if (x@model == "CIR") {
                  Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(phihat, 
                    paramfixed, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
                }
            }
            if (sum(x@random) == 2) {
                if (x@model == "OU") {
                  Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "EA", theta = c(paramfixed, 
                    phihat, sig), model = "OU", M = Mrep))
                }
                if (x@model == "CIR") {
                  Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, method = "milstein", theta = c(paramfixed, 
                    phihat, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x))))
                }
            }
        }
        vecttimes <- intersect(round(timessimu, 10), round(times, 10))
        N <- length(vecttimes)
        
        q <- rep(0, N)
        for (i in 1:N) {
            q[i] <- sum(Xtrue[numj, which(times == vecttimes[i])[1]] > Xnew[, which(timessimu == vecttimes[i])[1]])/Mrep
        }
        
        if (plot.valid == 1) {
            plotnumj <- numj
            op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            plot(c(0, timessimu), Xnew[1, ], type = "l", ylim = c(min(Xnew) * 0.8, max(Xnew) * 1.2), xlab = "", ylab = "")
            for (k in 1:Mrep) {
                lines(c(0, timessimu), Xnew[k, ])
            }
            lines(times, Xtrue[plotnumj, ], col = "red", lwd = 2)
            
            
            plot(1:N/N, sort(q), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
            abline(0, 1)
        }
    }
    return(list(quantiles = q, Xnew = Xnew, plotnumj = plotnumj))
})


#' Validation of the chosen model.
#' 
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Bayes.fit class
#' @param Mrep number of trajectories to be drawn
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.valid logical(1), if TRUE, the results are depicted grafically
#' @param numj optional number of series to be validated
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setMethod(f = "valid", signature = "Bayes.fit", definition = function(x, Mrep = 100, newwindow = FALSE, plot.valid = TRUE, numj, ...) {
    if (newwindow) {
        x11(width = 10)
    }
    original.settings <- par(no.readonly = TRUE)
    
    times <- round(x@times, 10)
    N <- length(times) - 1
    Xtrue <- x@X[-x@ind.4.prior, ]
    Tend <- max(times)
    del <- round(min(diff(times)), 10)
    timessimu <- round(seq(del, Tend, by = del), 10)
    
    vecttimes <- intersect(round(timessimu, 10), round(times, 10))
    
    M <- dim(Xtrue)[1]
    sigma <- sqrt(x@sigma2)
    
    if (missing(numj)) {
        
        Xnew <- as.list(1:M)  # for each phihat_j a new sample size Mrep
        
        if (length(x@random) == 2) {
            
            phihat <- rbind(apply(x@alpha, 2, mean), apply(x@beta, 2, mean))
            
            for (j in 1:M) {
                Xnew[[j]] <- mixedsde.sim(M = Mrep, Tend, N = N, model = x@model, random = x@random, density.phi = "normalnormal", param = c(phihat[1, 
                  j], 0, phihat[2, j], 0), sigma = sigma, X0 = Xtrue[j, 1], op.plot = 0)$X
            }
        } else {
            if (x@random == 1) {
                alphahat <- apply(x@alpha, 2, mean)
                betahat <- mean(x@beta)
                
                for (j in 1:M) {
                  Xnew[[j]] <- mixedsde.sim(M = Mrep, Tend, N = length(timessimu), model = x@model, random = x@random, fixed = betahat, 
                    density.phi = "normal", param = c(alphahat[j], 0), sigma = sigma, X0 = Xtrue[j, 1], op.plot = 0)$X
                }
            }
            if (x@random == 2) {
                betahat <- apply(x@beta, 2, mean)
                alphahat <- mean(x@alpha)
                
                for (j in 1:M) {
                  Xnew[[j]] <- mixedsde.sim(M = Mrep, Tend, N = length(timessimu), model = x@model, random = x@random, fixed = alphahat, 
                    density.phi = "normal", param = c(betahat[j], 0), sigma = sigma, X0 = Xtrue[j, 1], op.plot = 0)$X
                }
            }
        }
        
        N <- length(vecttimes)
        
        q <- matrix(0, M, N)
        for (j in 1:M) {
            for (i in 1:N) {
                q[j, i] <- sum(Xtrue[j, i] > Xnew[[j]][, which(timessimu == vecttimes[i])[1]])/Mrep
            }
        }
        
        if (plot.valid == 1) {
            op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            plotnumj <- sample(1:M, 1)
            
            plot(c(0, timessimu), Xnew[[plotnumj]][1, ], type = "l", ylim = range(c(range(Xnew[[plotnumj]]), Xtrue[plotnumj, ])), xlab = "t", 
                ylab = expression(X[t]))
            for (k in 1:Mrep) {
                lines(c(0, timessimu), Xnew[[plotnumj]][k, ])
            }
            lines(times, Xtrue[plotnumj, ], col = "red", lwd = 2)
            
            plot(1:N/N, sort(q[plotnumj, ]), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
            abline(0, 1)
        }
        
    }
    
    
    
    if (!missing(numj)) {
        
        
        if (length(x@random) == 2) {
            
            phihat <- c(mean(x@alpha[, numj]), mean(x@beta[, numj]))
            
            
            Xnew <- mixedsde.sim(M = Mrep, Tend, N = length(timessimu), model = x@model, random = x@random, density.phi = "normalnormal", 
                param = c(phihat[1], 0, phihat[2], 0), sigma = sigma, X0 = Xtrue[numj, 1], op.plot = 0)$X
            
        } else {
            if (x@random == 1) {
                alphahat <- mean(x@alpha[, numj])
                betahat <- mean(x@beta)
                
                Xnew <- mixedsde.sim(M = Mrep, Tend, N = length(timessimu), model = x@model, random = x@random, fixed = betahat, density.phi = "normal", 
                  param = c(alphahat, 0), sigma = sigma, X0 = Xtrue[numj, 1], op.plot = 0)$X
            }
            if (x@random == 2) {
                betahat <- mean(x@beta[, numj])
                alphahat <- mean(x@alpha)
                
                Xnew <- mixedsde.sim(M = Mrep, Tend, N = N, model = x@model, random = x@random, fixed = alphahat, density.phi = "normal", 
                  param = c(betahat, 0), sigma = sigma, X0 = Xtrue[numj, 1], op.plot = 0)$X
            }
        }
        
        N <- length(vecttimes)
        q <- rep(0, N)
        for (i in 1:N) {
            q[i] <- sum(Xtrue[numj, which(times == vecttimes[i])[1]] > Xnew[, which(timessimu == vecttimes[i])[1]])/Mrep
        }
        
        
        
        if (plot.valid == 1) {
            plotnumj <- numj
            op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            plot(c(0, timessimu), Xnew[1, ], type = "l", ylim = range(c(range(Xnew), Xtrue[plotnumj, ])), xlab = "t", ylab = expression(X[t]))
            for (k in 1:Mrep) {
                lines(c(0, timessimu), Xnew[k, ])
            }
            lines(times, Xtrue[plotnumj, ], col = "red", lwd = 2)
            
            
            plot(1:N/N, sort(q), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
            abline(0, 1)
        }
    }
    # set plot settings back
    par(original.settings)
    
    invisible(list(quantiles = q, Xnew = Xnew, plotnumj = plotnumj))
})

########################################################### PREDICTION

#' Prediction method
#' 
#' @description Prediction
#' @param x Freq.fit or Bayes.fit class
#' @param ... other optional parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 


setGeneric("pred", function(x, ...) {
    standardGeneric("pred")
})


#' Prediction method for the Freq.fit class object
#' 
#' @description Frequentist prediction
#' @param x Freq.fit class
#' @param invariant 1 if the initial value is from the invariant distribution, default X0 is fixed from Xtrue
#' @param level alpha for the predicion intervals, default 0.05
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.pred logical(1), if TRUE, the results are depicted grafically
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.
#' 
setMethod(f = "pred", signature = "Freq.fit", definition = function(x, invariant = 0, level = 0.05, newwindow = FALSE, plot.pred = TRUE, 
    ...) {
    if (newwindow) {
        x11(width = 10)
    }
    Xtrue <- x@X
    timestrue <- x@times
    T <- timestrue[length(timestrue)]
    sig <- sqrt(x@sigma2)
    
    
    if (dim(x@gridf)[1] == 1) {
        index <- x@index
        M <- length(index)
        N <- dim(Xtrue)[2] - 1
        delta <- T/N
        Xpred <- matrix(0, M, N + 1)
        times <- seq(0, T, by = delta)
        
        if (x@estim.fix == 1) {
            paramfixed <- x@estim.fixed
        }
        if (x@estim.fix == 0) {
            paramfixed <- x@fixed
        }
        
        
        if (sum(x@random) == 1) {
            phipred <- rep(0, M)
            if (x@bic == 0) {
                p <- x@estimf/sum(x@estimf)
                for (i in 1:M) {
                  phipred[i] <- discr(x@gridf, p)
                }
            }
            if (x@bic != 0) {
                phipred <- rnorm(M, x@mu, x@omega)
            }
            
            if (x@model == "OU") {
                indexpred <- 1:M
                for (j in 1:M) {
                  if (invariant == FALSE) {
                    Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N, delta = T/N, method = "EA", theta = c(phipred[j], paramfixed, 
                      sig), model = "OU")
                  }
                  if (invariant == TRUE) {
                    X0 <- phipred[j]/paramfixed + (sig/(sqrt(2 * paramfixed))) * rnorm(1)
                    Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "EA", theta = c(phipred[j], paramfixed, sig), model = "OU")
                  }
                }
            }
            if (x@model == "CIR") {
                indexpred <- which(phipred > 0)
                phipred <- phipred[indexpred]
                Mpred <- length(phipred)
                Xpred <- matrix(0, Mpred, N + 1)
                
                for (j in 1:Mpred) {
                  if (invariant == 0) {
                    Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j], 1], N = N, delta = T/N, method = "milstein", theta = c(phipred[j], 
                      paramfixed, sig), model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                  }
                  if (invariant == 1) {
                    X0 <- rgamma(1, 2 * phipred[j]/sig^2, scale = sig^2/(2 * paramfixed))
                    Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "milstein", theta = c(phipred[j], paramfixed, sig), 
                      model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                  }
                }
            }
        }
        
        if (sum(x@random) == 2) {
            
            phipred <- rep(0, M)
            if (x@bic == 0) {
                p <- x@estimf/sum(x@estimf)
                for (i in 1:M) {
                  phipred[i] <- discr(x@gridf, p)
                }
            }
            if (x@bic != 0) {
                phipred <- rnorm(M, x@mu, x@omega)
            }
            
            if (x@model == "OU") {
                indexpred <- which(phipred > 0)
                phipred <- phipred[indexpred]
                Mpred <- length(indexpred)
                Xpred <- matrix(0, Mpred, N + 1)
                
                for (j in 1:Mpred) {
                  if (invariant == 0) {
                    Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N, delta = T/N, method = "EA", theta = c(paramfixed, phipred[j], 
                      sig), model = "OU")
                  }
                  if (invariant == 1) {
                    X0 <- (sig/(sqrt(2 * phipred[j]))) * rnorm(1)
                    Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "EA", theta = c(paramfixed, phipred[j], sig), model = "OU")
                  }
                }
            }
            if (x@model == "CIR") {
                indexpred <- which(phipred > 0)
                phipred <- phipred[indexpred]
                Mpred <- length(phipred)
                Xpred <- matrix(0, Mpred, N + 1)
                
                for (j in 1:Mpred) {
                  if (invariant == 1) {
                    if (paramfixed == 0) {
                      print("no invariant distribution, please fix X0")
                    }
                    if (paramfixed != 0) {
                      X0 <- rgamma(1, 2 * paramfixed/sig^2, scale = sig^2/(2 * phipred[j]))
                      Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "milstein", theta = c(paramfixed, phipred[j], 
                        sig), model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                    }
                  }
                  if (invariant == 0) {
                    Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j], 1], N = N, delta = T/N, method = "milstein", theta = c(paramfixed, 
                      phipred[j], sig), model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                  }
                }
            }
        }
        
        if (plot.pred == TRUE) {
            op <- par(mfrow = c(1, 3), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            plot(sort(x@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(x@estimphi, phipred) * 0.8, max(x@estimphi, phipred) * 
                1.2), ylim = c(min(x@estimphi, phipred) * 0.8, max(x@estimphi, phipred) * 1.2), ylab = "", xlab = "")
            abline(0, 1)
            plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 
                1.5), main = "True trajectories")
            for (j in indexpred) {
                lines(timestrue, Xtrue[j, ], col = j)
            }
            
            plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
            for (j in 1:length(indexpred)) {
                lines(times, Xpred[j, ], col = 1)
            }
            l.bound = level/2
            u.bound = 1 - level/2
            PI <- matrix(0, 2, N + 1)
            for (k in 1:(N + 1)) {
                PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
            }
            lines(times, PI[1, ], col = "green", lwd = 2)
            lines(times, PI[2, ], col = "green", lwd = 2)
        }
        return(list(phipred = phipred, Xpred = Xpred, indexpred = indexpred))
    }
    
    if (dim(x@gridf)[1] == 2) {
        
        
        index <- x@index
        M <- length(index)
        N <- dim(Xtrue)[2] - 1
        delta <- T/N
        Xpred <- matrix(0, M, N + 1)
        times <- seq(0, T, by = delta)
        
        phipred <- matrix(0, 2, M)
        
        if (x@bic != 0) {
            mu <- x@mu
            omega <- x@omega
            mu1 <- mu[1]
            mu2 <- mu[2]
            omega1 <- omega[1]
            omega2 <- omega[2]
            phipred[1, ] <- rnorm(M, mu1, omega1)
            phipred[2, ] <- rnorm(M, mu2, omega2)
        }
        if (x@bic == 0) {
            gridf1 <- x@gridf[1, ]
            gridf2 <- x@gridf[2, ]
            
            
            marg1 <- ((max(gridf2) - min(gridf2))/length(gridf2)) * apply(x@estimf, 1, sum)
            marg2 <- ((max(gridf1) - min(gridf1))/length(gridf1)) * apply(x@estimf, 2, sum)
            p1 <- marg1/sum(marg1)
            p2 <- marg2/sum(marg2)
            phipred <- matrix(0, 2, sum(M))
            for (i in 1:M) {
                phipred[1, i] <- discr(gridf1, p1)
                phipred[2, i] <- discr(gridf2, p2)
            }
            
        }
        
        if (x@model == "OU") {
            indexpred <- which(phipred[2, ] > 0)
            phipred <- phipred[, indexpred]
            Mpred <- length(indexpred)
            Xpred <- matrix(0, Mpred, N + 1)
            for (j in 1:Mpred) {
                if (invariant == 0) {
                  Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N, delta = T/N, method = "EA", theta = c(phipred[, j], sig), model = "OU")
                }
                if (invariant == 1) {
                  X0 <- phipred[1, j]/phipred[2, j] + (sig/(sqrt(2 * phipred[2, j]))) * rnorm(1)
                  Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "EA", theta = c(phipred[, j], sig), model = "OU")
                }
            }
        }
        if (x@model == "CIR") {
            indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
            phipred <- phipred[, indexpred]
            Mpred <- length(indexpred)
            Xpred <- matrix(0, Mpred, N + 1)
            
            for (j in 1:Mpred) {
                if (invariant == 0) {
                  Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j], 1], N = N, delta = T/N, method = "milstein", theta = c(phipred[, 
                    j], sig), model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                }
                if (invariant == 1) {
                  X0 <- rgamma(1, 2 * phipred[1, j]/sig^2, scale = sig^2/(2 * phipred[2, j]))
                  Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N, method = "milstein", theta = c(phipred[, j], sig), model = "CIR", 
                    sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * sqrt(x)))
                }
            }
        }
        
        if (plot.pred == TRUE) {
            
            op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            
            
            plot(sort(x@estimphi[1, indexpred]), sort(phipred[1, ]), pch = 18, xlim = c(min(x@estimphi[1, ], phipred[1, ]) * 0.8, max(x@estimphi[1, 
                ], phipred[1, ]) * 1.2), ylim = c(min(x@estimphi[1, ], phipred[1, ]) * 0.8, max(x@estimphi[1, ], phipred[1, ]) * 1.2), 
                ylab = "", xlab = "", main = "First random effect")
            abline(0, 1)
            plot(sort(x@estimphi[2, indexpred]), sort(phipred[2, ]), pch = 18, xlim = c(min(x@estimphi[2, ], phipred[2, ]) * 0.8, max(x@estimphi[2, 
                ], phipred[2, ]) * 1.2), ylim = c(min(x@estimphi[2, ], phipred[2, ]) * 0.8, max(x@estimphi[2, ], phipred[2, ]) * 1.2), 
                ylab = "", xlab = "", main = "Second random effect")
            abline(0, 1)
            
            
            plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 
                1.5), main = "True trajectories")
            for (j in indexpred) {
                lines(times, Xtrue[j, ], col = j)
            }
            
            plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), 
                main = "Predictive trajectories")
            for (j in 1:length(indexpred)) {
                lines(timestrue, Xpred[j, ], col = 1)
            }
            
            l.bound = level/2
            u.bound = 1 - level/2
            PI <- matrix(0, 2, N + 1)
            for (k in 1:(N + 1)) {
                PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
            }
            lines(timestrue, PI[1, ], col = "green", lwd = 2)
            lines(timestrue, PI[2, ], col = "green", lwd = 2)
            
        }
        return(list(phipred = phipred, Xpred = Xpred, indexpred = indexpred))
    }
})

# setGeneric('Bayes.pred', function(estim, level = 0.05, newwindow = FALSE, plot.pred = TRUE, plot.legend = TRUE, burnIn, thinning,
# ...) { standardGeneric('Bayes.pred') })

######## 
#' Bayesian prediction method for a class object Bayes.fit
#' 
#' @description Bayesian prediction
#' @param x Bayes.fit class
#' @param invariant logical(1), if TRUE, the initial value is from the invariant distribution \eqn{X_t~N(\alpha/\beta, \sigma^2/2\beta)} for the OU and \eqn{X_t~\Gamma(2\alpha/\sigma^2, \sigma^2/2\beta)} for the CIR process, if FALSE (default) X0 is fixed from the data starting points
#' @param level alpha for the predicion intervals, default 0.05
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.pred logical(1), if TRUE, the results are depicted grafically
#' @param plot.legend logical(1), if TRUE, a legend is added to the plot
#' @param burnIn optional, if missing, the proposed value of the mixedsde.fit function is taken
#' @param thinning optional, if missing, the proposed value of the mixedsde.fit function is taken
#' @param only.interval logical(1), if TRUE, only prediction intervals are calculated, much faster than sampling from the whole predictive distribution
#' @param sample.length number of samples to be drawn from the predictive distribution, if only.interval = FALSE
#' @param cand.length number of candidates for which the predictive density is calculated, i.e. the candidates to be drawn from
#' @param trajectories logical(1), if TRUE, only trajectories are drawn from the point estimations instead of sampling from the predictive distribution, similar to the frequentist approach
#' @param ylim optional
#' @param xlab optional, default 'times'
#' @param ylab optional, default 'X'
#' @param col color for the prediction intervals, default 3
#' @param lwd linewidth for the prediction intervals, default 3 
#' @param ... optional plot parameters
#' @references 
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: an R package to fit mixed stochastic differential equations.


setMethod(f = "pred", signature = "Bayes.fit", definition = function(x, invariant = FALSE, level = 0.05, newwindow = FALSE, plot.pred = TRUE, 
    plot.legend = TRUE, burnIn, thinning, only.interval = TRUE, sample.length = 500, cand.length = 100, trajectories = FALSE, ylim, xlab = "times", 
    ylab = "X", col = 3, lwd = 2, ...) {
    if (newwindow) {
        x11(width = 10)
    }
    original.settings <- par(no.readonly = TRUE)
    
    if (missing(burnIn)) 
        burnIn <- x@burnIn
    if (missing(thinning)) 
        thinning <- x@thinning
    
    random <- x@random
    
    est <- chain2samples(x, burnIn, thinning)
    K <- length(est@sigma2)
    densPhi <- function(t, samples) mean(dnorm(t, samples$mu, sqrt(samples$omega)))
    M <- nrow(x@X)
    
    ##### 
    model <- x@model
    
    if (trajectories) {
        sigma2 <- mean(est@sigma2)
        
        
        if (length(random) == 2) {
            
            cand1 <- seq(mean(est@mu[, 1]) - 4 * sqrt(mean(est@omega[, 1])), mean(est@mu[, 1]) + 4 * sqrt(mean(est@omega[, 1])), length = 500)
            cand2 <- seq(mean(est@mu[, 2]) - 4 * sqrt(mean(est@omega[, 2])), mean(est@mu[, 2]) + 4 * sqrt(mean(est@omega[, 2])), length = 500)
            
            phi.pred <- matrix(0, M, 2)
            
            prob <- sapply(cand1, densPhi, samples = list(mu = est@mu[, 1], omega = est@omega[, 1]))
            phi.pred[, 1] <- replicate(M, discr(cand1, prob))
            prob <- sapply(cand2, densPhi, samples = list(mu = est@mu[, 2], omega = est@omega[, 2]))
            phi.pred[, 2] <- replicate(M, discr(cand2, prob))
            
            if (invariant) {
                if (model == "OU") {
                  X0 <- rnorm(M, phi.pred[, 1]/phi.pred[, 2], sqrt(sigma2/(2 * phi.pred[, 2])))
                }
                if (model == "CIR") {
                  X0 <- rgamma(M, 2 * phi.pred[, 1]/sigma2, scale = sigma2/(2 * phi.pred[, 2]))
                }
                
            } else {
                X0 <- x@X[, 1]
            }
            
        } else {
            
            cand <- seq(mean(est@mu) - 4 * sqrt(mean(est@omega)), mean(est@mu) + 4 * sqrt(mean(est@omega)), length = 500)
            prob <- sapply(cand, densPhi, samples = list(mu = est@mu, omega = est@omega))
            phi.pred <- replicate(M, discr(cand, prob))
            
            if (invariant) {
                if (model == "OU") {
                  if (random == 1) 
                    X0 <- rnorm(M, phi.pred/mean(est@beta), sqrt(sigma2/(2 * mean(est@beta))))
                  if (random == 2) 
                    X0 <- rnorm(M, mean(est@alpha)/phi.pred, sqrt(sigma2/(2 * phi.pred)))
                }
                if (model == "CIR") {
                  if (random == 1) 
                    X0 <- rgamma(M, 2 * phi.pred/sigma2, scale = sigma2/(2 * mean(est@beta)))
                  if (random == 2) 
                    X0 <- rgamma(M, 2 * mean(est@alpha)/sigma2, scale = sigma2/(2 * phi.pred))
                }
            } else {
                X0 <- x@X[, 1]
            }
            
        }
        
        
        Xpred <- matrix(0, M, length(x@times))
        Xpred[, 1] <- X0
        dt <- diff(x@times)
        if (model == "OU") {
            m.traj <- function(t, phi, Xn_1) phi[1]/phi[2] + (Xn_1 - phi[1]/phi[2]) * exp(-phi[2] * t)
            v.traj <- function(t, phi) phi[3] * (1 - exp(-2 * phi[2] * t))/2/phi[2]
            
            if (length(random) == 2) {
                for (i in 2:length(x@times)) {
                  for (j in 1:M) {
                    Xpred[j, i] <- rnorm(1, m.traj(dt[i - 1], phi.pred[j, ], Xpred[j, i - 1]), sqrt(v.traj(dt[i - 1], c(phi.pred[j, ], 
                      sigma2))))
                  }
                }
            } else {
                if (random == 1) {
                  beta <- mean(est@beta)
                  for (i in 2:length(x@times)) {
                    for (j in 1:M) {
                      Xpred[j, i] <- rnorm(1, m.traj(dt[i - 1], c(phi.pred[j], beta), Xpred[j, i - 1]), sqrt(v.traj(dt[i - 1], c(phi.pred[j], 
                        beta, sigma2))))
                    }
                  }
                } else {
                  alpha <- mean(est@alpha)
                  for (i in 2:length(x@times)) {
                    for (j in 1:M) {
                      Xpred[j, i] <- rnorm(1, m.traj(dt[i - 1], c(alpha, phi.pred[j]), Xpred[j, i - 1]), sqrt(v.traj(dt[i - 1], c(alpha, 
                        phi.pred[j], sigma2))))
                    }
                  }
                  
                }
            }
            
        } else {
            
            likeli.CIR <- function(x, t, phi, sigma2, Xn_1) dcCIR2(x, t, Xn_1, c(phi, sqrt(sigma2)))
            cand <- function(i) {
                he <- x@X[, i]
                seq(min(he) - abs(min(he)) * 0.5, max(he) + abs(max(he)) * 0.5, length = cand.length)
            }
            
            if (length(random) == 2) {
                for (i in 2:length(x@times)) {
                  for (j in 1:M) {
                    prob <- likeli.CIR(cand(i), dt[i - 1], phi.pred[j, ], sigma2, Xpred[j, i - 1])
                    Xpred[j, i] <- discr(cand(i), prob)
                  }
                }
                
            } else {
                if (random == 1) {
                  beta <- mean(est@beta)
                  for (i in 2:length(x@times)) {
                    for (j in 1:M) {
                      prob <- likeli.CIR(cand(i), dt[i - 1], c(phi.pred[j], beta), sigma2, Xpred[j, i - 1])
                      Xpred[j, i] <- discr(cand(i), prob)
                    }
                  }
                  
                } else {
                  # random == 2
                  alpha <- mean(est@alpha)
                  for (i in 2:length(x@times)) {
                    for (j in 1:M) {
                      prob <- likeli.CIR(cand(i), dt[i - 1], c(alpha, phi.pred[j]), sigma2, Xpred[j, i - 1])
                      Xpred[j, i] <- discr(cand(i), prob)
                    }
                  }
                }
            }
        }
        
        
        if (plot.pred == TRUE) {
            
            op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            if (missing(ylim)) 
                ylim <- range(c(range(x@X), range(Xpred)))
            
            plot(x@times, x@X[1, ], type = "l", ylim = ylim, xlab = xlab, ylab = ylab, ...)
            for (i in 2:M) lines(x@times, x@X[i, ], ...)
            
            for (i in 1:M) lines(x@times, Xpred[i, ], col = col, lwd = lwd, ...)
            for (i in 2:M) lines(x@times, x@X[i, ], ...)
            
            if (plot.legend) 
                legend("bottomright", c("data", "drawn trajectories"), lty = 1, col = c(1, col), lwd = c(1, lwd), cex = 0.7, box.lty = 0, 
                  inset = 0.01)
        }
        
        return(new(Class = "Bayes.pred", phi.pred = as.matrix(phi.pred), Xpred = Xpred, estim = out(x)))
        
    } else {
        # end if(trajectories)
        
        if (length(random) == 2) {
            
            cand1 <- seq(mean(est@mu[, 1]) - 4 * sqrt(mean(est@omega[, 1])), mean(est@mu[, 1]) + 4 * sqrt(mean(est@omega[, 1])), length = 500)
            cand2 <- seq(mean(est@mu[, 2]) - 4 * sqrt(mean(est@omega[, 2])), mean(est@mu[, 2]) + 4 * sqrt(mean(est@omega[, 2])), length = 500)
            
            phi.pred <- matrix(0, K, 2)
            
            prob <- sapply(cand1, densPhi, samples = list(mu = est@mu[, 1], omega = est@omega[, 1]))
            phi.pred[, 1] <- replicate(K, discr(cand1, prob))
            prob <- sapply(cand2, densPhi, samples = list(mu = est@mu[, 2], omega = est@omega[, 2]))
            phi.pred[, 2] <- replicate(K, discr(cand2, prob))
            
            if (invariant) {
                if (model == "OU") {
                  densX0 <- function(t) mean(dnorm(t, phi.pred[, 1]/phi.pred[, 2], sqrt(est@sigma2/(2 * phi.pred[, 2]))))
                  cand <- seq(min(x@X[, 1]) - abs(max(x@X[, 1])) * 0.5, max(x@X[, 1]) + abs(max(x@X[, 1])) * 0.9, length = 1000)
                }
                if (model == "CIR") {
                  densX0 <- function(t) mean(dgamma(t, 2 * phi.pred[, 1]/est@sigma2, scale = est@sigma2/(2 * phi.pred[, 2])))
                  cand <- seq(max(0.001, min(x@X[, 1]) - abs(max(x@X[, 1])) * 0.5), max(x@X[, 1]) + abs(max(x@X[, 1])) * 0.9, length = 1000)
                }
                prob <- sapply(cand, densX0)
                X0 <- replicate(K, discr(cand, prob))
            } else {
                if (length(unique(x@X[, 1])) == 1) 
                  X0 <- x@X[1, 1]
                if (length(unique(x@X[, 1])) > 1) 
                  X0 <- rep(x@X[, 1], length = K)
            }
            
        } else {
            
            cand <- seq(mean(est@mu) - 4 * sqrt(mean(est@omega)), mean(est@mu) + 4 * sqrt(mean(est@omega)), length = 500)
            prob <- sapply(cand, densPhi, samples = list(mu = est@mu, omega = est@omega))
            phi.pred <- replicate(K, discr(cand, prob))
            
            if (invariant) {
                if (model == "OU") {
                  if (random == 1) 
                    densX0 <- function(t) mean(dnorm(t, phi.pred/est@beta, sqrt(est@sigma2/(2 * est@beta))))
                  if (random == 2) 
                    densX0 <- function(t) mean(dnorm(t, est@alpha/phi.pred, sqrt(est@sigma2/(2 * phi.pred))))
                  cand <- seq(min(x@X[, 1]) - abs(max(x@X[, 1])) * 0.5, max(x@X[, 1]) + abs(max(x@X[, 1])) * 0.9, length = 1000)
                }
                if (model == "CIR") {
                  if (random == 1) 
                    densX0 <- function(t) mean(dgamma(t, 2 * phi.pred/est@sigma2, scale = est@sigma2/(2 * est@beta)))
                  if (random == 2) 
                    densX0 <- function(t) mean(dgamma(t, 2 * est@alpha/est@sigma2, scale = est@sigma2/(2 * phi.pred)))
                  cand <- seq(max(0.001, min(x@X[, 1]) - abs(max(x@X[, 1])) * 0.5), max(x@X[, 1]) + abs(max(x@X[, 1])) * 0.9, length = 1000)
                }
                prob <- sapply(cand, densX0)
                X0 <- replicate(K, discr(cand, prob))
            } else {
                if (length(unique(x@X[, 1])) == 1) 
                  X0 <- x@X[1, 1]
                if (length(unique(x@X[, 1])) > 1) 
                  X0 <- rep(x@X[, 1], length = K)
                
            }
            
        }
        
        
        cand <- function(i) {
            he <- x@X[, i + 1]
            seq(min(he) - abs(min(he)) * 0.5, max(he) + abs(max(he)) * 0.5, length = cand.length)
        }
        lt <- length(x@times) - 1
        l.bound = level/2
        u.bound = 1 - level/2
        
        dcCIR2.vec <- function(x, t, x0, theta, log = FALSE) {
            c <- 2 * theta[, 2]/((1 - exp(-theta[, 2] * t)) * theta[, 3]^2)
            ncp <- 2 * c * x0 * exp(-theta[, 2] * t)
            df <- 4 * theta[, 1]/theta[, 3]^2
            lik <- (dchisq(2 * x * c, df = df, ncp = ncp, log = TRUE) + log(2 * c))
            if (!log) 
                lik <- exp(lik)
            lik
        }
        
        if (length(random) == 2) {
            
            if (model == "OU") {
                m <- function(t) phi.pred[1]/phi.pred[, 2] + (X0 - phi.pred[, 1]/phi.pred[, 2]) * exp(-phi.pred[, 2] * t)
                v <- function(t) est@sigma2 * (1 - exp(-2 * phi.pred[, 2] * t))/2/phi.pred[, 2]
                likeli <- function(x, t) dnorm(x, m(t), sqrt(v(t)))
            } else {
                likeli <- function(x, t) dcCIR2.vec(x, t, X0, cbind(phi.pred, sqrt(est@sigma2)))
            }
            
            
        } else {
            if (random == 1) {
                
                if (model == "OU") {
                  m <- function(t) phi.pred/est@beta + (X0 - phi.pred/est@beta) * exp(-est@beta * t)
                  v <- function(t) est@sigma2 * (1 - exp(-2 * est@beta * t))/2/est@beta
                  likeli <- function(x, t) dnorm(x, m(t), sqrt(v(t)))
                } else {
                  likeli <- function(x, t) dcCIR2.vec(x, t, X0, cbind(phi.pred, est@beta, sqrt(est@sigma2)))
                }
                
            } else {
                # random == 2
                
                if (model == "OU") {
                  m <- function(t) est@alpha/phi.pred + (X0 - est@alpha/phi.pred) * exp(-phi.pred * t)
                  v <- function(t) est@sigma2 * (1 - exp(-2 * phi.pred * t))/2/phi.pred
                  likeli <- function(x, t) dnorm(x, m(t), sqrt(v(t)))
                } else {
                  likeli <- function(x, t) dcCIR2.vec(x, t, X0, cbind(est@alpha, phi.pred, sqrt(est@sigma2)))
                }
                
            }
            
        }
        
        pred <- function(i) {
            densXi <- function(c) mean(likeli(c, x@times[i + 1] - x@times[1]))
            ca <- cand(i)
            prob <- sapply(ca, densXi)
            
            if (!only.interval) {
                samp.X <- replicate(sample.length, discr(ca, prob))
            } else {
                samp.X <- NULL
            }
            
            VF <- cumsum(prob)/sum(prob)
            qu.l <- ca[which(VF >= l.bound)[1]]
            qu.u <- ca[which(VF >= u.bound)[1]]
            return(list(samp.X = samp.X, qu.l = qu.l, qu.u = qu.u))
        }
        
        qu.l <- numeric(lt)
        qu.u <- numeric(lt)
        if (!only.interval) 
            Xpred <- matrix(0, sample.length, lt)
        for (i in 1:lt) {
            he <- pred(i)
            if (!only.interval) 
                Xpred[, i] <- he$samp.X
            qu.u[i] <- he$qu.u
            qu.l[i] <- he$qu.l
        }
        
        cr <- apply(qu.l <= t(x@X[, -1]) & qu.u >= t(x@X[, -1]), 1, mean)
        
        if (plot.pred == TRUE) {
            
            op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, 
                cex.axis = 0.7)
            if (missing(ylim)) 
                ylim <- range(c(min(qu.l), max(qu.u), range(x@X)))
            
            plot(x@times[-1], qu.l, type = "l", ylim = ylim, xlab = xlab, ylab = ylab, col = col, lwd = lwd, ...)
            lines(x@times[-1], qu.u, col = col, lwd = lwd, ...)
            
            for (i in 1:nrow(x@X)) lines(x@times, x@X[i, ])
            lines(x@times[-1], qu.l, col = col, lwd = lwd, ...)
            lines(x@times[-1], qu.u, col = col, lwd = lwd, ...)
            
            if (plot.legend) 
                legend("bottomright", c("data", "prediction intervals"), lty = 1, col = c(1, col), lwd = c(1, lwd), cex = 0.7, box.lty = 0, 
                  inset = 0.01)
            
            plot(x@times[-1], cr, ylim = c(min(cr) * 0.9, max(c(cr), 1)), type = "l", xlab = xlab, ylab = "coverage rates")
            abline(h = 1 - level, col = 2, lty = 2)
            if (plot.legend) 
                legend("bottomright", paste((1 - level) * 100, "%", sep = ""), lty = 2, col = 2, cex = 0.7, box.lty = 0, inset = 0.01)
            
        }
        # set plot settings back
        par(original.settings)
        
        if (!only.interval) {
            return(new(Class = "Bayes.pred", phi.pred = as.matrix(phi.pred), Xpred = Xpred, coverage.rate = cr, qu.l = qu.l, qu.u = qu.u, 
                estim = out(x)))
        } else {
            return(new(Class = "Bayes.pred", phi.pred = as.matrix(phi.pred), coverage.rate = cr, qu.l = qu.l, qu.u = qu.u, estim = out(x)))
        }
    }
}) 

# Functions for EM estimation of network tomography
# Based on methods of Cao et al. (JASA, 2000)


#' Simple initialization for phi in model of Cao et al. (2000)
#'
#' Uses a crude estimator to get a starting point for phi in the model of
#' Cao et al. (2000).
#'
#' @param Y matrix (n x k) of observed link loads over time
#' @param A routing matrix (m x k)
#' @param lambda0 numeric vector (length k) of initial guesses for lambda
#' @param c power parameter in model of Cao et al. (2000)
#' @return numeric starting value for phi
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
phi_init <- function(Y, A, lambda0, c) {
    exp(mean(log(diag(cov(Y))/ (A %*% lambda0^c) )))
}

#' Q function for locally IID EM algorithm of Cao et al. (2000)
#'
#' Computes the Q function (expected log-likelihood) for the EM algorithm of
#' Cao et al. (2000) for their locally IID model.
#'
#' @param logtheta numeric vector (length k+1) of log(lambda) (1:k) and log(phi)
#'      (last entry)
#' @param c power parameter in model of Cao et al. (2000)
#' @param M matrix (n x k) of conditional expectations for OD flows,
#'      one time per row
#' @param rdiag numeric vector (length k) containing diagonal of conditional
#'      covariance matrix R
#' @param epsilon numeric nugget to add to diagonal of covariance for numerical
#'      stability
#' @return numeric value of Q function; not vectorized in any way
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
Q_iid <- function(logtheta, c, M, rdiag, epsilon) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])
    phi <- exp(logtheta[length(logtheta)])

    # Calculate first part of Q :
    #  -1/2 sum_t (m_t - lambda)' \Sigma^{-1} (m_t - lambda)
    Q <- -1/2*sum(apply( M, 1, function(m) sum((m-lambda)^2 /
                                               (phi*lambda^c+epsilon) ) ))

    # Calculate remaining part:
    ## -T/2 (I*log(\phi) + \sum_i c*log(\lambda_i) + 
    # \sum_i r_ii/ \phi / lambda_i^c)
    # With nugget, revised to:
    # -T/2 ( \sum_i log(\phi * \lambda_i^c + \varepsilon) +
    # \sum_i r_ii/(\phi * \lambda_i^c + \varepsilon) )
    Q <- Q - nrow(M)/2*(sum(log( phi*lambda^c + epsilon )) +
                        sum(rdiag/(phi * lambda^c + epsilon) ) )

    return(Q);
}

#' Compute conditional expectations for EM algorithms of Cao et al. (2000)
#'
#' Computes conditional expectation of OD flows for E-step of EM algorithm from
#' Cao et al. (2000) for their locally IID model.
#'
#' @param yt numeric vector (length m) of link loads from single time
#' @param lambda numeric vector (length k) of mean OD flows from last M-step
#' @param phi numeric scalar scale for covariance matrix of xt
#' @param A routing matrix (m x k) for network being analyzed
#' @param c power parameter in model of Cao et al. (2000)
#' @param epsilon numeric nugget to add to diagonal of covariance for numerical
#'      stability
#' @return numeric vector of same size as lambda with conditional expectations
#'      of x
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
m_estep <- function(yt, lambda, phi, A, c, epsilon) {
    sigma <- diag_mat(phi*lambda^c + epsilon)
    
    m <- lambda
    m <- m  +  sigma %*% t(A) %*% solve(A %*% sigma %*% t(A), yt - A %*% lambda)
    
    return(m)
}

#' Compute conditional covariance matrix for EM algorithms of Cao et al. (2000)
#'
#' Computes conditional covariance of OD flows for E-step of EM algorithm from
#' Cao et al. (2000) for their locally IID model.
#'
#' @param lambda numeric vector (length k) of mean OD flows from last M-step
#' @param phi numeric scalar scale for covariance matrix of xt
#' @param A routing matrix (m x k) for network being analyzed
#' @param c power parameter in model of Cao et al. (2000)
#' @param epsilon numeric nugget to add to diagonal of covariance for numerical
#'      stability
#' @return conditional covariance matrix (k x k) of OD flows given parameters
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
R_estep <- function(lambda, phi, A, c, epsilon) {
    sigma <- diag_mat(phi*lambda^c + epsilon)
    AdotSigma <- A %*% sigma
    
    R <- diag(phi*lambda^c + epsilon)
    R <- R - t(AdotSigma) %*% solve( AdotSigma %*% t(A), AdotSigma )
    
    return(R)
}

#' Compute analytic gradient of Q-function for locally IID EM algorithm of Cao
#' et al. (2000)
#'
#' Computes gradient of Q-function with respect to log(c(lambda,phi)) for EM
#' algorithm from Cao et al. (2000) for their locally IID model.
#'
#' @param logtheta numeric vector (length k+1) of log(lambda) (1:k) and log(phi)
#'      (last entry)
#' @param c power parameter in model of Cao et al. (2000)
#' @param M matrix (n x k) of conditional expectations for OD flows,
#'      one time per row
#' @param rdiag numeric vector (length k) containing diagonal of conditional
#'      covariance matrix R
#' @param epsilon numeric nugget to add to diagonal of covariance for numerical
#'      stability
#' @return numeric vector of same length as logtheta containing calculated
#'      gradient
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
grad_iid <- function(logtheta, c, M, rdiag, epsilon) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])
    phi <- exp(logtheta[length(logtheta)])

    grad <- rep(NA, length(logtheta))

    # Calculate gradient for phi
    grad[length(grad)] <- (-nrow(M)/2*(ncol(M)/phi -
                                       1/phi/phi*sum(rdiag/lambda^c)) +
                           1/2/phi/phi* sum(apply(M, 1, function(m)
                                                  sum((m-lambda)^2/lambda^c))) )

    # Calculate gradient for each lambda
    grad[-length(grad)] <- (-nrow(M)/2*c*(1/lambda - rdiag/phi/lambda^(c+1)) +
                            (colSums(M)-nrow(M)*lambda)/phi/lambda^c +
                            c/2/phi*sapply(1:ncol(M), function(j)
                                           sum((M[,j]-lambda[j])^2 /
                                               lambda[j]^(c+1))))

    # Jacobian adjustment
    grad <- grad*exp(logtheta)

    # Return gradient
    return(grad)
}

#' Run EM algorithm to obtain MLE for locally IID model of Cao et al. (2000)
#'
#' Runs EM algorithm to compute MLE for the locally IID model of Cao et al.
#' (2000). Uses numerical optimization of Q-function for each M-step with
#' analytic computation of its gradient.
#'
#' @param Y matrix (h x k) of observations in local window; columns correspond
#'      to OD flows, and rows are individual observations
#' @param A routing matrix (m x k) for network being analyzed
#' @param lambda0 initial vector of values (length k) for lambda; \code{ipfp} is
#'      a good way to obtain this
#' @param phi0 initial value for covariance scale phi; initializes automatically
#'      using \code{phi_init} if NULL, but you can likely do better
#' @param c power parameter in model of Cao et al. (2000)
#' @param maxiter maximum number of EM iterations to run
#' @param tol tolerance (in relative change in Q function value) for stopping EM
#'      iterations
#' @param epsilon numeric nugget to add to diagonal of covariance for numerical
#'      stability
#' @param method optimization method to use (in optim calls)
#' @param checkActive logical check for deterministically known OD flows
#' @return list with 3 elements: \code{lambda}, the estimated value of lambda;
#'      \code{phi}, the estimated value of phi; and \code{iter}, the number of
#'      iterations run
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
locally_iid_EM <- function(Y, A, lambda0, phi0=NULL, c=2,
                           maxiter = 1e3, tol=1e-6, epsilon=0.01,
                           method="L-BFGS-B", checkActive=FALSE) {
    if (checkActive) {
        # Check for inactive (deterministically-known) OD flows
        activeLink <- which(apply(Y, 2, function(col) max(col) > 0))
        activeOD <- apply(Y, 1, getActive, A=A)
        activeOD <- which(apply(activeOD, 1, any))

        # Subset A and lambda
        A <- A[activeLink,activeOD]
        Y <- Y[,activeLink]
        lambda0 <- lambda0[activeOD]
        lambda <- lambda[activeOD]
    } else {
        activeOD <- seq(ncol(A))
    }
    
    # Determine number of latent variables
    p <- ncol(A)

    # Setup initial parameters
    lambda <- lambda0

    if (is.null(phi0)) {
        phi0 <- phi_init(Y, A, lambda0, c)
        phi <- phi0
    } else {
        phi <- phi0
    }

    # Calculate bounds
    lower <- rep(log(.Machine$double.eps)/c+log(max(Y)), p)
    upper <- rep(log(max(Y)), p)

    # Run first iteration
    # E step - calculation of m_t = E(x_t | y_t, \theta) and 
    #          R = Var(x_t | y_t, \theta)
    tryCatch(M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi, A=A, c=c,
                          epsilon=epsilon)),
             error = function(e) {
                 print(phi); print(lambda); print("M"); stop(e);
             },
             finally = NULL )
    R <- tryCatch(R_estep(lambda, phi=phi, A, c, epsilon),
                  error = function(e) {
                      print(phi); print(lambda); print("M"); stop(e);
                  },
                  finally = NULL )
    # M step - Optimize Q wrt theta
    mstep <- optim(log(c(lambda0,phi0)),
                   function(...) -Q_iid(...),# gr=function(...) -grad_iid(...),
                   c=c, M=M, rdiag=diag(R), epsilon=epsilon,
                   lower=lower, upper=upper,
                   method="L-BFGS-B")
    theta <- exp(mstep$par)
    lambda <- theta[-length(theta)]
    phi <- theta[length(theta)]

    # Calculate initial ll
    ll <- -mstep$value

    # Run EM iterations
    for (iter in 1:maxiter) {
        # E step - calculation of m_t = E(x_t | y_t, \theta) and 
        #  R = Var(x_t | y_t, \theta)
        tryCatch(M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi, A=A, c=c,
                              epsilon=epsilon)),
                 error = function(e) {
                     print(phi); print(lambda); print("M"); stop(e);
                 },
                 finally = NULL )
        R <- tryCatch( R_estep(lambda, phi=phi, A, c, epsilon),
                      error = function(e) {
                          print(phi); print(lambda); print("M"); stop(e);
                      },
                      finally = NULL )

        # M step - Optimize Q wrt theta
        mstep <- optim(par=log(c(lambda0,phi0)), fn=function(...) -Q_iid(...),
                       gr=function(...) -grad_iid(...), c=c, M=M, rdiag=diag(R),
                       epsilon=epsilon, lower=lower, upper=upper,
                       method="L-BFGS-B" )
        theta <- exp(mstep$par)
        lambda <- theta[-length(theta)]
        phi <- theta[length(theta)]
        ll.new <- -mstep$value

        # cat(ll, ll.new, abs(ll.new-ll)/abs(ll.new+ll)*2, "\n")

        # Check for convergence
        if (abs(ll.new-ll)/abs(ll.new+ll)*2 < tol) {
            ll <- ll.new
            break
        } else {
            ll <- ll.new
        }
    }
    tmp <- numeric(p)
    tmp[ activeOD ] <- lambda
    lambda <- tmp
    return(list(lambda=lambda, phi=phi, iter=iter))
}

#' Q function for smoothed EM algorithm of Cao et al. (2000)
#'
#' Computes the Q function (expected log-likelihood) for the EM algorithm of
#' Cao et al. (2000) for their smoothed model.
#'
#' @param logtheta numeric vector (length k+1) of log(lambda) (1:k) and log(phi)
#'      (last entry)
#' @param c power parameter in model of Cao et al. (2000)
#' @param M matrix (n x k) of conditional expectations for OD flows,
#'      one time per row
#' @param rdiag numeric vector (length k) containing diagonal of conditional
#'      covariance matrix R
#' @param eta0 numeric vector (length k+1) containing value for log(c(lambda,
#'      phi)) from previous time (or initial value)
#' @param sigma0 covariance matrix (k+1 x k+1) of log(c(lambda, phi)) from
#'      previous time (or initial value)
#' @param V evolution covariance matrix (k+1 x k+1) for log(c(lambda, phi))
#'      (random walk)
#' @param eps.lambda numeric small positive value to add to lambda for numerical
#'      stability; typically 0
#' @param eps.phi numeric small positive value to add to phi for numerical
#'      stability; typically 0
#' @return numeric value of Q function; not vectorized in any way
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
Q_smoothed <- function(logtheta, c, M, rdiag, eta0, sigma0, V,
                       eps.lambda, eps.phi) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])+eps.lambda
    phi <- exp(logtheta[length(logtheta)])+eps.phi
    logtheta <- c(log(lambda),log(phi))

    # Calculate first part of Q (from locally iid model)
    Q <- Q_iid(logtheta, c, M, rdiag)

    # Add log-normal prior component
    Q <- (Q - 1/2*t(logtheta-eta0) %*% chol2inv(chol(V+sigma0))
          %*% (logtheta-eta0))

    return(Q);
}

#' Compute analytic gradient of Q-function for smoothed EM algorithm of Cao et
#' al. (2000)
#'
#' Computes gradient of Q-function with respect to log(c(lambda,phi)) for EM
#' algorithm from Cao et al. (2000) for their smoothed model.
#'
#' @param logtheta numeric vector (length k+1) of log(lambda) (1:k) and log(phi)
#'      (last entry)
#' @param c power parameter in model of Cao et al. (2000)
#' @param M matrix (n x k) of conditional expectations for OD flows,
#'      one time per row
#' @param rdiag numeric vector (length k) containing diagonal of conditional
#'      covariance matrix R
#' @param eta0 numeric vector (length k+1) containing value for log(c(lambda,
#'      phi)) from previous time (or initial value)
#' @param sigma0 covariance matrix (k+1 x k+1) of log(c(lambda, phi)) from
#'      previous time (or initial value)
#' @param V evolution covariance matrix (k+1 x k+1) for log(c(lambda, phi))
#'      (random walk)
#' @param eps.lambda numeric small positive value to add to lambda for numerical
#'      stability; typically 0
#' @param eps.phi numeric small positive value to add to phi for numerical
#'      stability; typically 0
#' @return numeric vector of same length as logtheta containing calculated
#'      gradient
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
grad_smoothed <- function(logtheta, c, M, rdiag, eta0, sigma0, V,
                          eps.lambda, eps.phi) {
    # Get parameter values
    lambda <- exp(logtheta[-length(logtheta)])+eps.lambda
    phi <- exp(logtheta[length(logtheta)])+eps.phi
    logtheta <- c(log(lambda),log(phi))

    # Calculate first part of grad (from iid model)
    grad <- grad_iid(logtheta, c, M, rdiag)

    # Add gradient from normal prior
    grad <- grad - chol2inv(chol(V+sigma0)) %*% (logtheta - eta0)

    # Return gradient
    return(grad)
}

#' Run EM algorithm to obtain MLE (single time) for smoothed model of Cao et al.
#' (2000)
#'
#' Runs EM algorithm to compute MLE for the smoothed model of Cao et al.
#' (2000). Uses numerical optimization of Q-function for each M-step with
#' analytic computation of its gradient. This performs estimation for a single
#' time point using output from the previous one.
#'
#' @param Y matrix (h x k) of observations in local window; columns correspond
#'      to OD flows, and rows are individual observations
#' @param A routing matrix (m x k) for network being analyzed
#' @param eta0 numeric vector (length k+1) containing value for log(c(lambda,
#'      phi)) from previous time (or initial value)
#' @param sigma0 covariance matrix (k+1 x k+1) of log(c(lambda, phi)) from
#'      previous time (or initial value)
#' @param V evolution covariance matrix (k+1 x k+1) for log(c(lambda, phi))
#'      (random walk)
#' @param c power parameter in model of Cao et al. (2000)
#' @param maxiter maximum number of EM iterations to run
#' @param tol tolerance (in relative change in Q function value) for stopping EM
#'      iterations
#' @param eps.lambda numeric small positive value to add to lambda for numerical
#'      stability; typically 0
#' @param eps.phi numeric small positive value to add to phi for numerical
#'      stability; typically 0
#' @param method optimization method to use (in optim calls)
#' @return list with 5 elements: \code{lambda}, the estimated value of lambda;
#'      \code{phi}, the estimated value of phi; \code{iter}, the number of
#'      iterations run; \code{etat}, log(c(lambda, phi)); and sigmat, the
#'      inverse of the Q functions Hessian at its mode
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @export
#' @family CaoEtAl
smoothed_EM <- function(Y, A, eta0, sigma0, V, c=2, maxiter = 1e3, tol=1e-6,
                        eps.lambda=0, eps.phi=0, method="L-BFGS-B") {
    # Determine number of latent variables
    p <- ncol(A)

    # Setup initial parameters
    theta0 <- exp(eta0)

    lambda0 <- theta0[-length(eta0)] + 1e-2
    # lambda0 <- (lambda0 + lambda_init(Y))/2
    # lambda0 <- lambda_init(Y)
    lambda <- lambda0

    phi0 <- theta0[length(eta0)] + 1e-2
    phi <- phi0

    # Remove final column from Y to avoid colinearity
    Y <- as.matrix(Y[,-ncol(Y)])

    # Run first iteration
    # E step - calculation of m_t = E(x_t | y_t, \theta) and 
    #  R = Var(x_t | y_t, \theta)
    M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi,
                 A=A, c=c))
    R <- R_estep(lambda, phi=phi, A, c)

    # M step - Optimize Q wrt theta
    mstep <- optim(par=c(log(lambda0),log(phi0)), fn=function(...)
                   -Q_smoothed(...), gr=function(...) -grad_smoothed(...), c=c,
                   M=M, rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
                   eps.lambda=eps.lambda, eps.phi=eps.phi, method="BFGS")
    theta <- exp(mstep$par)
    lambda <- theta[-length(theta)]+eps.lambda
    phi <- theta[length(theta)]+eps.phi

    # Calculate initial ll
    ll <- -mstep$value

    # Run EM iterations
    for (iter in 1:maxiter) {
        # E step - calculation of m_t = E(x_t | y_t, \theta) and 
        #  R = Var(x_t | y_t, \theta)
        M <- t(apply(Y, 1, m_estep, lambda=lambda, phi=phi,
                     A=A, c=c))
        R <- R_estep(lambda, phi=phi, A, c)

        # M step - Optimize Q wrt theta
        mstep <- optim(par=c(log(lambda0),log(phi0)), fn=function(...)
                       -Q_smoothed(...), gr=function(...) -grad_smoothed(...),
                       c=c, M=M, rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
                       eps.lambda=eps.lambda, eps.phi=eps.phi, method="BFGS")
        theta <- exp(mstep$par)
        lambda <- theta[-length(theta)]+eps.lambda
        phi <- theta[length(theta)]+eps.phi
        ll.new <- -mstep$value

        # cat(ll, ll.new, abs(ll.new-ll)/abs(ll.new+ll)*2, "\n")

        # Check for convergence
        if (abs(ll.new-ll)/abs(ll.new+ll)*2 < tol) {
            ll <- ll.new
            break
        } else {
            ll <- ll.new
        }
    }
    mstep <- optim(par=log(theta), fn=function(...) -Q_smoothed(...),
                   gr=function(...) -grad_smoothed(...), c=c, M=M,
                   rdiag=diag(R), eta0=eta0, sigma0=sigma0, V=V,
                   eps.lambda=eps.lambda, eps.phi=eps.phi, method="BFGS",
                   hessian=TRUE)

    return(list(lambda=lambda, phi=phi, iter=iter, etat=c(log(lambda),log(phi)),
                sigmat=chol2inv(chol(mstep$hessian))))
}

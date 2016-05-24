##' @title (private) HOA test for regression coefficients in an NBP GLM model
##' @param y an n vector of counts
##' @param s an n vector of effective library sizes
##' @param x an n by p design matrix
##' @param phi an n vector of dispersion parameters
##' @param beta0 a p vector specifying null hypothesis: non NA
##' components are hypothesized values of beta, NA components are free
##' components
##' @param tol.mu convergence criteria
##' @param print.level a number, print level
##' @return test statistics and  p-values of HOA, LR, and Wald tests
hoa.hd  = function(y, s, x, phi, beta0, tol.mu=1e-3/length(y),
  print.level=1) {

  if (print.level>0)
    message("HOA test for regression coefficients.");

  ## NA values to be returned if things go wrong
  res.na =list(
    lambda.star=NA, pstar=NA,
    lambda=NA, p=NA,
    z.wald=NA, p.wald=NA,
    u=NA, p.score=NA);

  nobs = length(y);

  if (length(s)==1)  s = rep(s, nobs);

  if (length(phi)==1) phi = rep(phi, nobs);
  dim(phi)=NULL;

  ## The shape parameters 
  kappa = 1/phi;

  ## Index to the parameter of interest
  idh = !is.na(beta0);
  ## Dimension of the hypothesis
  nh = sum(idh);

  ## Indices to nuisance parameters
  idn = is.na(beta0);
  nn = sum(idn);

  ## log likelihood function
  el = function(mu, y) {
    sum(y*log(mu) - (y+kappa) * log(mu+kappa));
  }

  ## Find MLE under the full model
  fit = irls.nb.1(y, s, x, phi, tol.mu=tol.mu);

  ## If any mu is near 0, add 0.5 to the total count (proportionally
  ## distributed to columns)
  if (any(fit$mu<tol.mu*2)) {
    y = y + 0.5*s/sum(s);
    ## y = y + s/sum(s);
    fit = irls.nb.1(y, s, x, phi, tol.mu=tol.mu);
  }

  ## Likelihood quantities under the full model
  beta.hat = fit$beta;
  mu.hat = fit$mu;
  v.hat = mu.hat + phi * mu.hat^2;
  l.hat = el(mu.hat, y);
  i.hat = t(x) %*% diag(mu.hat^2/v.hat) %*% x;
  j.hat = t(x) %*% diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa)/(mu.hat+kappa)^2) - (y-mu.hat)*mu.hat/v.hat) %*% x;
  l.hat = el(mu.hat, y);

  if (print.level>2) {
    print(list(beta.hat=beta.hat, mu.hat=mu.hat, v.hat=v.hat,
               i.hat = i.hat, j.hat=j.hat, l.hat=l.hat));
  }

  ## Wald test statistic and p-value
  w = t(beta.hat[idh] - beta0[idh]) %*% solve(solve(i.hat)[idh, idh]) %*% (beta.hat[idh] - beta0[idh]);
  ## p.wald = 1 - pchisq(w, df=nh);
  p.wald = pchisq(w, df=nh, lower.tail=FALSE);
  
  ## Find MLE under the hypothesis
  fit0 = irls.nb.1(y, s, x, phi, beta0);

  ## Likelihood quantities under the null model
  beta.tilde = fit0$beta;
  mu.tilde = fit0$mu;
  v.tilde = mu.tilde + phi * mu.tilde^2;
  l.tilde = el(mu.tilde, y);
  i.tilde = t(x) %*% diag(mu.tilde^2/v.tilde) %*% x;
  j.tilde = (t(x) %*% diag(mu.tilde^2 * (y/mu.tilde^2 - (y+kappa)/(mu.tilde+kappa)^2) - (y-mu.tilde)*mu.tilde/v.tilde) %*% x);
  D1.tilde = t(x) %*% matrix(mu.tilde * (y-mu.tilde)/v.tilde, nobs, 1);

  if (print.level>2) {
    print(list(beta.tilde=beta.tilde, mu.tilde=mu.tilde, v.tilde=v.tilde,
               j.tilde=j.tilde, l.tilde=l.tilde));
  }

  ## Score test statistic and p-value
  i.tilde.inv = solve(i.tilde);
  u = t(D1.tilde) %*% i.tilde.inv %*% D1.tilde;
  p.score = pchisq(u, df=nh, lower.tail=FALSE);

  ## Likelihood ratio test statistic and p-value
  lambda = 2*(l.hat - l.tilde);
  ## p = 1 - pchisq(lambda, df=nh);
  p = pchisq(lambda, df=nh, lower.tail=FALSE);

  ## return(list(LRT=p, Wald=p.wald));

  if (l.hat <= l.tilde + sqrt(.Machine$double.eps)) {
    return(
  list(beta.hat = beta.hat,
       mu.hat = mu.hat,
       beta.tilde = beta.tilde,
       mu.tilde = mu.tilde,
       lambda.star=0, pstar=1,
       lambda = 0, p=1,
       w = 0, p.wald = 1,
       u = 0, p.score = 1));
  }

  ## Approximations to sample space derivatives
  S.hat = t(x) %*% diag(mu.tilde * mu.hat /v.tilde) %*% x;
  p.hat = mu.hat/(mu.hat+kappa);
  p.tilde = mu.tilde/(mu.tilde+kappa);
  q.hat = t(x)  %*% diag(mu.hat*log(p.hat/p.tilde)) %*% matrix(1, nobs, 1);

  ##
  if (print.level>2) {
    print(list(lambda=lambda, S.hat=S.hat, q.hat=q.hat));
  }

  S.hat.inv = solve(S.hat);
  i.hat.inv = solve(i.hat);
  j.hat.inv = solve(j.hat);
  
  j.tilde.nn = j.tilde[idn, idn];
  dim(j.tilde.nn) = c(nn,nn);

  I = sqrt(det(i.tilde))* sqrt(det(i.hat)) / det(S.hat) * sqrt(det(j.tilde.nn));

  j.tilde.tilde=i.tilde %*% S.hat.inv %*% j.hat %*% i.hat.inv %*% S.hat;
  j.tilde.tilde.nn = j.tilde.tilde[idn, idn];
  dim(j.tilde.tilde.nn) = c(nn, nn);

  II = det(j.tilde.tilde.nn)^{-0.5};

  III = (t(D1.tilde) %*% S.hat.inv %*% i.hat %*% j.hat.inv %*% S.hat %*% i.tilde.inv %*% D1.tilde)^(nh/2);
  IV = lambda^(nh/2 - 1) * (t(D1.tilde) %*% S.hat.inv %*% q.hat);

  adj = I * II * III / IV;

  if (print.level>2) {
    print(list(I=I, II=II, III=III, IV=IV, adj=adj));
  }

  lambda.star = lambda * ( 1 - 1/lambda * log(adj))^2;
  ## pstar = 1 - pchisq(lambda.star, df=nh);
  pstar = pchisq(lambda.star, df=nh, lower.tail=FALSE);
  
  list(beta.hat = beta.hat,
       mu.hat = mu.hat,
       beta.tilde = beta.tilde,
       mu.tilde = mu.tilde,
       lambda.star=lambda.star, pstar=pstar,
       lambda = lambda, p=p,
       w = w, p.wald = p.wald,
       u = u, p.score = p.score
       );
  ## list(r=r, rstar=rstar, p=p, pstar=pstar);
}


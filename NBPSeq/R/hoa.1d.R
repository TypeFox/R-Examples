## One-dimensional HOA test
## source('nbp-mle.R');
## source('nb.glm.R');

##' @title (private) One-dimensional HOA test for a regression coefficient in an NB GLM model
##' @param y an n vector of counts
##' @param s an n vector of effective library sizes
##' @param x an n by p design matrix
##' @param phi an n vector of dispersion parameters
##' @param beta0 a p vector specifying null hypothesis: non NA
##' components are hypothesized values of beta, NA components are free
##' components
##' @param tol.mu convergence criteria
##' @param alternative "less" means phi < 0.
##' @param print.level a number, print level
##' @return test statistics and p-values of HOA, LR, and Wald tests
hoa.1d = function(y, s, x, phi, beta0, tol.mu=1e-3/length(y),
  alternative="two.sided",
  print.level=1) {

  if (print.level>0)
    message("HOA test for regression coefficients.");

  ## NA values to be returned if things go wrong
  res.na =list(
    rstar=NA, pstar=NA,
    r=NA, p=NA,
    w=NA, p.wald=NA,
    u=NA, p.score=NA,
    NPadj=NA, INFadj=NA,
    alternative=alternative);

  ## Determine the alternative
  alt = pmatch(alternative, c("two.sided", "less", "greater"));
  if (is.na(alt)) {
    warning('alternative should be "two.sided", "less" or "greater"');
    return(res.na);
  }
  
  ## The shape (size) parameters 
  dim(phi)=NULL;
  kappa = 1/phi;

  nobs = length(y);
  if (length(s)==1) {
    s = rep(s, nobs);
  }

  ## Index to the parameter of interest
  idh = !is.na(beta0);
  ## Dimension of the hypothesis
  nh = sum(idh);

  if (sum(idh)!=1) stop("The hypothesis should be one dimensional.");

  ## Indices to nuisance parameters
  idn = is.na(beta0);
  nn = sum(idn);

  ## log likelihood function
  el = function(mu, y) {
    sum(y*log(mu) - (y+kappa) * log(mu+kappa));
  }

  ## Likelihood quantities under the full model
  ## fit = irls.nb.1(y, s, x, phi, tol.mu=1e-5);
  fit = irls.nb.1(y, s, x, phi, tol.mu=tol.mu);

  ## If any mu is near 0, add 0.5 to the total count (proportionally
  ## distributed to columns)
  if (any(fit$mu<tol.mu*2)) {
    y = y + 0.5*s/sum(s);
    ## y = y + s/sum(s);
    fit = irls.nb.1(y, s, x, phi, tol.mu=tol.mu);
  }

  beta.hat = fit$beta;
  mu.hat = fit$mu;
  v.hat = drop(mu.hat + phi * mu.hat^2);
  l.hat = el(mu.hat, y);
  i.hat = t(x) %*% diag(mu.hat^2/v.hat) %*% x;
  j.hat = t(x) %*% diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa)/(mu.hat+kappa)^2) - (y-mu.hat)*mu.hat/v.hat) %*% x;
  l.hat = el(mu.hat, y);

  if (print.level>2) {
    print(list(beta.hat=beta.hat, mu.hat=mu.hat, v.hat=v.hat,
               i.hat = i.hat, j.hat=j.hat, l.hat=l.hat));
  }

  ## Wald test statistic and p-value
  w = (beta.hat[idh] - beta0[idh]) / sqrt(solve(i.hat)[idh, idh]);

  ## Determine which tail to compute the probability for
  if (alt==2 | (alt==1 & w<0)){
    lower.tail=TRUE;
  } else {
    lower.tail=FALSE
  }

  p.wald =  pnorm(w, lower.tail=lower.tail);
  
  ## If estimate the variance under the null,
  ## w0 = t(beta.hat[idh] - beta0[idh]) %*% solve(solve(i.hat)[idh, idh]) %*% (beta.hat[idh] - beta0[idh]);
  ## p.wald0 = 1 - pchisq(w0, df=nh);

  ## Likelihood quantities under the null model
  ## fit0 = irls.nb.1(y, s, x, phi, beta0, tol.mu=1e-5);
  fit0 = irls.nb.1(y, s, x, phi, beta0, tol.mu=tol.mu);
  beta.tilde = fit0$beta;
  mu.tilde = fit0$mu;
  v.tilde = mu.tilde + phi * mu.tilde^2;
  l.tilde = el(mu.tilde, y);
  j.tilde = (t(x) %*% diag(mu.tilde^2 * (y/mu.tilde^2 - (y+kappa)/(mu.tilde+kappa)^2) - (y-mu.tilde)*mu.tilde/v.tilde) %*% x)[idn, idn];
  dim(j.tilde) = c(nn, nn);

  if (print.level>2) {
    print(list(beta.tilde=beta.tilde, mu.tilde=mu.tilde, v.tilde=v.tilde,
               j.tilde=j.tilde, l.tilde=l.tilde));
  }

  if (l.hat <= l.tilde + sqrt(.Machine$double.eps)) {
    ## This part can be improved, but it may not worth it
    return(
  list(beta.hat = beta.hat,
       mu.hat = mu.hat,
       beta.tilde = beta.tilde,
       mu.tilde = mu.tilde,
                r=0, rstar=0, rstar0=0, NPadj=0, INFadj=0, p=0.5, pstar=0.5,
                w=w, p.wald=p.wald));
  }

  ## Score test statistic and p-value
  ## D1 =
  ## u = t(D1) %*% solve(j.hat(idh, idh)) %*% D1;
  u = NA;
  p.score = NA;

  ## Deviance
  ## lambda = 2*(l.hat - l.tilde);

  ## Likelihood ratio test p-value

  ## Directed deviance
  r = sign(beta.hat[idh] - beta.tilde[idh])*sqrt(2*(l.hat - l.tilde));
  p = pnorm(r, lower.tail=lower.tail);

  ## Approximations to sample space derivatives
  S.hat = t(x) %*% diag(mu.tilde * mu.hat /v.tilde) %*% x;
  p.hat = mu.hat/(mu.hat+kappa);
  p.tilde = mu.tilde/(mu.tilde+kappa);
  q.hat = t(x)  %*% diag(mu.hat*log(p.hat/p.tilde)) %*% matrix(1, nobs, 1);

  ##
  if (print.level>2) {
    print(list(r=r, S.hat=S.hat, q.hat=q.hat));
  }

  ## If there are all 0's in one group, det(j.hat) can be negative!
  ## if (det(j.hat) <= 0){

  eps = 1e-2;
  if (abs(r) < eps){
    ## If r is too close to 0, INFadj will not be accurate (and there is no need for that)
    ## How about the NPadj
    z = 1;
    rstar = r;
  } else {
    z = r*det(i.hat)/sqrt(det(j.hat))/ det(S.hat) * sqrt(det(j.tilde)) / (solve(S.hat) %*% q.hat)[idh];
    rstar = r - (1/r) * log(z);
  }
    
  if (TRUE) {
    ## Compute nuisance parameter adjustment and NP adjustment
    ## See Pierce and Bellio (2006) for details

    ## Approximations to sample space derivatives
    el.th.thhat = S.hat %*% solve(i.hat) %*% j.hat;
    del.el.thhat =  t(q.hat) %*% solve(i.hat) %*% j.hat;

    ## Nuisance parameter adjustment
    el.nu.nu = el.th.thhat[idn, idn];
    dim(el.nu.nu) = c(nn, nn);
    
    C = det(el.nu.nu) / sqrt(det(j.tilde %*% j.hat[idn, idn]));
    NPadj = (1/r) * log(C);

    ## Information adjustemnt
    ## el.nu.thhat = t(matrix(el.th.thhat[idn,]));   # strange to need this
    elprof.psihat = del.el.thhat[idh] - el.th.thhat[idn, idh] %*% solve(el.nu.nu) %*% del.el.thhat[idn];
    info.hat.inv = solve(j.hat);
    info.prof = 1 / info.hat.inv[idh, idh];
    dim(info.prof) = c(nh, nh);
    
    util = elprof.psihat/sqrt(det(info.prof));
    INFadj = (1/r) * log(util/r);

    ## this is the final computation of rstar
    ## rstar0 should be identical to rstar comptued earlier
    rstar0 =r + NPadj + INFadj;

    if (print.level>2) {
      print(list(util=util, C = C, INFadj=INFadj, NPadj=NPadj, rstar0=rstar0));
    }

  }
  
  pstar = pnorm(rstar, lower.tail=lower.tail);

  if (alt==1) {
    ## Two-sided alternative
    pstar = min(pstar*2, 1);
    p = min(p*2, 1);
    p.wald = min(p.wald*2, 1);
    p.score = min(p.score*2, 1);
  }
  
  list(beta.hat = beta.hat,
       mu.hat = mu.hat,
       beta.tilde = beta.tilde,
       mu.tilde = mu.tilde,
       rstar=rstar, pstar=pstar,
       r = r, p=p,
       w = w, p.wald = p.wald,
##     u = u, p.score = p.score,
       NPadj=NPadj, INFadj=INFadj,
       alternative=alternative
       );
}


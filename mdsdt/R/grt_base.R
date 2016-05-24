#' @import ellipse mnormt polycor

#' @importFrom graphics abline box lines mtext par plot points text title
#' @importFrom stats coef constrOptim dnorm lm pchisq pnorm qnorm vcov xtabs

# S3 grt object
# 
# Constructor for a grt fit object, containing information about estimated parameters and likelihoods
#  dists: Matrix giving means, standard deviations, and correlation
#           for each stimuls type
#  fit:   Optional information from fitting program
#    obs:      Observed frequencies
#    fitted:   Fitted frequencies
#    estimate: Estimated parameter vector
#    expd2:    Hessian (expected second derivative) matrix at estimate
#    map:      Parameter map
#    loglik:   Log likelihood at estimate
#    code:     Convergence code (as nlm)
#    iter:     Number of iterations required
#  rcuts: Optional cutpoints for the rows
#  ccuts: Optional cutpoints for the columns
grt <- function (dists, fit=NULL, rcuts = 0, ccuts = 0) {
  if (is.null(colnames(dists))) {
    colnames(dists) <- c('mu_r','sd_r','mu_c','sd_c','rho')
  }
  structure(list(dists=dists,fit=fit, rowcuts=rcuts, colcuts=ccuts),
            class = 'grt')
}

#' Fit full Gaussian GRT model
#' 
#' Fit the mean and covariance of a bivariate Gaussian distribution for each stimulus class, subject to given constraints. 
#' Standard case uses confusion matrix from a 2x2 full-report identification experiment, but will also work in designs with N levels of confidence associated with each dimension (e.g. in Wickens, 1992).
#' 
#' @param freq Can be entered in two ways: 1) a 4x4 confusion matrix containing counts, 
#' with each row corresponding to a stimulus and each column corresponding to a response. 
#' row/col order must be a_1b_1, a_1b_2, a_2b_1, a_2b_2. 
#' 2) A three-way 'xtabs' table with the stimuli as the third index and the 
#' NxN possible responses as the first two indices.
#' @param PS_x if TRUE, will fit model with assumption of perceptual separability on the x dimension (FALSE by default)
#' @param PS_y if TRUE, will fit model with assumption of perceptual separability on the y dimension (FALSE by default)
#' @param PI 'none' by default, imposing no restrictions and fitting different correlations for all distributions. 
#' If 'same_rho', will constrain all distributions to have same correlation parameter. 
#' If 'all', will constain all distribution to have 0 correlation. 
#' @param method The optimization method used to fit the Gaussian model. Newton-Raphson gradient descent by default, but
#' may also specify any method available in \code{\link[stats]{optim}}, e.g. "BFGS".
#' @return An S3 \code{grt} object
#' @examples 
#' # Fit unconstrained model
#' data(thomas01b); 
#' grt_obj <- fit.grt(thomas01b);
#' 
#' # Use standard S3 generics to examine
#' print(grt_obj);
#' summary(grt_obj);
#' plot(grt_obj);
#' 
#' # Fit model with assumption of perceptual separability on both dimensions
#' grt_obj_PS <- fit.grt(thomas01b, PS_x = TRUE, PS_y = TRUE);
#' summary(grt_obj_PS);
#' plot(grt_obj_PS);
#' 
#' # Compare models 
#' GOF(grt_obj, teststat = 'AIC');
#' GOF(grt_obj_PS, teststat = 'AIC');
#' @export
fit.grt <- function(freq, PS_x = FALSE, PS_y = FALSE, PI = 'none', method=NA) {
  if (length(freq) == 16) {
    return(two_by_two_fit.grt(freq, PS_x, PS_y, PI))
  } else {
    n_by_nmap <- create_n_by_n_mod(freq,PS_x, PS_y, PI);
    return(n_by_n_fit.grt(freq, pmap=n_by_nmap, method=method))
  }
}

#' Print the object returned by fit.grt
#' @param x An object returned by fit.grt 
#' @param ... further arguments passed to or from other methods, as in the generic print function
#' @export
print.grt <- function (x, ...) {
  cat('Row cuts:',format(x$rowcuts,digits=3),'\n')
  cat('Col cuts:',format(x$colcuts,digits=3),'\n')
  cat('Distributions:\n')
  print(round(x$dists,3))
  invisible(x)
}

#' Summarize the object returned by fit.grt
#' @param object An object returned by fit.grt
#' @param ... additional arguments affecting the summary produced, as in the generic summary function
#' @export
summary.grt <- function(object, ...) {
  print.grt(object)
  if (!is.null(fit <- object$fit)){
    cat('\n')
    cat('Standard errors:\n')
    print(round(distribution.se(object),3))
    cat('\n')
    cat('Fit statistics:\n')
    cat(c('Log likelihood: ',round(fit$loglik,2),'\n'))
    cat(c('AIC: ',object$AIC,'\n'))
    cat(c('AIC.c: ',object$AIC.c,'\n'))
    cat(c('BIC: ',object$BIC,'\n'))
    cat('\nConvergence code',fit$code,
        'in',fit$iter,'iterations\n')
  }
  invisible(object)
}

#' Plot the object returned by fit.grt
#' 
#' @param x a grt object, as returned by fit.grt
#' @param level number between 0 and 1 indicating which contour to plot (defaults to .5)
#' @param xlab optional label for the x axis (defaults to NULL)
#' @param ylab optional label for the y axis (defaults to NULL)
#' @param marginals Boolean indicating whether or not to plot marginals (only available for 2x2 model; defaults to FALSE)
#' @param main string to use as title of plot (defaults to empty string)
#' @param plot.mu Boolean indicating whether or not to plot means (defaults to T)
#' @param ... Arguments to be passed to methods, as in generic plot function
#' @export
plot.grt <- function(x, level = .5, xlab=NULL, ylab=NULL, marginals=F, main = "", plot.mu=T,...) {
#                     connect=NULL, names=NULL, clty=1,ccol='Black',llty=1,lcol='Black', ...) {
  origPar <- par(no.readonly=TRUE); 
  lim.sc=1
  connect=NULL;
  names=NULL;
  clty=1;
  ccol='Black';
  llty=1;
  lcol='Black';
  if (length(x$fit$obs) == 16) {
    two_by_two_plot.grt(x, xlab, ylab, level = level, 
                        marginals=marginals, main = main, 
                        plot.mu = plot.mu);
  } else {
    #require(mvtnorm);
    #main=deparse(substitute(x))
    xc <- x$colcuts
    yc <- x$rowcuts
    dd <- x$dists
    mx <- dd[,3]; my <- dd[,1]
    sx <- dd[,4]; sy <- dd[,2]
    rho <- dd[,5]
    min.ax <- min(min(mx-lim.sc*sx),min(my-lim.sc*sy))
    max.ax <- max(max(mx+lim.sc*sx),max(my+lim.sc*sy))
    X <- c(min.ax,max.ax)
    Y <- c(min.ax,max.ax)
    if (is.null(xlab)) xlab <- if(is.null(x$fit)) 'A' else
      names(dimnames(x$fit$obs)[2])
    if (is.null(ylab)) ylab <- if(is.null(x$fit)) 'B' else
      names(dimnames(x$fit$obs)[1])
    
    # axes=F, box(which="plot") added 1.24.14
    par(fig=c(0,1,0,1),mar=c(2.5,2.5,2.5,2.5))
    plot(X,Y,type='n',main=main,xlab="",ylab="",axes=F,...)
    mtext(text=xlab,side=1,line=1)
    mtext(text=ylab,side=2,line=1)
    box(which="plot")
    for (i in 1:length(yc)) abline(h=yc[i],lty=llty,col=lcol)
    for (i in 1:length(xc)) abline(v=xc[i],lty=llty,col=lcol)
    for (i in 1:dim(dd)[1]) {
      v <- matrix(c(sx[i]^2,rep(sx[i]*sy[i]*rho[i],2),sy[i]^2),2)
      lines(ellipse(v,centre=c(mx[i],my[i]),level=level))
      if(plot.mu){
        points(mx[i],my[i],pch=3)
      }
    } 
    if (!is.null(connect[1])) 
      lines(mx[c(connect,connect[1])],my[c(connect,connect[1])],
            lty=clty,col=ccol)
    # names changed 1.24.14
    if (!is.null(names)) text(mx,my,names)
  }
  # Restore user's original graphics params
  par(origPar);
}

#' Test report independence 
#' 
#' Test report independence for each stimulus response distribution
#'
#' @param x four-by-four confusion matrix 
#' @return data frame containing z-scores and p-values for all four tests
#' @details If p value is sufficiently low, we're justified in rejecting the null hypothesis of sampling within that factor. p values come from a chi-squared test on the confusion matrix, as explaned in a footnote of Thomas (2001).
#' @examples
#' data(thomas01a)
#' riTest(thomas01a)
#' @source
#' Ashby, F. G., & Townsend, J. T. (1986). Varieties of perceptual independence. Psychological review, 93(2), 154.
#'
#' Thomas, R. D. (2001).Perceptual interactions of facial dimensions in speeded classification and identification. Perception \& Psychophysics, 63(4), 625--650.
#'
#' Silbert, N. H., & Thomas, R. D. (2013). Decisional separability, model identification, and statistical inference in the general recognition theory framework. Psychonomic bulletin & review, 20(1), 1-20.
#' @export
riTest <- function(x) {
  if(!checkConfusionMatrix(x)) {
    return(FALSE)
  }
  
  stimulus <- c("(A1,B1)", "(A1,B2)", "(A2,B1)", "(A2,B2)")
  statistic <- rep(NA,4)
  for ( i in 1:4 ) { 
    x1 <- matrix(x[i,], 2,2, byrow=T)
    ex1 <- c(apply(x1,1,sum)*apply(x1,2,sum),
             rev(apply(x1,1,sum)) * apply(x1,2,sum))
    ex1 <- matrix(ex1[c(1,4,3,2)],2,2,byrow=TRUE) / sum(x1)
    statistic[i] <- sum( (x1 - ex1)^2/ ex1 )
  }
  
  return(data.frame(stimulus=stimulus, chi.2=round(statistic,3), 
                    p.value=round(1-pchisq(statistic, 1),3)))
}

#' Test marginal response invariance
#' 
#' Tests marginal response invariance at both levels on each dimension
#'
#' @param x four-by-four confusion matrix 
#' @return data frame containing z-scores and p-values for all four tests
#' @details If the p value for either level of the x dimension is significant, 
#' we are justified in rejecting the null hypothesis of perceptual separability on the x dimension. 
#' Similarly for the y dimension. 
#' 
#' The estimator is derived in a footnote of Thomas (2001).
#' @examples
#' data(thomas01a)
#' mriTest(thomas01a)
#' @source
#' Ashby, F. G., & Townsend, J. T. (1986). Varieties of perceptual independence. Psychological review, 93(2), 154.
#'
#' Thomas, R. D. (2001).Perceptual interactions of facial dimensions in speeded classification and identification. Perception \& Psychophysics, 63(4), 625--650.
#'
#' Silbert, N. H., & Thomas, R. D. (2013). Decisional separability, model identification, and statistical inference in the general recognition theory framework. Psychonomic bulletin & review, 20(1), 1-20.
#' @export
mriTest <- function(x) {
  if(!checkConfusionMatrix(x)) {
    return(FALSE)
  }
  
  response <- c("(A1,-)", "(A2,-)", "(-,B1)", "(-,B2)")
  statistic <- rep(NA,4)
  
  for ( A in 1:2 ) {
    rw <- 2*(A-1)+1
    rA.sAB1 <- sum( x[rw,  rw:(rw+1)] )
    rA.sAB2 <- sum( x[rw+1,rw:(rw+1)] )
    nAB1 <- sum(x[rw,])
    nAB2 <- sum(x[rw+1,])
    
    p.s <- (rA.sAB1 + rA.sAB2)/(nAB1 + nAB2)
    statistic[A] <- ((rA.sAB1/nAB1 - rA.sAB2/nAB2)/
                       sqrt(p.s*(1-p.s)*(1/nAB1+1/nAB2)) )
  }
  
  for ( B in 1:2 ) {
    rB.sA1B <- sum( x[B,c(B,B+2)] )
    rB.sA2B <- sum(x[B+2,c(B,B+2)] )
    nA1B <- sum(x[B,])
    nA2B <- sum(x[B+2,])
    
    p.s <- (rB.sA1B + rB.sA2B)/(nA1B + nA2B)
    statistic[B+2] <- ((rB.sA1B/nA2B - rB.sA2B/nA1B)/
                         sqrt(p.s*(1-p.s)*(1/nA1B+1/nA2B)) )
  }
  return(data.frame(response=response, z=round(statistic,3), 
                    p.value= round(2*(pmin(1-pnorm(statistic),pnorm(statistic))),3) ))
}


two_by_two_fit.grt <- function(freq, PS_x = FALSE, PS_y = FALSE, PI = 'none') {
  if(!checkConfusionMatrix(freq)) return(FALSE); # Make sure confusion matrix valid
  delta <- 1/10000; # Tolerance
  freq = freq + 1;   # protection against zeros, which cause the algorithm to explode
  w = 1;  # initialize weight for adjusting step size/preventing oscillation
  
  # initialize predicted probability matrix
  prob <- matrix(data=0, nrow = 4, ncol = 4)
  
  # use observed relative frequencies as first estimates of probs
  for (i in 1:4) {
    prob[i,] = freq[i,]/sum(freq[i,]);
  }
  
  # Get good initial estimates
  initial = initial_point(prob, PS_x, PS_y, PI);
  xpar=initial$xpar; ypar=initial$ypar; rpar=initial$rpar; 
  ps_old=initial$ps_old; rows=initial$rows; npar = length(xpar)+length(ypar)+length(rpar);
  d = matrix(data = 0, nrow = npar, ncol = 1); # Store gradient at estimate
  v <- array(0, dim = c(4,4,3)); # Store estimate variances
  E = matrix(data = 0, nrow = npar, ncol = npar); # Store information matrix
  
  # calculate prob estimates and co-variances for first step
  temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_old);
  prob = temp$prob;
  v = temp$v;
  
  # initial computation of log-likelihood gradient d & information matrix E
  for (i in 1:npar) {
    g <- if(sum(i==xpar)) 1 else (if(sum(i==ypar)) 2 else 3);
    ir <- rows[i,];
    d[i] = sum(sum((freq[ir,]/prob[ir,])
                   *v[ir,, g]));   
    for (j in 1:npar) {
      h <- if(sum(j==xpar)) 1 else (if(sum(j==ypar)) 2 else 3);
      jr <- rows[j,];
      ijr <- ir == 1 & jr == 1;
      if (sum(ijr) > 0) {
        sum_f <- if(is.null(dim(ijr))) sum else rowSums; # rowSums doesn't reduce to sum when d=1...
        K = sum_f(freq[ijr,]) / sum_f(prob[ijr,]);
        L = (sum_f(v[ijr,,g] * v[ijr,,h] / prob[ijr,]) -
               sum_f(v[ijr,,g])* sum_f(v[ijr,,h]) / sum_f(prob[ijr,]));
        E[i,j] = -t(K)*L; 
      }
    }
  }
  Ei = solve(E, diag(npar));  
  ps_new = ps_old - Ei %*% d;
  # iterate!
  
  it = 1;
  df = abs( ps_new - ps_old ) / abs(ps_new);
  dfp_new = t(df) %*% df;
  dfp_old = dfp_new;
  while (dfp_new > delta) {    
    ps_old = ps_new;
    
    # calculate prob estimates and co-variances
    temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_old);
    prob = temp$prob;
    v = temp$v;
    
    # log-likelihood gradient, information matrix
    for (i in 1:npar) {
      g <- if(sum(i==xpar)) 1 else (if(sum(i==ypar)) 2 else 3);
      ir <- rows[i,];
      d[i] = sum(sum((freq[ir,]/prob[ir,])
                     *v[ir,, g]));   
      for (j in 1:npar) {
        h <- if(sum(j==xpar)) 1 else (if(sum(j==ypar)) 2 else 3);
        jr <- rows[j,];
        ijr <- ir == 1 & jr == 1;
        if (sum(ijr) > 0) {
          sum_f <- if(is.null(dim(ijr))) sum else rowSums; # rowSums doesn't reduce to sum when d=1...
          K = sum_f(freq[ijr,]) / sum_f(prob[ijr,]);
          L = (sum_f(v[ijr,,g] * v[ijr,,h] / prob[ijr,]) -
                 sum_f(v[ijr,,g])* sum_f(v[ijr,,h]) / sum_f(prob[ijr,]));
          E[i,j] = -t(K)*L; 
        }
      }
    }
    Ei = solve(E, diag(npar));
    
    if (dfp_new > dfp_old) { #w halving procedure
      w = .5*w;
    }
    
    ps_new = ps_old - w * Ei %*% d;
    df = abs(ps_new-ps_old) / abs(ps_new);
    dfp_old = dfp_new;
    dfp_new = t(df) %*% df;
    it = it + 1;
  }
  
  parameters <- make_parameter_mat(xpar, ypar, rpar, ps_new);
  temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_new);
  prob = temp$prob; 
  
  # calculate various fit statistics and put them in output structure
  loglike <- sum(sum(freq * log(prob)));
  info_mat <- -solve(E, diag(npar));  
  nll = -loglike;
  aic = 2*npar - 2*loglike;
  bic = npar*log(sum(freq)) - 2*loglike;
  icomp = -loglike + (npar/2)*log(tr(info_mat)/npar)- .5*log(det(info_mat));
  fit <- list(obs=freq2xtabs(freq),fitted=freq2xtabs(prob), estimate=ps_new,
              expd2=E, map=create_n_by_n_mod(freq, PS_x, PS_y, PI, from_2x2 = TRUE), iter=it, 
              loglik=nll);#, aic = aic, bic = bic, icomp = icomp)
  output = grt(parameters, fit, 0, 0)
  output[['AIC']] = GOF(output,'AIC')
  output[['AIC.c']] = GOF(output,'AIC.c')
  output[['BIC']] = GOF(output,'BIC')
  return(output)
}



two_by_two_plot.grt <- function(obj, xlab, ylab, level = .5, 
                                marginals=F, main = "", plot.mu = T) {
  bin_width= .05; # determines smoothness of marginal plots
  fit_params = get_fit_params(obj)
  xlims = c(min(c(fit_params[[1]][1],fit_params[[2]][1],
                  fit_params[[3]][1],fit_params[[4]][1])-2.5),
            max(c(fit_params[[1]][1],fit_params[[2]][1],
                  fit_params[[3]][1],fit_params[[4]][1])+2.5))
  ylims = c(min(c(fit_params[[1]][2],fit_params[[2]][2],
                  fit_params[[3]][2],fit_params[[4]][2])-2.5),
            max(c(fit_params[[1]][2],fit_params[[2]][2],
                  fit_params[[3]][2],fit_params[[4]][2])+2.5))
  
  x = seq(xlims[1],xlims[2],by=bin_width)
  y = seq(ylims[1],ylims[2],by=bin_width)
  xra = xlims[2] - xlims[1]
  yra = ylims[2] - ylims[1]
  if (is.null(xlab)) {
    xlab <- "A";
  }
  if (is.null(ylab)) {
    ylab <- "B";
  }

  if(marginals){
    ex = .25
    xlb = ""
    ylb = ""
    par(fig=c(ex,1,ex,1), mar=c(.05, .05, 1.25, 1.25),pty="m",xaxt="n",yaxt="n")
  }else{
    xlb = xlab
    ylb = ylab
    par(fig=c(0,1,0,1), mar=c(2.5,2.5,2.5,2.5),pty="m",xaxt="n",yaxt="n")
  }
  # Make plotting frame
  plot(0,0,type="n",xlim=xlims,ylim=ylims,xlab=xlb,ylab=ylb,main=main,axes=F)
  box(which="plot",mai=rep(0,4))
  # Plot decision bounds
  abline(h=0);
  abline(v=0);
  # Plot Gaussian contours 
  for (i in 1:4) {
    cond = fit_params[[i]]
    cov <- matrix(data=c(1, cond[3], cond[3], 1), nrow = 2)
    mu = cond[1:2]
    par(new = TRUE)
    lines(ellipse(cov,centre=mu,level=level))
    if(plot.mu){
      points(cond[1], cond[2], pch = '+')
    }
    title(xlab = xlb, ylab = ylb)
  }
  
  # Add labels inset at 10% of the total x and y range
  labs = dimnames(obj$fit$obs)$Stim

  # If we're using the default, we need to get subscripts in there...
  newLabs <- vector("expression", 4)
  if(isDefaultLabel(labs)) {
    for(i in 1:4){
      first = strsplit(substring(labs[i], 0, 3), "_")[[1]];
      second = strsplit(substring(labs[i], 4, 6), "_")[[1]];
      firstIndex = as.numeric(first[2]);
      secondIndex = as.numeric(second[2]);
      exp = eval(bquote(expression(.(first[1])[.(firstIndex)]~.(second[1])[.(secondIndex)])))
      newLabs[i] = exp;
    }
  } else {
    newLabs = labs;
  }
  text(xlims[1]+(xra * .15), ylims[1]+(yra * .1), newLabs[1])
  text(xlims[1]+(xra * .15), ylims[2]-(yra * .1), newLabs[2])
  text(xlims[2]-(xra * .15), ylims[1]+(yra * .1), newLabs[3])
  text(xlims[2]-(xra * .15), ylims[2]-(yra * .1), newLabs[4])

  if(marginals){
    # compute marginals
    margx = margy = list(aa=NULL,ab=NULL,ba=NULL, bb=NULL);
    for (i in 1:4) {
      cond = fit_params[[i]];
      margx[[i]] = (1 / sqrt(2*pi)) * exp(-.5*((x-cond[1]))^2);
      margy[[i]] = (1 / sqrt(2*pi)) * exp(-.5*((y-cond[2]))^2);
    }
    
    # Plot X marginals
    par(fig=c(ex,1,0,ex), mar=c(.05, .05, 0.05, 1.25), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
    plot(x,margx$aa,type='l', lty=1, xlab = xlab, ylab = NULL);
    lines(x,margx$ab,type='l',lty=2);
    lines(x,margx$ba,type='l',lty=1);
    lines(x,margx$bb,type='l',lty=2);
    
    # Plot Y marginals
    par(fig=c(0,ex,ex,1), mar=c(.05, .05, 1.25, 0.05), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
    plot(margy$aa,y,type='l', lty=1, xlab = NULL, ylab = ylab);
    lines(margy$ab,y, type='l',lty=1);
    lines(margy$ba,y,type='l',lty=2);
    lines(margy$bb,y,type='l',lty=2);
  }
  
}

isDefaultLabel <- function(label) {
  return(label[1] == "a_1b_1" 
         & label[2] == "a_1b_2" 
         & label[3] == "a_2b_1" 
         & label[4] == "a_2b_2")
}

#' Conduct goodness of fit tests
#' 
#' Includes a number of common goodness of fit measures to compare different GRT models.
#' 
#' @param grtMod a \code{grt} object
#' @param teststat a string indicating which statistic to use in the test. 
#' May be one of the following:
#' \itemize{
#' \item{'X2'}{for a chi-squared test} 
#' \item{'G2'}{for a likelihood-ratio G-test}
#' \item{'AIC'}{for Akaike information criterion score}
#' \item{'AIC.c'}{for the AIC with finite sample size correction}
#' \item{'BIC'}{for Bayesian information criterion score}}
#' @param observed optional, to provide a matrix of observed frequencies if no fit conducted.
#' @examples 
#' data(thomas01a)
#' fit1 <- fit.grt(thomas01a)
#' fit2 <- fit.grt(thomas01a, PI = 'same_rho')
#' 
#' # Take the model with the lower AIC
#' GOF(fit1, teststat = 'AIC')
#' GOF(fit2, teststat = 'AIC')
#' @export
GOF <- function(grtMod,teststat='X2',observed=NULL){
  if (!identical(class(grtMod),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- grtMod$fit) && is.null(observed))
    stop('Must have fitted model, observed frequencies, or both')
  # added AIC, AICc, BIC 1.27.14
  statlist <- c('X2','G2','AIC','AIC.c','BIC')
  #print(c(statlist,teststat))
  #print(statlist %in% teststat)
  #teststat <- toupper(teststat)
  test <- pmatch(teststat,statlist)
  if (is.na(test)) stop('Test statistic unrecognized')
  teststat <- statlist[test]
  df <- if (is.null(observed)) length(ff$estimate) else 0
  if (is.null(observed)) observed <- ff$obs
  if (!is.null(ff)) ex <- ff$fitted else{
    nk <- apply(observed,3,sum)
    ex <- array(dim=dim(observed))
    for (k in 1:dim(observed)[3])
      ex[,,k] <- bsd.freq(grtMod$rowcuts,grtMod$colcuts,grtMod$dists[k,],nk[k])
  }
  df <- length(observed) - dim(observed)[3] - df
  if (test == 1){
    if(df == 0) {
      warning("degrees of freedom = 0, so the resulting statistic is uninterpretable. Try using AIC or BIC instead.")
    }
    tstat <- sum((observed-ex)^2/ex)}
  if (test == 2){
    ex <- ex[observed>0]; observed <- observed[observed>0]
    tstat <- 2*sum(observed*log(observed/ex))
  }
  if (test == 3){
    map = grtMod$fit$map
    k = 0
    for(i in 1:ncol(map)){
      k = k + sum(unique(map[,i])>0)
    }
    tstat <- 2*grtMod$fit$loglik + 2*k
  }
  if (test == 4){
    map = grtMod$fit$map
    k = 0
    for(i in 1:ncol(map)){
      k = k + sum(unique(map[,i])>0)
    }
    n <- sum(observed)
    aicStat <- 2*grtMod$fit$loglik + 2*k
    tstat <- aicStat + (2*k*(k+1))/(n-k-1)
  }
  if (test == 5){
    map = grtMod$fit$map
    k = 0
    for(i in 1:ncol(map)){
      k = k + sum(unique(map[,i])>0)
    }
    n <- sum(observed)
    tstat <- 2*grtMod$fit$loglik + log(n)*k  }
  if (test < 3){
    structure(tstat,names=teststat,df=df,class='bsdGOF')    
  }else{
    structure(round(tstat,1),names=teststat)
  } 
}

# Overall fitting function.  Fitting either uses Newton-Raphson iteration
# or one of the method provided by the R function constrOptim 
# xx: The frequency table.  It can be entered in two ways
#     1)  A three-way 'xtabs' table with the stimulis as the third index
#     2)  A data frame contiaing the three indices with condition last,
#         and frequencies as the variable 'x' (see 'form' if not this way)
# pmap:     Parameter-mapping array (default: complete parameterization)
# form:     A formula to convert a data frame (default x~.)
# p0:       Initial point---calculated if not given
# verbose:  Whether to print results (default FALSE)
# trace:    Print steps in minimization (default FALSE)
# method:   NA for method of scoring (default) or a method used by
#           constrOptim
# maxiter:  Maximum number of iterations
# stepsize: Size of iteration step (usually 1) in method of scoring
# maxch:    Maximum change for Newton Raphson convergence
# ...:      Arguments passed to minimization or likelihood routines
# Returns a grt object
n_by_n_fit.grt <- function (xx, pmap=NA, formula=x~., p0=NA, method=NA,
                        verbose=FALSE, trace=FALSE, maxiter=100, stepsize=1.0, maxch=1.0e-5, ...) {
  if (identical(class(xx)[1],'data.frame')) xx <- xtabs(formula,xx)
  if (!any(class(xx) == 'table'))
    stop('First argument must be data frame or contingency table')
  #  if (!(is.na(method) || (method == 0) || (method == 1)))
  #     stop('Method must be NA, 0, or 1')
  xx = xx + 1 # protection against zeros, which causes problems with the algorithm
  dxx <- dim(xx)
  if (length(dxx) != 3) stop('Table must have dimension 3')
  KK <- dxx[3];
  if (is.na(pmap)[1]) pmap <- matrix(c(rep(0:(KK-1),4),1:KK),4)
  colnames(pmap) <- c('mu','sigma','nu','tau','rho')
  if (is.null(rownames(pmap))) rownames(pmap) <- dimnames(xx)[[3]]
  bsd.valid.map(pmap,KK)
  # Create initial vector if required
  if (is.na(p0[1])) p0 <- bsd.initial(xx,pmap)
  # Construct index arrays
  imap <- bsd.imap(pmap,dxx)
  if (verbose) {
    cat('Parameter mapping vector\n'); print(pmap)
    cat('Initial parameter vector\n'); # print(round(p0,3))
    cat('Row cutpoints', round(p0[imap$xi],4),'\n')
    cat('Col cutpoints', round(p0[imap$eta],4),'\n')
    cat('Parameters by groups\n')
    print(round(bsd.map2array(p0,pmap,imap),4))
  }
  # Do the minimization
  if (is.na(method)){
    if (verbose) cat('Fitting by Newton-Raphson iteration\n')
    found <- FALSE
    pold <- p0
    for (iter in 1:maxiter) {
      #print(pold)
      #print(pmap)
      #print(imap)
      grtMod <- bsd.llike(pold,xx,pmap,imap,d.level=2,...)
      #print(attr(grtMod,'ExpD2'))
      #print(attr(grtMod,'gradient'))
      dlt <- solve(attr(grtMod,'ExpD2'),attr(grtMod,'gradient'))
      if (trace){
        cat('Iteration number',iter,'\n')
        cat('Row cutpoints', round(pold[imap$xi],4),'\n')
        cat('Col cutpoints', round(pold[imap$eta],4),'\n')
        print(round(bsd.map2array(pold,pmap,imap),4))
        cat('Value',grtMod[1],'\n')
      }
      s <- stepsize
      repeat{
        pp <- pold + s*dlt
        if (all(c(diff(pp[imap$xi]),diff(pp[imap$eta]))>0)) break
        s <- s/2
        warning('Reduced stepsize to',s,call.=FALSE)
        if (s < 0.001)
          stop('Stepsize reduction too large: check problem definition')
      }
      if (max(abs(dlt)) < maxch) {
        found <- TRUE
        iterations <- iter
        break
      }
      else iterations <- maxiter
      pold <- pp
    }
    code <- if (found) 0 else 4
  }
  else {
    if (verbose) cat('Fitting using optim\n')
    ixxi <- imap$xi; lxi <- length(ixxi)
    ixeta <- imap$eta; leta <- length(ixeta)
    ixvar <- c(imap$sigma,imap$tau) ; lvar <- length(ixvar)
    ixrho <- imap$rho; lrho <- length(ixrho)
    lp <- length(p0)
    if(length(ixxi)>1 && length(ixeta)>1){
      nconst <- lxi+leta+lvar+2*lrho-2
    }else{
      nconst <- 2*lrho
      b <- 1
    }
    cm <- matrix(0,nrow=nconst,ncol=lp)
    if(length(ixxi)>1){
      for(i in 1:(lxi-1)) cm[i,ixxi[i:(i+1)]] <- c(-1,1)
      b <- lxi-1
    }
    if(length(ixeta)>1){
      for(j in 1:(leta-1)) cm[j+b,ixeta[j:(j+1)]] <- c(-1,1)
      b <- b + leta-1
    }
    if (lvar > 0){
      for (i in 1:lvar) cm[i+b,ixvar[i]] <- 1
      b <- b + lvar
    }
    if (lrho > 0) for (i in 1:lrho){
      if(length(ixxi)>1 && length(ixeta)>1){
        cm[(b+1):(b+2),ixrho[i]] <- c(1,-1)
        b <- b+2
      }else{
        cm[b:(b+1),ixrho[i]] <- c(1,-1)
        b <- b+2
      }
    }
    cv <- c(rep(0,nconst-2*lrho),rep(-1,2*lrho))
    ffit <- constrOptim(p0,bsd.like,bsd.grad,cm,cv,method=method,
                        x=xx, pmap=pmap, imap=imap,...)
    pp <- ffit$par
    code <- ffit$convergence
    iterations <- ffit$counts
    grtMod <- bsd.llike(pp,xx,pmap,imap,d.level=2)
    found <- code < 4
  }
  if (!found) warning('Minimization routine encountered difficultites',
                      call.=FALSE)
  # Assemble results and print if required
  names(pp) <- names(p0)
  xi <- pp[imap$xi]
  eta <- pp[imap$eta]
  dists <- bsd.map2array(pp,pmap,imap)
  nk <- apply(xx,3,sum)
  muh <- array(dim=dim(xx),dimnames=dimnames(xx))
  for (k in 1:KK) muh[,,k] <- bsd.freq(xi,eta,dists[k,],nk[k])
  fit <- list(obs=xx,fitted=muh,estimate=pp,
              expd2=attr(grtMod,'ExpD2'),map=pmap,
              loglik=grtMod[1],code=code,iter=iterations)
  if (verbose) {
    if (found) cat('\nConvergence required',iterations,'iterations\n\n')
    else cat (iterations, 'iterations used without reaching convergence\n\n')
    cat('Parameter estimate vector\n'); # print(round(pp,3))
    cat('Row cutpoints', round(xi,4),'\n')
    cat('Col cutpoints', round(eta,4),'\n')
    print(round(dists,4))
    cat('Log likelihood =',grtMod[1],'\n')
  }
  #print(dists);
  output = grt(dists,fit=fit,rcuts = xi,ccuts = eta)
  output[['AIC']] = GOF(output,'AIC')
  output[['AIC.c']] = GOF(output,'AIC.c')
  output[['BIC']] = GOF(output,'BIC')
  return(output)
}

create_n_by_n_mod <- function(freq=NULL, PS_x=F, PS_y=F, PI="none", from_2x2 = FALSE) {
  # Each row is distribution, cols are y_mean, y_std, x_mean, x_std, rho
  if(from_2x2 | length(freq)==16){
    nst = 4
  }else{
    nst = dim(freq)[3]
  }
  # number of stimulus levels per dimension
  # assumes equal numbers of levels on each dimension
  nspd = sqrt(nst) 
  map <- matrix(data = 0, nrow = nst, ncol= 5)
  if (PS_y) {
    for(si in seq(2,nspd)){
      for(ri in seq(si,nst,nspd)){
        map[ri,1:2] <- c(si-1,si-1)
      }
    }
  } else for (i in 1:nst) map[i,1:2] <- c(i-1,i-1);
  if (PS_x) {
    ci = seq(nspd+1,nst,nspd)
    for(si in seq(1,nspd-1)){
      qi = ci[si]
      for(ri in seq(qi,qi+nspd-1)){
        map[ri,3:4] <- c(si,si)
      }
    }
  } else for (i in 1:nst) map[i,3:4] <- c(i-1,i-1);
  if (PI == 'same_rho') {
    for (i in 1:nst) map[i,5] <- 1;
  } else if (PI == 'none') {
    for (i in 1:nst) map[i,5] <- i;
  } 
  if (from_2x2) {
    map[,c(1,3)] = map[,c(1,3)] + 1;
    map[,c(2,4)] = c(0,0,0,0);
  }
  colnames(map) <- c("nu","tau","mu","sigma","rho")
  #print(map)
  return(map)
}


# Get parameter map
parameter.map <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- bb$fit)) NULL else ff$map
}


# Estimated parameters
coef.grt <- function(bb)
  if (is.null(ff <- bb$fit)) NULL else ff$estimate


# Covariance matrix of parameters
vcov.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  vcv <- solve(-ff$expd2)
  rownames(vcv) <- colnames(vcv) <- names(ff$estimate)
  vcv
}


# Parameters by stimuli
distribution.parameters <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  bb$dists
}


# Standard errors by stimuli
distribution.se <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- bb$fit)) return(NULL)
  pmap <- ff$map
  dimx <- dim(ff$obs)
  imap <- bsd.imap(pmap,dimx)
  sds <- bsd.map2array(sqrt(diag(vcov(bb))),pmap,imap,0,0)
  rownames(sds) <- rownames(bb$dists)
  colnames(sds) <- colnames(bb$dists)
  sds
}


# Log likelihood
logLik.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  sum(lfactorial(apply(ff$obs,3,sum))) - sum(lfactorial(ff$obs)) +
    ff$loglik
}


# Pearson residuals
residuals.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  xx <- ff$obs
  ex <- ff$fitted
  (xx-ex)/sqrt(ex)
}


# Fitted values
fitted.grt <- function(bb)
  if (is.null(ff <- bb$fit)) NULL else ff$fitted



print.bsdGOF <- function(gof){
  df <- attr(gof,'df')
  cat(names(gof),'(',df,') = ',gof,', p = ',
      round(pchisq(gof,df,lower.tail=F),5), '\n',sep='')
}


# Wald test of a linear hypothesis m p = c 
#   b:   fitted grt model containing estimates of p
#   m:   contrast vector or matrix with contrast vectors as rows
#   c:   a vector of numerical values (default zero)
#   set: set of parameters to test (means, sds, correlations,
#           distribution parameters, cutpoints, or all parameters)
linear.hypothesis <- function(b,m,c=0,set='means'){
  if (is.null(ff <- b$fit)) stop('Must test a fitted model')
  imap <- bsd.imap(ff$map,dim(ff$obs))
  if(is.na(set <- pmatch(set,
                         c('means','sds','correlations','distributions','cutpoints','all'))))
    stop('Test set unrecognized')
  set <- switch(set,
                '1'=c(imap$mu,imap$nu),
                '2'=c(imap$sigma,imap$tau),
                '3'=imap$rho,
                '4'=-c(imap$xi,imap$eta),
                '5'=c(imap$xi,imap$eta),
                '6'=1:(length(ff$estimate)))
  p <- ff$estimate[set]
  varp <- vcov(b)[set,set]
  if (!is.matrix(m)) m <- t(m)
  if (length(p) != dim(m)[2])
    stop('Size of hypothesis matrix must agree with number of parameters')
  df <- dim(m)[1]
  h <- m %*% p
  v <- m %*%varp %*% t(m)
  W <- t(h) %*% solve(v)%*% h
  structure(W,names='W',df=df,class='bsdGOF')
}


#' Compare nested GRT models
#' 
#' Conducts a likelihood-ratio G-test on nested GRT models. Currently only accepts pairs of nested models, not arbitrary sequences.
#' 
#' @param object A fitted GRT model returned by fit.grt
#' @param ... A larger GRT model, with model1 nested inside
#' @export
anova.grt <- function(object, ...){
  model1 <- object;
  model2 <- list(...)[[1]];
  g21 <- GOF(model1,teststat='G2')
  df1 <- attr(g21,'df')
  g22 <- GOF(model2,teststat='G2')
  df2 <- attr(g22,'df')
  DG2 <- round(g21-g22,3)
  ddf <- df1-df2
  p.val <- round(pchisq(DG2,ddf,lower.tail=F),4)
  table <- matrix(c(round(g21,3),round(g22,3),df1,df2,'',
                    DG2,'',ddf,'',p.val),2)
  # Want to get model var names passed in by user
  arg <- deparse(substitute(object))
  dots <- substitute(list(...))[-1]
  modelNames <- c(arg, sapply(dots, deparse))
  dimnames(table) <- list(modelNames,
                          c('G2', 'df', 'DG2', 'df','p-val'))
  as.table(table)
}


# Test parameters for equality
test.parameters <- function(bb,set='means'){
  if (is.null(ff <- bb$fit)) stop('Must work with fitted model')
  imap <- bsd.imap(ff$map,dim(ff$obs))
  if(is.na(set <- pmatch(set, c('means','sds','correlations','cutpoints'))))
    stop('Test set unrecognized')
  c0 <- if (set==2) 1 else 0
  set <- switch(set,
                '1'=c(imap$mu,imap$nu),
                '2'=c(imap$sigma,imap$tau),
                '3'=imap$rho,
                '4'=c(imap$xi,imap$eta))
  p <- ff$estimate[set]
  vp <- vcov(bb)[set,set]
  ix2 <- lower.tri(vp)
  ix1 <- col(vp)[ix2]
  ix2 <- row(vp)[ix2]
  tv <- c(p-c0,p[ix1]-p[ix2])
  dv <- diag(vp)
  se <- sqrt(c(dv,dv[ix1]+dv[ix2]-2*vp[cbind(ix1,ix2)]))
  z <- tv/se
  structure(cbind(tv,se,z,2*pnorm(abs(z),lower.tail=FALSE)),
            dimnames=list(c(paste(names(p),'-',c0),
                            paste(names(p)[ix1],'-',names(p)[ix2])),
                          c('Est','s.e.','z','p')),
            class=c('bsd.test','matrix'))
}

print.bsd.test <- function(mx,digits=3){
  class(mx) <- 'matrix'
  print(round(mx,digits))
}

# Likelihood calculation routines

# Calculate minus the log-likelihood for a bivariate detection model
# pv:      Parameter vector
# x:       Data as a table
# pmap:    A 5 by KK array pmam mapping tables onto parameter vectors
# imap:    List of parameter indices by type in parameter vector
# d.level: Derivitive level to return
#           0: Likelihood only
#           1: Likelihood and first derivative
#           2: Likelihood, first derivative, and expected second derivatives
# diag:    Diagnostic print code: 0, 1, 2, or 3
# fracw:   Value use to space overlapping criteria
bsd.llike <- function (pv,x,pmap,imap,d.level=0,diag=0,fracw=10) {
  tt <- dim(x); II <- tt[1]; JJ <- tt[2]; KK <- tt[3]
  lpv <- length(pv)
  llike <- 0
  if (d.level > 0){
    gradient <- numeric(lpv)
    grx <- numeric(lpv-II-JJ+2)
    if (d.level == 2){
      ExpD2 <- matrix(0,lpv,lpv)
      d2eta <- 1 + (lpv+1)*(0:(II+JJ-1))
      d2xi <- d2eta[imap$xi];   d2xi1 <- (d2xi-1)[-1]
      d2eta <- d2eta[imap$eta]; d2eta1 <- (d2eta-1)[-1]
      nk <- apply(x,3,sum)
    }
  }
  aa <- bsd.map2array(pv,pmap,imap)
  xi <- pv[imap$xi];
  eta <- pv[imap$eta];
  if (diag > 0){
    cat('Call of bsd.llike\n')
    cat('xi ',xi,'\neta',eta,'\n')
    print(aa)
  }
  # Check the validity of parameters
  if (any(c(diff(xi),diff(eta),pv[imap$sigma],pv[imap$tau])<0)
      || any(abs(pv[imap$rho])>1)) {
    warning('Criteria out of order or invalid distribution',
            call.=FALSE)
    llike <- Inf
    if (d.level > 0) attr(llike,'gradient') <- gradient
    return(llike)
  }
  # Loop over tables
  for (k in 1:KK) {
    bf <- bsd.freq(xi,eta,aa[k,],if(d.level==0) 1 else NULL)
    gamma <- if (d.level==0) bf else bf$pi
    llike <- llike - sum(x[,,k]*log(gamma))
    if (diag > 1) {
      cat('Table level ',k,'\nProbabilities\n')
      if(is.list(bf)){
        print(bf$pi)  
      }
      cat('Likelihood contribution ',-sum(x[,,k]*log(gamma)),
          '\nCumulated log-likelihood',llike,'\n')
    }
    # Calculate gradients if required
    if (d.level > 0){
      xg <- x[,,k]/gamma
      t1 <- rbind(bf$d1,0)
      t2 <- rbind(bf$d2,0)
      if (diag > 1) {print(bf)
                     cat('x_ijk/pi_ijk\n'); print(xg)
      }
      vxi <-  D2(t1)[-II,]/aa[k,2]
      veta <- D2(t2)[-JJ,]/aa[k,4]
      if (diag > 2) {
        cat('vxi\n');  print(vxi)
        cat('veta\n'); print(veta)
      }
      # NHS: added 'length != 1' bit; not sure if it's working right...
      if(length(xi)!=1 && length(eta)!=1){
        gd <- c(rowSums(vxi*(xg[-II,]-xg[-1,])),
                rowSums(veta*t(xg[,-JJ]-xg[,-1])), grx)
      }else{
        gd <- c(sum(vxi*(xg[-II,]-xg[-1,])),
                sum(veta*t(xg[,-JJ]-xg[,-1])), grx)
      }
      ipm <- pmap[k,1]
      ips <- pmap[k,2]
      ipn <- pmap[k,3]
      ipt <- pmap[k,4]
      ipr <- pmap[k,5]
      if (ipm > 0){
        ixa <- imap$mu[ipm]
        vmu <- -D12(t1)/aa[k,2]
        if (diag > 2) {cat('vmu\n'); print(vmu)}
        gd[ixa] <- sum(xg*vmu)
      }
      if (ips > 0) {
        ixb <- imap$sigma[ips]
        vsigma <- -D12(((c(xi,0)-aa[k,1])/aa[k,2]^2)*t1)
        if (diag > 2) {cat('vsigma\n'); print(vsigma)}
        gd[ixb] <- sum(xg*vsigma)
      }
      if (ipn > 0) {
        ixk <- imap$nu[ipn]
        vnu <- t(-D12(t2))/aa[k,4]
        if (diag > 2) {cat('vnu\n'); print(vnu)}
        gd[ixk] <- sum(xg*vnu)
      }
      if (ipt > 0) {
        ixl <- imap$tau[ipt]
        vtau <- t(-D12(((c(eta,0)-aa[k,3])/aa[k,4]^2)*t2))
        if (diag > 2) {cat('vtau\n'); print(vtau)}
        gd[ixl] <- sum(xg*vtau)
      }
      if (ipr > 0) { 
        ixr <- imap$rho[ipr]
        vrho <- D12(rbind(cbind(bf$phi,0),0))
        if (diag > 2) {cat('vrho\n'); print(vrho)}
        gd[ixr] <- sum(xg*vrho)
      }
      gradient <- gradient + gd
      if (diag > 0) {
        cat('First derivative contribution\n')
        names(gd) <- names(pv)
        print(round(gd,4))
      }
      # Calculate expected second derivatives if requested
      # if(length(xi)!=1) and if(length(eta)!=1) clauses added 1.25.14 -NHS
      if (d.level == 2){
        oog <- 1/gamma
        if(length(eta)!=1){
          veta <- t(veta)          
        }
        hs <- matrix(0,lpv,lpv)
        
        if(length(xi)!=1){
          hs[d2xi] <- rowSums(vxi^2*(oog[-II,]+oog[-1,]))
          tt <- vxi[-(II-1),]*vxi[-1,]*oog[2:(II-1),]
        }else{
          hs[d2xi] <- sum(vxi^2*(oog[-II,]+oog[-1,]))
          tt <- vxi[-(II-1)]*vxi[-1]*oog[2:(II-1),]
        }
        hs[d2xi1] <- -if (is.matrix(tt)) rowSums(tt) else sum(tt)
        if(length(eta)!=1){
          hs[d2eta] <- colSums(veta^2*(oog[,-JJ]+oog[,-1]))
          tt <- veta[,-(JJ-1)]*veta[,-1]*oog[,2:(JJ-1)]          
        }else{
          hs[d2eta] <- sum(veta^2*(oog[,-JJ]+oog[,-1]))
          tt <- veta[-(JJ-1)]*veta[-1]*oog[,2:(JJ-1)]
        }
        hs[d2eta1] <- -if (is.matrix(tt)) colSums(tt) else sum(tt)
        
        if(length(xi)!=1 && length(eta)!=1){
          hs[imap$xi,imap$eta] <- 
            vxi[,-JJ]*(veta[-II,]*oog[-II,-JJ] - veta[-1,]*oog[-1,-JJ]) -
            vxi[,-1]*(veta[-II,]*oog[-II,-1] - veta[-1,]*oog[-1,-1])          
        }else{
          hs[imap$xi,imap$eta] <- 
            vxi[-JJ]*(veta[-II]*oog[-II,-JJ] - veta[-1]*oog[-1,-JJ]) -
            vxi[-1]*(veta[-II]*oog[-II,-1] - veta[-1]*oog[-1,-1])
        }
        if (ipm > 0){
          tt <- vmu*oog
          hs[ixa,ixa] <- sum(vmu*tt)
          if(length(xi)!=1){
            hs[imap$xi,ixa] <- rowSums(vxi*(tt[-II,] - tt[-1,]))
          }else{
            hs[imap$xi,ixa] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixa] <- colSums(veta*(tt[,-JJ] - tt[,-1]))            
          }else{
            hs[imap$eta,ixa] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ips > 0) hs[ixa,ixb] <- sum(tt*vsigma)
          if (ipn > 0) hs[ixa,ixk] <- sum(tt*vnu)
          if (ipt > 0) hs[ixa,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixa,ixr] <- sum(tt*vrho)
        }
        if (ips > 0){
          tt <- vsigma*oog
          hs[ixb,ixb] <- sum(tt*vsigma)
          hs[imap$xi,ixb] <- rowSums(vxi*(tt[-II,] - tt[-1,]))
          hs[imap$eta,ixb] <- colSums(veta*(tt[,-JJ] - tt[,-1]))
          if (ipn > 0) hs[ixb,ixk] <- sum(tt*vnu)
          if (ipt > 0) hs[ixb,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixb,ixr] <- sum(tt*vrho)
        }
        if (ipn > 0){
          tt <- vnu*oog
          hs[ixk,ixk] <- sum(tt*vnu)
          if(length(xi)!=1){
            hs[imap$xi,ixk] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixk] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixk] <- colSums(veta*(tt[,-JJ] - tt[,-1])) 
          }else{
            hs[imap$eta,ixk] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ipt > 0) hs[ixk,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixk,ixr] <- sum(tt*vrho)
        }
        if (ipt > 0){
          tt <- vtau*oog
          hs[ixl,ixl] <- sum(tt*vtau)
          if(length(xi)!=1){
            hs[imap$xi,ixl] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixl] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixl] <- colSums(veta*(tt[,-JJ] - tt[,-1])) 
          }else{
            hs[imap$eta,ixl] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ipr > 0) hs[ixl,ixr] <- sum(tt*vrho)
        }
        if (ipr > 0){
          tt <- vrho*oog
          hs[ixr,ixr] <- sum(tt*vrho)
          if(length(xi)!=1){
            hs[imap$xi,ixr] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixr] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixr] <- colSums(veta*(tt[,-JJ] - tt[,-1]))            
          }else{
            hs[imap$eta,ixr] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
        }
        ExpD2 <- ExpD2 - nk[k]*hs
        if (diag > 2) {
          cat('Expected second derivative contribution\n')
          print(round(-nk[k]*hs,3))
        }
      }
    }
  }
  if (diag > 0) cat('Log-likelihood:',llike,'\n')
  if (d.level > 0){
    attr(llike,'gradient') <- -gradient
    if (d.level == 2) {
      hs <- ExpD2
      diag(hs) <- 0
      attr(llike,'ExpD2') <- ExpD2+t(hs)
    }
  }
  llike
}


# Differencing of first and second subscript of matrix
D1 <- function(x){
  dx <- dim(x)[1]
  rbind(x[1,],x[-1,]-x[-dx,])
}

D2 <- function(x){
  dx <- dim(x)[2]
  cbind(x[,1],x[,-1]-x[,-dx])
}

D12 <- function(x){
  r <- dim(x)[1]; c <- dim(x)[2]
  x <- rbind(x[1,],x[-1,]-x[-r,])
  cbind(x[,1],x[,-1]-x[,-c])
}

bsd.like <- function(p,...) bsd.llike(p,d.level=0,...)
bsd.grad <- function(p,...) attr(bsd.llike(p,d.level=1,...),'gradient')

# Calculate the cell frequencies for a bivariate SDT model.
# xi and eta: Row and column criteria
# m:          A vector of the distributional parameters
#               (mu_r, sigma_r, mu_c, sigma_c, rho).
# n:          Sample size or NULL
# When a sample size n is given, the function returns expected frequencies;
# when it is NULL, the function returns a list containing the probabilities pi,
# the densities phi, and the row and column derivative terms (the latter
# as its transpose).
bsd.freq <- function (xi,eta,m,n=NULL) {
  #require(mvtnorm)
  fracw <- 10
  nrow <- length(xi) +1
  ncol <- length(eta)+ 1
  Xi  <- c(-Inf, (xi-m[1])/m[2], Inf)
  Eta <- c(-Inf, (eta-m[3])/m[4], Inf)
  rho <- m[5]
  pii <- matrix(nrow=nrow, ncol=ncol)
  cx <- matrix(c(1,rho,rho,1),2)
  for (i in 1:nrow) for (j in 1:ncol) 
    pii[i,j] <- sadmvn(lower = c(Xi[i],Eta[j]), upper = c(Xi[i+1],Eta[j+1]), mean = rep(0,2), varcov=cx)
  if (is.null(n)){
    Xis <- Xi[2:nrow]
    Etas <- Eta[2:ncol]
    phi <- matrix(0, nrow=nrow-1, ncol=ncol-1)
    for (i in 1:(nrow-1)) for (j in 1:(ncol-1))
      phi[i,j] <- dmnorm(c(Xis[i],Etas[j]),varcov=cx)
    list(pi = pii, phi = phi,
         d1 = dnorm(Xis)*pnorm(outer(-rho*Xis,c(Etas,Inf),'+')/sqrt(1-rho^2)), 
         d2 = dnorm(Etas)*pnorm(outer(-rho*Etas,c(Xis,Inf),'+')/sqrt(1-rho^2)))
  }
  else pii*n
}


# Checks that a parameter map has the correct form
bsd.valid.map <- function(map,K){
  dm <- dim(map)
  if (!is.matrix(map)) stop('Map must be a matrix')
  if (dm[1] != K) stop('Map must have same number of rows as conditions')
  if (dm[2] != 5) stop('Map must have 5 columns')
  for (i in 1:4){
    u <- unique(map[,i])
    if (min(u) != 0) stop(paste('Values in map column',i,'must start at 0'))
    if (max(u) != length(u)-1) 
      stop(paste('Map column',i,'must be dense series'))
  }
  u <- unique(map[,5])
  nu <- min(u); xu <- max(u)
  if (((nu==0) && (xu!=length(u)-1)) || ((nu==1) & (xu!=length(u))))
    stop('Map column 5 must start at 0 or 1 and be dense series')
  TRUE
}


# Changes notation from cutpoints to differences of cutpoints and back
#bsd.todiff <- function(p,ixxi,ixeta){
#  c(p[ixxi[1]],diff(p[ixxi]),p[ixeta[1]],diff(p[ixeta]),
#     p[(max(ixeta)+1):length(p)])
#  }
#bsd.tocuts <- function(p,ixxi,ixeta){
#  c(cumsum(p[ixxi]),cumsum(p[ixeta]),p[(max(ixeta)+1):length(p)])
#  }


# Constructs the parameter indices
# pmap is parameter map and dimx is dimension of data
bsd.imap <- function(pmap,dimx){
  II <- dimx[1]; JJ <- dimx[2]
  # Handle fact that 2x2 case doesn't use cutpoints (6.30.14 -- RDH)
  if(II > 2) {
    ixxi = 1:(II-1);
    ib <- II+JJ-1;
    ixeta = II:(ib-1);
  } else {
    ixxi = NULL;  
    ib <- 1;
    ixeta <- NULL;
  }
  inp <- max(pmap[,1]); ixmu  <- ib:(ib+inp-1);     ib <- ib+inp
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,2])>0){
    inp <- max(pmap[,2]); ixsigma <- ib:(ib+inp-1); ib <- ib+inp    
  }else{
    ixsigma <- NULL#inp <- 1; ixsigma <- ib; ib <- ib+inp
  }
  inp <- max(pmap[,3]); ixnu  <- ib:(ib+inp-1);     ib <- ib+inp
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,4])>0){
    inp <- max(pmap[,4]); ixtau <- ib:(ib+inp-1); ib <- ib+inp    
  }else{
    ixtau <- NULL#inp <- 1; ixtau <- ib; ib <- ib+inp
  }
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,5])>0){
    inp <- max(pmap[,5]); ixrho <- ib:(ib+inp-1)    
  }else{
    ixrho <- NULL
  }
  list(xi=ixxi,eta=ixeta,mu=ixmu,sigma=ixsigma,nu=ixnu,
       tau=ixtau,rho=ixrho)
}


# Takes a parameter vector 'p' and a map of parameters are returns a
# table of parameter values.
# 'm0' and 'x0' are mean and variance values for parameter with map 0.
bsd.map2array <- function(p,pmap,imap,m0=0,s0=1){
  KK <- dim(pmap)[1]
  aa <- matrix(m0,KK,5)
  rownames(aa) <- rownames(pmap)
  colnames(aa) <- colnames(pmap)
  for (k in 1:KK) {
    if (pmap[k,1] != 0) aa[k,1] <- p[imap$mu[pmap[k,1]]]
    aa[k,2] <- if (pmap[k,2] != 0) p[imap$sigma[pmap[k,2]]] else s0
    if (pmap[k,3] != 0) aa[k,3] <- p[imap$nu[pmap[k,3]]]
    aa[k,4] <- if (pmap[k,4] != 0) p[imap$tau[pmap[k,4]]] else s0
    if (pmap[k,5] != 0) aa[k,5] <- p[imap$rho[pmap[k,5]]]
  }
  aa}


# Construct an initial vector
# The value delta is added to all frequencies to avoid problems with zeros
bsd.initial <- function(xx,pmap,delta=0.5){
  pnames <- c('mu','sigma','nu','tau','rho')
  dxx <- dim(xx)
  II <- dxx[1]; JJ <- dxx[2];  KK <- dxx[3];
  ixx <- 1:(II-1)
  ixy <- 1:(JJ-1)
  xx <- xx + delta
  ptx <- prop.table(margin.table(xx,c(3,1)),1)
  pty <- prop.table(margin.table(xx,c(3,2)),1)
  xi <- rep(0,II); eta <- rep(0,JJ)
  ni <- nj <- 0
  for (k in 1:KK){
    ptx[k,] <- qnorm(cumsum(ptx[k,]))
    if (pmap[k,1]+pmap[k,2]==0) {xi <- xi+ptx[k,]; ni <- ni+1}
    pty[k,] <- qnorm(cumsum(pty[k,]))
    if (pmap[k,3]+pmap[k,4]==0) {eta <- eta+pty[k,]; nj <- nj+1}
  }
  # NHS 1.20.14
  # added 'as.matrix' terms to maintain KK rows in ptx, pty
  # code still doesn't work with 2x2 data
  # si matrix (below) ends up with NAs in cols 2 and 4
  ptx <- as.matrix(ptx[,ixx]); pty <- as.matrix(pty[,ixy])
  xi <- (xi/ni)[ixx]; eta <- (eta/nj)[ixy]
  pv <- c(xi,eta)
  nv <- c(paste('xi',ixx,sep=''),paste('eta',ixy,sep=''))
  np <- apply(pmap,2,max)
  si <- matrix(0,KK,4)
  for (k in 1:KK){
    si[k,1:2] <- coef(lm(ptx[k,]~xi)) 
    si[k,3:4] <- coef(lm(pty[k,]~eta))
  }
  if(II>2){
    si[,1] <- -si[,1]/si[,2]
    si[,2] <- 1/si[,2]    
  }else{
    si[,1] <- -si[,1]
    si[,2] <- 1
  }
  if(JJ>2){
    si[,3] <- -si[,3]/si[,4]
    si[,4] <- 1/si[,4]    
  }else{
    si[,3] <- -si[,3]
    si[,4] <- 1
  }
  for (i in 1:4) if (np[i] > 0) {
    pv <- c(pv,tapply(si[,i],pmap[,i],mean)[-1])
    nv <- c(nv, paste(pnames[i],1:np[i],sep=''))
  }
  if (np[5]>0){
    r <- rep(0,KK)
    for (k in 1:KK) r[k] <- polychor(xx[,,k])
    rr <- tapply(r,pmap[,5],mean)
    if (min(pmap[,5]) == 0) rr <- rr[-1]
    pv <- c(pv,rr)
    nv <- c(nv, paste('rho',1:np[5],sep=''))
  }
  names(pv) <- nv
  pv
}


# Some test functions

test.bsd.llike <- function(pv,xx,pmap,n=1,d.level=0,diag=d.level+1){
  dxx <- dim(xx)
  imap <- bsd.imap(pmap,dxx)
  cat('Parameter mapping vector\n');   print(pmap)
  cat('Parameters by groups\n')
  print(bsd.map2array(pv,pmap,imap))
  bsd.llike(pv,xx,pmap,imap,d.level=d.level,diag=diag)
}

test.bsd.freq <- function(xi=c(-.6,.15,.65), eta=c(-.5,.25,.75),
                          m=c(0,1,.2,1.2,-.3),n=NULL) {
  print('Calling bsd.freq')
  bsd.freq(xi,eta,m,n)
}

# Take trace of matrix
tr <- function (m) {
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
    stop("m must be a square matrix")
  return(sum(diag(m)))
}

estimate_prob_and_var <- function(xpar,ypar,rpar,ps_old){
  prob <- matrix(data=0, nrow = 4, ncol = 4)
  v <- array(0, dim = c(4,4,3)); 
  for (i in 1:4){
    # Bookkeeping
    x_i <- if(length(xpar)==2) ceiling(i/2) else i;
    alpha = ps_old[xpar[x_i]];
    y_i <- if(length(ypar)==2) ((i-1) %% 2) + 1 else i;
    kappa = ps_old[ypar[y_i]];
    if (is.null(rpar)) {
      rho = 0;
    } else if (length(rpar) == 1) {
      rho = ps_old[rpar];
    } else {
      rho = ps_old[rpar[i]];
    }
    prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
    v[i,,] = vcalc(alpha,kappa,rho);
  }
  return(list(prob=prob,v=v));
}

make_parameter_mat <- function(xpar, ypar, rpar, ps_new){
  if (length(xpar) == 2) {
    mu = c(ps_new[1], ps_new[1], ps_new[2], ps_new[2]); offset = 2;
  } else {
    mu = c(ps_new[1],ps_new[2],ps_new[3],ps_new[4]); offset = 4;
  }
  if (length(ypar) == 2) {
    nu = c(ps_new[offset+1], ps_new[offset+2], ps_new[offset+1], ps_new[offset+2]);
    offset = offset + 2;
  } else {
    nu = c(ps_new[offset+1],ps_new[offset+2], ps_new[offset+3], ps_new[offset+4]);
    offset = offset + 4;
  }
  if (is.null(rpar)) {
    rho = rep(0,4);
  } else if (length(rpar) == 1) {
    rho = rep(ps_new[offset+1],4);
  } else {
    rho = c(ps_new[offset+1], ps_new[offset+2], ps_new[offset+3], ps_new[offset+4]);
  }
  sigma = rep(1.0,4);
  tau = rep(1.0,4);
  return(cbind(mu,sigma,nu,tau,rho));
}
# initialize various scalars and arrays 

# parameters in order:
# mu_x_** mu_y_** rho_**
# where ** = aa, ab, ba, bb
initial_point <- function(prob, PS_x, PS_y, PI) {
  nx =0; ny = 0; nr = 0;xpar = NULL;ypar=NULL;rpar=NULL;
  # Figure out how many params we need
  if (PS_x) {
    xpar=1:2; nx=2; 
    rows=matrix(data=rbind(c(1,1,0,0),c(0,0,1,1)),nrow=2,ncol = 4);
  } else {
    xpar=1:4; nx=4;
    rows=matrix(data=diag(4),nrow=4,ncol=4);
  }
  if (PS_y) {
    ypar = nx + 1:2; ny=2;
    rows=rbind(rows,rbind(c(1,0,1,0),c(0,1,0,1)));
  } else {
    ypar = nx + 1:4; ny=4;
    rows=rbind(rows,diag(4));
  }
  if (PI == 'same_rho') {
    rpar = nx+ny+1; nr = 1;
    rows=rbind(rows,c(1,1,1,1));
  } 
  if (PI == 'none') {
    rpar = nx+ny+1:4; nr=4;
    rows=rbind(rows, diag(4));
  }  
  npar = nx + ny + nr;
  rows = matrix(data=as.logical(rows),ncol=4,nrow=npar);
  param_estimate = matrix(data = 0, nrow= npar, ncol = 1); # For param estimates
  # initial estimates: y means
  if (PS_x) {
    for (i in c(1,3)) { 
      param_estimate[xpar[ceiling(i/2)]] = -qnorm(.5*(prob[i,1]   + prob[i,2]) 
                                                  + .5*(prob[i+1,1] + prob[i+1,2]));}
  } else {
    for (i in 1:4) {
      param_estimate[xpar[i]] = -qnorm(prob[i,1] + prob[i,2]);}
  }
  # initial estimates: y means
  if (PS_y) {
    for (i in c(1,2)) {
      param_estimate[ypar[i]] = -qnorm(.5*(prob[i,1]   + prob[i,3])
                                       + .5*(prob[i+2,1] + prob[i+2,3])); }
  } else {
    for (i in 1:4) {
      param_estimate[ypar[i]] = -qnorm(prob[i,1] + prob[i,3]);}
  }  
  # initial estimates: correlation  
  if (PI=='same_rho') {
    r = cos(pi/(1+sqrt((.25*sum(prob[,4])*.25*sum(prob[,1]))
                       /(.25*sum(prob[,3])*.25*sum(prob[,2])))));
    if (r <= -1) { r = -.95; }
    else if (r >= 1) { r = .95; }
    param_estimate[rpar] = r;
  } else if (PI == 'none') {
    for (i in 1:4) {
      r = cos(pi/(1+sqrt((prob[i,4]*prob[i,1])/(prob[i,3]*prob[i,2]))));
      if (r <= -1) { r = -.95; }
      else if (r >= 1) { r = .95; }
      param_estimate[rpar[i]] = r;
    }
  }
  return(list(xpar=xpar, ypar=ypar, rpar=rpar, ps_old=param_estimate, rows=rows))
}  

create_two_by_two_mod <- function(PS_x, PS_y ,PI) {
  mod <- matrix(data = 0, nrow = 1, ncol = 7);
  if (PS_x) mod[1] = 1;
  if (PS_y) mod[2] = 1;
  if (PI == 'all') {
    mod[3:6] = rep(1,times=4);
  } else { 
    if (PI == 'same_rho') mod[7] = 1;
  }
  return(mod);
}

defaultNames <- c("a_1b_1", "a_1b_2", "a_2b_1", "a_2b_2")

# In 2x2 case, it's typical to use a 4x4 frequency matrix w/ each row being a stim
# and each col being the freqency of responding "aa", "ab", "ba", "bb", respectively, 
# to that stim. Wickens' code for nxn case requires data in xtabs format.
freq2xtabs <- function(freq) {
  xdim = dim(freq)[1]; 
  ydim = dim(freq)[2];
  d = as.data.frame(matrix(rep(x=0,times=xdim*ydim*4), nrow = xdim*ydim, ncol = 4));
  stimuli = if(length(rownames(freq)) > 0) rownames(freq) else defaultNames; 
  names(d) <- c("Stim", "L1", "L2", "x");
  for (i in 1:4) {
    for (j in 1:4) {
      d[4*(i-1) + j,] = c(stimuli[i], floor((j+1) / 2), ((j-1) %% 2) + 1, freq[i,j]);
    }
  }
  d$Stim <- ordered(d$Stim,levels=stimuli);
  d$L1  <- ordered(d$L1);
  d$L2 <- ordered(d$L2);
  d$x <- as.numeric(d$x);
  return(xtabs(x~L1+L2+Stim, d));
}

# calculate predicted response probabilities for the given stimulus
# b/c we assume decisional sep, each response is a quadrant of Cartesian plane
prcalc <- function(mean, cov) {
  pr <- matrix(data = 0, nrow = 1, ncol = 4);
  pr[1,1] = sadmvn(lower = c(-Inf, -Inf), upper = c(0, 0), mean, cov);
  pr[1,2] = sadmvn(lower = c(-Inf, 0), upper = c(0, +Inf), mean, cov);
  pr[1,3] = sadmvn(lower = c(0, -Inf), upper = c(+Inf, 0), mean, cov);
  pr[1,4] = sadmvn(lower = c(0, 0), upper = c(+Inf, +Inf), mean, cov);
  return(pr);
}

# calculate v-matrix elements
vcalc <- function(ap,kp,rh) {        
  ve <- matrix(data = 0, ncol = 3, nrow = 4)
  d_mx <- matrix(data = 0, ncol = 2, nrow = 2)
  d_my <- matrix(data = 0, ncol = 2, nrow = 2)
  
  d_mx_arg = (rh*ap-kp)/sqrt(1-rh^2);
  d_mx[1,1] = -dnorm(-ap)*pnorm( d_mx_arg );
  d_mx[1,2] = -dnorm(-ap);
  d_mx[2,1] = 0;
  d_mx[2,2] = 0;
  
  ve[1,1] =  d_mx[1,1];
  ve[2,1] =  d_mx[1,2] - d_mx[1,1];
  ve[3,1] =  d_mx[2,1] - d_mx[1,1];
  ve[4,1] =  d_mx[2,2] - d_mx[2,1] - d_mx[1,2] + d_mx[1,1];
  
  d_my_arg = (rh*kp-ap)/sqrt(1-rh^2);
  d_my[1,1] = -dnorm(-kp)*pnorm( d_my_arg );
  d_my[1,2] = 0;
  d_my[2,1] = -dnorm(-kp);
  d_my[2,2] = 0;
  
  ve[1,2] = d_my[1,1];
  ve[2,2] = d_my[1,2] - d_my[1,1];
  ve[3,2] = d_my[2,1] - d_my[1,1];
  ve[4,2] = d_my[2,2] - d_my[1,2] - d_my[2,1] + d_my[1,1];
  
  S_aa = matrix(c(1, rh, rh, 1), ncol = 2, nrow = 2);
  d_rho = dmnorm(c(-ap, -kp),mean = c(0, 0), varcov = S_aa);
  
  ve[1,3] = d_rho;
  ve[2,3] = -d_rho;
  ve[3,3] = -d_rho;
  ve[4,3] = d_rho;
  
  return(ve);
}

get_fit_params <- function(grt_obj) {
  d = distribution.parameters(bb=grt_obj);
  return(list(aa=d[1,c(1,3,5)], ab = d[2,c(1,3,5)], 
              ba = d[3, c(1,3,5)], bb = d[4,c(1,3,5)]));
}
checkConfusionMatrix <- function(x) {
  dimx <- dim(x)[1]
  if( dimx != dim(x)[2]){ 
    cat("Confusion matrix must have an equal number of rows and columns!\n")
    return(FALSE)
  }
  
  if(max(x)<=1 & min(x) >=0) {
    if(all( apply(x, 1, sum) == rep(1,dimx))) {
      return(TRUE)
    } else {
      cat("The rows of confusion probability matrix must sum to one!\n")
      return(FALSE)
    }
  } else {
    return(TRUE)}
}
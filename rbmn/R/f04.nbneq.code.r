
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn2gema <- function(nbn)
#TITLE computes a /gema/ from  a /nbn/
#DESCRIPTION from a /nbn/ object defining a normal
# Bayesian network, computes the vector \code{mu} and
# the matrix \code{li} such that if the vector \code{E}
# is a vector of i.i.d. centred and standardized normal, then
# \code{mu + li %*% E} has the same distribution as the
# input /nbn/.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object for which the generating matrices.>>
#[INPUTS]
#VALUE
# a list with the two following components:
#     \code{mu} and \code{li}.
#EXAMPLE
# identical(nbn2gema(rbmn0nbn.02),rbmn0gema.02);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Doesn't work for rbmn0nbn.04 !!!
#AUTHOR J.-B. Denis
#CREATED 12_01_14
#REVISED 12_01_17
#--------------------------------------------
{
  # checking
  # < to be done >
  # number of nodes
  nn <- length(nbn);
  nam <- names(nbn);
  # initializing the matrices
  mu <- rep(0,nn);
  li <- matrix(0,nn,nn);
  names(mu) <- nam;
  dimnames(li) <- list(nam,NULL);
  # checking or putting a topological order
  if (!order4nbn(nbn,r.bc(nn))) {
    tord <- order4nbn(nbn);
    nbn <- nbn[tord];
  }
  # processing
  for (ii in r.bc(nn)) {
    iino <- nam[ii];
    pare <- nbn[[ii]]$parents;
    muco <- nbn[[ii]]$mu;
    sico <- nbn[[ii]]$sigma;
    if (length(pare)==0) {
      mu[iino] <- muco;
      #
      li[iino,] <- sico * ((1:nn)==ii);
    } else {
      regr <- nbn[[ii]]$regcoef;
      #
      mu[iino] <- muco + sum(regr * mu[pare]);
      #
      li[iino,] <- sico * ((1:nn)==ii) +
                     matrix(regr,nrow=1) %*% li[pare,];
    }
  }
  # returning
  list(mu=mu,li=li);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
gema2nbn <- function(gema)
#TITLE computes a /nbn/ from  a /gema/
#DESCRIPTION from a /gema/ object defining a normal
# Bayesian network, computes more standard
# /nbn/ where each node is defined from its parents.
#DETAILS
# using general formulae rather a sequential algorithm
# as done in the original \code{gema2nbn} implementation.
#KEYWORDS 
#INPUTS
#{gema} <<Initial \code{gema} object.>>
#[INPUTS]
#VALUE
# the corresponding /nbn/.
#EXAMPLE
# print8nbn(gema2nbn(rbmn0gema.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 12_01_16
#REVISED 12_01_17
#--------------------------------------------
{
  # zero limit
  eps <- 10^-10;
  # number of nodes
  nn <- nrow(gema$li);
  # checking
  if (nn < 1) {
    stop(paste(nn,"is not a sufficient number of variables"));
  }
  if (sum(gema$li[outer(r.bc(nn),r.bc(nn),"<")]^2) > eps) {
    cat("For eps =",eps,"the 'gema$li' matrix is not triangular\n");
    stop("Bad argument 'gema'");
  }
  dd <- svd(gema$li)$d;
  if (abs(dd[nn]) < eps) {
    stop("With respect to 'eps', the 'gema$li' matrix is singular");
  }
  # the naming          
  nam <- dimnames(gema$li)[[1]];
  dimnames(gema$li) <- list(nam,nam);
  # getting the conditional standard deviations
  sigma <- diag(gema$li);
  # getting the conditional expectations
  s1 <- diag(diag(gema$li),nn) %*% solve(gema$li);
  lambda <- s1 %*% matrix(gema$mu,nn);
  # getting the linking coefficients
  rho <- diag(1,nrow=nn) - s1;
  rho <- rho * (rho > eps);
  # initializing the /nbn/
  res <- vector("list",nn);
  names(res) <- nam;
  # filling each node
  for (ii in r.bc(nn)) {
    #pare <- which(abs(gema$li[ii,]) > eps);
    pare <- which(abs(rho[ii,]) > eps);
    pama <- nam[setdiff(pare,ii)];
    coef <- rho[ii,];
    names(coef) <- nam;
    if (length(pama)>0) {
      res[[ii]]$parents <- pama;
      res[[ii]]$regcoef <- coef[pama];
    }
    res[[ii]]$sigma <- sigma[ii];
    res[[ii]]$mu <- lambda[ii];
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
gema2mn <- function(gema)
#TITLE computes a /mn/ from a /gema/
#DESCRIPTION from a /gema/ object defining a normal
# Bayesian network, computes the expectation and
# variance matrix.
#DETAILS
#KEYWORDS 
#INPUTS
#{gema} <<Initial \code{gema} object.>>
#[INPUTS]
#VALUE
# a list with the following components:
#     \code{mu} and \code{gamma}.
#EXAMPLE
# print8mn(gema2mn(rbmn0gema.04));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#AUTHOR J.-B. Denis
#CREATED 12_01_14
#REVISED 12_01_16
#--------------------------------------------
{
  # checking
  # < to be done >
  # getting the four other matrices
  gamma <- gema$li %*% t(gema$li);
  # returning
  list(mu=gema$mu,gamma=gamma);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
mn2gema <- function(mn)
#TITLE computes a /gema/ from a /mn/
#DESCRIPTION proposes generating matrices of a Bayesian network
# from a /mn/ object defining a multinormal
# distribution by expectation and variance,
# under the assumption that the nodes are in topological order.
#DETAILS
#KEYWORDS 
#INPUTS
#{mn} <<Initial \code{mn} object.>>
#[INPUTS]
#VALUE
# a list with the /gema/ components
#     \code{$mu} and \code{$li}.
#EXAMPLE
# print8gema(mn2gema(rbmn0mn.04));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#AUTHOR J.-B. Denis
#CREATED 12_04_27
#REVISED 12_04_27
#--------------------------------------------
{
  # checking
  # < to be done >
  # getting the triangular matrix
  li <- chol(mn$gamma);
  # returning
  list(mu=mn$mu,li=t(li));
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn2nbn <- function(nbn,norder) 
#TITLE computes the /nbn/ changing its topological order
#DESCRIPTION returns the proposed /nbn/ with a
# new topological order without modifying 
# the joint distribution of all variables.\cr
# This allows to directly find regression
# formulae within the Gaussian Bayesian networks.
#DETAILS
# BE aware that for the moment, no check is made about the topological
# order and if it is not, the result is FALSE!
#KEYWORDS 
#INPUTS
#{nbn}<< The /nbn/ to transform.>>
#{norder}<< The topological order to follow.
#   It can be indicated by names or numbers.
# When not all nodes are included, the resulting
# /nbn/ is restricted to these nodes after marginalization.>>
#[INPUTS]
#VALUE
# The resulting /nbn/.
#EXAMPLE
# print8mn(nbn2mn(rbmn0nbn.01,algo=1));
# print8mn(nbn2mn(rbmn0nbn.01,algo=2));
# print8mn(nbn2mn(rbmn0nbn.01,algo=3));
# print8mn(nbn2mn(nbn2nbn(rbmn0nbn.02,c(1,2,4,5,3))));
# print8mn(nbn2mn(nbn2nbn(rbmn0nbn.02,c(4,1,2,3,5))));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE check the topological order
#AUTHOR J.-B. Denis
#CREATED 13_01_28
#REVISED 13_01_28
#--------------------------------------------
{
  # constants
  nn <- length(nbn);
  # checking
  if (length(norder)!=length(unique(norder))) {
    stop("duplicate nodes in 'norder'!");
  }
  if (!is.character(norder)) {
    if (!all(norder <= nn)) {
      stop("'norder' not consistent with 'nbn' [1]");
    }
    norder <- names(nbn)[norder];
  } else {
    if (!setequal(union(norder,names(nbn)),names(nbn))) {
      stop("'norder' not consistent with 'nbn' [2]");
    }
  }
  # going to the joint distribution
  mn1 <- gema2mn(nbn2gema(nbn));
  # ordering/restricting it
  mn2 <- vector("list",0);
  mn2$mu <- mn1$mu[norder];
  mn2$gamma <- mn1$gamma[norder,norder];
  # transforming it to a nbn
  res <- gema2nbn(mn2gema(mn2));
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rmatrix4nbn <- function(nbn,stdev=TRUE) 
#TITLE regression matrix of a /nbn/
#DESCRIPTION returns a dimnamed matrix
# indicating with \code{rho} an arc from row to column nodes
# (0 everywhere else) where \code{rho} is the regression
# coefficient. Also conditional standard deviations can be
# introduced as diagonal elements but \code{mu} coefficient are
# lost... It is advisable to normalize the /nbn/ first.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn}<< The initial \code{nbn} object.>>
#[INPUTS]
#{stdev} <<Indicates if the standard deviations must placed
#          in the diagonal positions.>>
#VALUE
# A dimnamed matrix
#EXAMPLE
# rmatrix4nbn(rbmn0nbn.02);
# (rmatrix4nbn(rbmn0nbn.02,FALSE)>0)*1;
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_22
#REVISED 13_04_22
#--------------------------------------------
{
  # checking
  # To be done
  # getting the parentship matrix
  nbno <- length(nbn);
  nbna <- names(nbn);
  res <- matrix(0,nbno,nbno);
  dimnames(res) <- list(from=nbna,to=nbna);
  for (nn in r.bf(nbn)) {
    if (length(nbn[[nn]]$parents) > 0) {
      res[match(nbn[[nn]]$parents,nbna),nn] <- nbn[[nn]]$regcoef;
    }
    if (stdev) {
      res[nn,nn] <- nbn[[nn]]$sigma;
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn4rmatrix <- function(rmatrix) 
#TITLE a /nbn/ from a regression matrix 
#DESCRIPTION 
# reverse of \code{rmatrix4nbn} but the standard 
# deviations must be included.
#DETAILS
# \code{mu}s are put to nought
#KEYWORDS 
#INPUTS
#{rmatrix}<< The regression coefficient matrix with the
#            standard deviations in the diagonal.>>
#[INPUTS]
#VALUE
# A /nbn/ object
#EXAMPLE
# print8nbn(nbn4rmatrix(rmatrix4nbn(rbmn0nbn.02)));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_22
#REVISED 13_04_24
#--------------------------------------------
{
  # checking
  # To be done
  # getting the parentship matrix
  nbno <- nrow(rmatrix);
  nbna <- dimnames(rmatrix)[[1]];
  res <- vector("list",nbno);
  names(res) <- nbna;
  for (nn in r.bc(nbno)) {
    res[[nn]]$mu <- 0;
    res[[nn]]$sigma <- rmatrix[nn,nn];
    rmatrix[nn,nn] <- 0;
    ou <- which(rmatrix[,nn]!=0);
    if (length(ou) > 0) {
      res[[nn]]$parents <- nbna[ou];
      res[[nn]]$regcoef <- rmatrix[ou,nn];
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn2rr <- function(nbn)
#TITLE computes standard matrices from  a /nbn/
#DESCRIPTION from a /nbn/ object defining a normal
# Bayesian network, returns a list comprising (i) \code{mm}
# the vector of the mean of the different nodes when
# the parents are nought, (ii) \code{ss} the vector of the conditional
# standard deviations and (iii) \code{rr} the matrix of the 
# regression coefficients of the direct parents (\code{rr[i,j]} contains
# the regression coefficient of the node \code{j} for its parents \code{i}
# or zero when \code{i} is not a parent of \code{j}.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object.>>
#[INPUTS]
#VALUE
# the resulting list with the three components:
#     \code{mm}, \code{ss} and \code{rr}.
#EXAMPLE
# nbn2rr(rbmn0nbn.01);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_04_19
#REVISED 13_04_19
#--------------------------------------------
{
  # checking
  # < to be done >
  # number of nodes
  nn <- length(nbn);
  nam <- names(nbn);
  # initializing the matrices
  mm <- rep(0,nn);
  rr <- matrix(0,nn,nn);
  names(mm) <- nam;
  ss <- mm;
  dimnames(rr) <- list(nam,nam);
  # filling the three components
  for (ii in r.bc(nn)) {
    iino <- nam[ii];
    mm[iino] <- nbn[[ii]]$mu;
    ss[iino] <- nbn[[ii]]$sigma;
    for (jj in r.bf(nbn[[ii]]$parents)) {
      jjno <- nbn[[ii]]$parents[jj];
      rr[jjno,iino] <- nbn[[ii]]$regcoef[jj];
    }
  }
  # returning
  list(mm=mm,ss=ss,rr=rr);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn2mn <- function(nbn,algo=3)
#TITLE computes the joint distribution of a /nbn/
#DESCRIPTION 
# Computes the joint distribution of a /nbn/ with three
# possible algorithms according to \code{algo}.
#DETAILS
# To be explained if it works
#KEYWORDS 
#INPUTS
#{nbn} <<The \code{nbn} object to be converted.>>
#[INPUTS]
#{algo} <<either \code{1}: transforming the \code{nbn} into
# a \code{gema} first before getting the \code{mn} form;
# or \code{2}: one variable after another is added
# to the joint distribution following a topological order;
# or \code{3}: variances are computed  through
# the differents paths of the DAG.
#VALUE
# the resulting /mn/ object
#EXAMPLE
# print8mn(nbn2mn(rbmn0nbn.05));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_04_19
#REVISED 13_04_21
#--------------------------------------------
{
  # checking
  algo <- algo[1];
  if (!(algo %in% 1:3)) {
    warning("'nbn2mn' got a bad 'algo' argument: forced to '3'");
    algo <- 3;
  }
  if (algo==1) {
    # getting the gema expression
    gg <- nbn2gema(nbn);
    # getting the mn expression
    res <- gema2mn(gg);
  }
  #
  if (algo==2) {
    nn <- length(nbn);
    # degenerated case
    if (nn == 0) {
      mu <- numeric(0); names(mu) <- character(0);
      gamma <- matrix(NA,0,0);
      dimnames(gamma) <- list(character(0),character(0));
      res <- list(mu=mu,gamma=gamma)
    } else {
      # getting a topological order
      topo <- order4nbn(nbn);
      inio <- names(nbn);
      nbn <- nbn[topo];
      # initializing with the first
      mu <- nbn[[1]]$mu;
      names(mu) <- names(nbn)[1];
      gamma <- matrix(nbn[[1]]$sigma^2,1,1);
      dimnames(gamma) <- list(names(mu),names(mu));
      # following for all remaining nodes
      for (ii in r.bd(2,nn)) {
        mun <- nbn[[ii]]$mu;
        gan <- c(rep(0,ii-1),nbn[[ii]]$sigma^2);
        if (length(nbn[[ii]]$parents) > 0) {
          # parents have to be considered
          ou <- match(nbn[[ii]]$parents,names(nbn)[r.bc(ii-1)]);
          coe <- rep(0,ii-1);
          coe[ou] <- nbn[[ii]]$regcoef;
          mun <- mun + sum(mu*coe);
          mu <- c(mu,mun);
          gg <- cbind(gamma,rep(0,ii-1));
          gg <- rbind(gg,gan);
          coef <- cbind(diag(nrow=ii-1,ncol=ii-1),rep(0,ii-1));
          coef <- rbind(coef,c(coe,1));
          gamma <- coef %*% gg %*% t(coef);
        } else {
          # it is a root node
          mu <- c(mu,mun);
          gamma <- cbind(gamma,rep(0,ii-1));
          gamma <- rbind(gamma,gan);
        }
        # naming
        names(mu) <- names(nbn)[r.bc(ii)];
        dimnames(gamma) <- list(names(mu),names(mu));
      }
      # coming back to the initial order
      mu <- mu[inio];
      gamma <- gamma[inio,inio];
      # building the object
      res <- list(mu=mu,gamma=gamma);
    }
  }
  #
  if (algo==3) {
    # getting the  conditional expectation matrix
    nn <- length(nbn);
    cmu <- rep(NA,nn);
    for (ii in r.bc(nn)) {
      cmu[ii] <- nbn[[ii]]$mu;
    }
    # (1) computing the rr matrix
    deco <- nbn2rr(nbn);
    rr <- deco$rr;
    # (2) computing all the paths
    pp <- uu <- rr;
    for (ii in r.bc(nn-1)) {
      uu <- uu%*%rr;
      if (sum(uu)==0) { break}
      pp <- pp + uu;
    }
    # (3) adding the diagonal of ones
    pp <- pp + diag(nrow=nn,ncol=nn);
    # (3b) getting the marginal expectation
    mu <- t(pp) %*% cmu;
    names(mu) <- names(nbn);
    # (4) multiplying with the st.dev.
    pp <- diag(deco$ss,nrow=nn) %*% pp;
    # (5) getting the variance
    va <- t(pp) %*% pp;
    dimnames(va) <- list(names(nbn),names(nbn));
    # constituting the /mn/ object
    res <- list(mu=mu,gamma=va);
  }
  #
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

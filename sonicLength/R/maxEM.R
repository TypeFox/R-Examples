maxEM <- function(slmat, theta.var=FALSE, phi.update=NULL,
                  phi.deriv=NULL, lframe = NULL, glm.frm=NULL, iter.control=NULL, ...) {
  ## Purpose: maximum likelihood estimates of theta and phi
  ## ----------------------------------------------------------------------

  
  ## Arguments: slmat - t[l,ijk] = 0/1 - a matrix whose rows may
  ##                    correspond to unique lengths with rownames
  ##                    indicating those lengths, or a data.frame will
  ##                    be provided to guide estimation of phi

  ##            theta.var - logical, return variance of theta estimates?
  
  ##            phi.update - function of a matrix like slmat which
  ##                            returns a vector of estimates of phi

  
  ##            phi.deriv - function of theta and phi that returns
  ##                        derivatives of phi wrt beta (its
  ##                        parameters)

  #             lframe - a data.frame which will be used to estimate
  #                     phi, the supplied function phi.update must
  #                     know what to do with this arg or it will be
  #                     ignored
  #             glm.frm - a formula to use in fitting phi

  ##            iter.control - a list, such as returned by
  ##                           maxEM.iter.control()

  ##            ... - possibly other args to pass to phi.update
  
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:  1 May 2011, 14:09

  mc <- match.call()
  if (is.null(phi.update)) phi.update <- phi.update.default
  if (theta.var && is.null(phi.deriv)) phi.deriv <- phi.deriv.default
  if (is.null(iter.control)) iter.control <- maxEM.iter.control()
  ## phi NOT Constrained - beware of logLik that this might yield:
  phi.new <- phi.old <- prop.table(rowSums(slmat+1/length(slmat))) 
  theta.new <- theta.old <- Yplus <- colSums(slmat)
  ## do one EM step to start
  ## phi.update assures that constraints on phi are met
  Y <- Ey.given.x(slmat,theta.new,phi.new)
  theta.new <- theta.old <- colSums(Y)
  phi.new <- phi.old <- prop.table(phi.update(Y, lframe=lframe , glm.frm = glm.frm )+iter.control$phi.min)
  llk.old <-  -Inf
  i <- 0 
  not.done <- TRUE
  while ( not.done ){
    mres <- mstep( slmat, theta.new, phi.new )
    theta.new <- mres$theta
    Y <- Ey.given.x(slmat,theta.new,phi.new)
    llk <- attr(Y,"logLik")
    ## check increase in loglik - use one full EM step if llk doesn't increase
    fail <- 0L
    if (llk>llk.old){
      llk.old <- llk
      theta.old <- theta.new
    } else {
      fail <- fail+1
      ##      cat('reset theta\n')
      theta.new <- theta.old
      Y <- Ey.given.x(slmat, theta.new,phi.new)
    }
    ## a little padding is needed to guard against numerical zeroes in phi 
    phi.new <- prop.table(phi.update(Y, lframe=lframe)+iter.control$phi.min)
    Y <- Ey.given.x(slmat, theta.new,phi.new)
    llk <- attr(Y,'logLik')
    if (llk>llk.old){
      llk.old <- llk
      phi.old <- phi.new
    } else {
      ## cat('reset phi\n')
      phi.new <- phi.old
      fail <- fail+1
    }
    ## print(llk,digits=20)

    adjs <- mres$grad
    max.abs <- max(abs(adjs))
    max.rel <- max(abs(adjs/theta.old))
    i <- i+1
    not.done <- (fail<2) && ((i < iter.control$min.reps) || (max.abs > iter.control$max.abs.le) || (max.rel > iter.control$max.rel.le)) # defaults 3, 0.1, 1e-05
    if ( i >= iter.control$max.reps ) warning( "iteration limit reached")
  }
  theta.var <-
    if (theta.var){
      dphi.db <- phi.deriv( theta.new, phi.new, phi.update, lframe=lframe, ... )
      blks <- info.blocks( theta.new, phi.new, dphi.db ) 
      var.theta( blks , theta.new )
    } else {
      NULL
    }
  list(theta=theta.new, phi=phi.new, var.theta=theta.var, iter=i, call=mc)
}





maxEM.iter.control <- function(min.reps=3, max.reps=2000, max.abs.le=0.1, max.rel.le=1e-05, phi.min=.Machine$double.eps)
{
  ## Purpose: iteration settings
  ## ----------------------------------------------------------------------
  ## Arguments: min.reps - default 3,
  ##            max.reps - default 2000,
  ##            max.abs.le - default 0.1,
  ##            max.rel.le - default 1e-05
  ##            phi.min - default sqrt(.Machine$double.eps)
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:  1 May 2011, 14:32

  if (!is.numeric(phi.min) || exp(-phi.min)==1.0 )
    stop(' phi.min too small ')  
  
  if (!is.numeric(min.reps) || (min.reps>max.reps) )
    stop(' min.reps must be less than max.reps ')
  
  if (!is.numeric(max.reps) || (0>max.reps) )
    stop(' max.reps must be positive ')
  
  if (!is.numeric(max.abs.le) || (max.abs.le<0) )
    stop(' max.abs.le must be gt 0 ')

  if (!is.numeric(max.rel.le) || (max.rel.le<0) )
    stop(' max.rel.le must be gt 0 ')
  
  if (!is.numeric(min.reps) || (min.reps>max.reps) )
    stop(' min.reps must be less than max.reps ')
  
  list( min.reps=min.reps, max.reps=max.reps, max.abs.le=max.abs.le,
        max.rel.le=max.rel.le, phi.min=phi.min)
  
}


phi.update.default <- function(obj, return.fit=FALSE, glm.frm=NULL, ... )
{
  ## Purpose: estimate phi
  ## ----------------------------------------------------------------------
  ## Arguments: obj - a matrix like formals(maxEM)$slmat
  ##            return.fit - logical, return just the fitted glm()??
  ##            ... - other args will be ignored
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:  1 May 2011, 14:56
  tmp.frame <- pad.tab(rowSums(obj))
   if ( is.null(glm.frm) ) glm.frm <- y~bs(x,knots=c(50,100)) 
  fit <- glm( glm.frm ,tmp.frame,family=quasipoisson,control=list(epsilon=1e-12))
  if (return.fit) {
    fit
  } else {
    prop.table(exp(predict(fit)[tmp.frame$orig]))
  }
}

phi.update.lframe <- function(obj, return.fit=FALSE, lframe, glm.frm=NULL, ...)
{
  ## Purpose: estimate phi flexibly
  ## ----------------------------------------------------------------------
  ## Arguments: obj - a matrix like formals(maxEM)$slmat
  ##            return.fit - logical, return just the fitted glm()??

  ##            lframe - data.frame with columns y, x, strata, orig and
  ##                     rownames are a superset of the rownames of
  ##                     slmat

  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 26 May 2011, 12:03

  ## if obj is NOT supplied assume that lframe$y is all set up,
  ## otherwise get the values from obj
  if (!missing(obj)) lframe[rownames(obj),'y'] <- rowSums(obj)
  if ( is.null(glm.frm) ){
    lframe$strata <- factor(lframe$strata)
    glm.frm <-
      if (length(levels(lframe$strata))==1){
        y ~ bs( x, knots= c(50,100) )
      } else {
        y ~ factor(strata):bs( x, knots= c(50,100) )
      }
  }
  fit <- glm( glm.frm, lframe, family=quasipoisson,control=list(epsilon=1e-12) )
  if (return.fit) {
    fit
  } else {
    prop.table(exp(predict(fit)[lframe$orig]))
  }
}

phi.deriv.default <- function(theta,phi, phi.update=phi.update.default, ... )
{
  ## Purpose: deriv of phi w.r.t coef(fit)
  ## ----------------------------------------------------------------------
  ## Arguments: theta - like maxEM$theta
  ##            phi -   like maxEM$phi
  ##            phi.update - function to update phi
  ##            ... - possibly other args for phi.update
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 13 May 2011, 13:30

  fit <- phi.update( exp( -phi %o% theta ), return.fit=TRUE, ... )
  xphi <- (phi ) * model.matrix( fit, subset=fit$data$orig )
  res <- xphi - phi %o% colSums( xphi )
  res[, apply(res,2,sd) > .Machine$double.eps ]
}

EdL2.dlambda <- function(theta,phi)
{
  ## Purpose: Fisher Info
  ## ----------------------------------------------------------------------
  ## Arguments: theta - like maxEM$theta
  ##            phi -   like maxEM$phi
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 14 May 2011, 08:35
  lambda <- phi %o% theta
  expnl <- exp(-lambda)
  expnl/(1-expnl)
}

info.blocks <- function( theta, phi, dphi.dbeta )
{
  ## Purpose: Information Matrix in compact form
  ## ----------------------------------------------------------------------
  ## Arguments: theta - like maxEM$theta
  ##            phi -   like maxEM$phi
  ##            dphi.dbeta - the partial of ph w.r.t. beta
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 14 May 2011, 08:40

  info.lambda <- EdL2.dlambda( theta, phi )
  ## off diagonal elements are zero for theta
  info.theta <- as.vector( phi^2 %*% info.lambda )
  ## WAS: dphi.dbeta <- phi.deriv( theta, phi )
  nc <- ncol( dphi.dbeta )
  info.theta.beta <- ( t( dphi.dbeta * phi ) %*% info.lambda ) * rep(theta,each=nc)
  dimnames(info.theta.beta) <- NULL
  xxp <- dphi.dbeta[ , rep(1:nc,nc)] * dphi.dbeta[ , rep(1:nc,each=nc)]
  info.beta <- t( xxp ) %*% ( info.lambda %*% theta^2 )
  dim( info.beta ) <- c( nc, nc )
  list(theta=info.theta, theta.beta=info.theta.beta,beta=info.beta)
}

var.theta <- function( blocks, theta)
{
  ## Purpose: parts of the inverse the result of info.blocks()
  ## ----------------------------------------------------------------------
  ## Arguments: blocks - a list as returned by info.blocks
  ##            theta - the current estimate
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 14 May 2011, 09:36

  ## this function uses a block inverse approach, which saves
  ## computation as the block for theta is a diagonal matrix. Further,
  ## since the desired results are just the diagonals of var(theta) or
  ## var(p), these can be built up from the expression for the block
  ## inverse
  Q <- with( blocks,  beta - theta.beta %*% ( t( theta.beta ) / theta ) )
  Q.inv.b <-  solve(Q, blocks$theta.beta)
  diag.b.Q.inv.b <- colSums( blocks$theta.beta*Q.inv.b )
  res <- with( blocks,
              1/theta + diag.b.Q.inv.b / theta^2
              )
  ## get the variance of the proportion theta/sum(theta). the matrix
  ## of derviatives of p wrt theta is the sum of an identity matrix
  ## times a constant (say P1) and a matrix whose columns are all the
  ## same (say P2). so multiplications by the derivatives can be
  ## simplified.

  diag.element <- 1/sum(theta)
  column <- - diag.element^2 * theta 
  ## (P1+P2) %*% Dinv %*% t(P1+P2)
  res.p <- diag.element^2/blocks$theta # P1 %*% Dinv %*% t(P1)
  res.p <- res.p + column^2 * sum( 1 / blocks$theta )
  res.p <- res.p + 2 * column * diag.element / blocks$theta

  ## (P1+P2) %*% Dinv %*% ( t(B) %*% Qinv %*% B ) Dinv %*% t(P1+P2)

  P1.etc.P1 <- diag.element^2*diag.b.Q.inv.b / blocks$theta^2
  dvec.B <- with(blocks, theta.beta %*% (1/theta) )
  P1.etc.P2 <- diag.element* column * (t(dvec.B) %*% Q.inv.b)/blocks$theta 
  P2.etc.P2 <- column^2 * t(dvec.B) %*% solve(Q,dvec.B)
  res.p <- res.p + P1.etc.P1 + 2* P1.etc.P2 + P2.etc.P2
  dim(res.p) <- NULL
  list(theta=res,prop=res.p)
}

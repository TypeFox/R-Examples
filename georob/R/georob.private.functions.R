####################################
#                                  #
#   Hilfsfunktionen fuer georob    #
#                                  #
####################################

##  ##############################################################################

covariances.fixed.random.effects <- 
  function(
    Valphaxi.objects, 
    Aalphaxi, Palphaxi, Valphaxi.inverse.Palphaxi,
    rweights, XX, TT, TtT, names.yy,
    nugget, eta, 
    expectations, family = c( "gaussian", "long.tailed" ),
    cov.bhat, full.cov.bhat,
    cov.betahat, 
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    control.pmm,
    verbose
  )
{
  
  ##  ToDos:
  
  ##  function computes the covariance matrices of 
  ##  - bhat
  ##  - betahat
  ##  - bhat and betahat
  ##  - delta.b = b - bhat
  ##  - delta.b and betahat
  ##  - residuals ehat = y - X betahat - bhat
  ##  - residuals ehat.p.bhat = y - X betahat = ehat + bhat
  ##  - auxiliary matrix to compute covariance between kriging predictions of
  ##    y and y
  
  ## 2011-10-13 A. Papritz
  ## 2011-12-14 AP modified for replicated observations
  ## 2012-02-23 AP checking new variant to compute covariances of betahat and bhat
  ## 2012-04-27 AP scaled psi-function
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-11-04 AP unscaled psi-function
  ## 2013-02-05 AP covariance matrix of zhat
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-19 AP correcting error when computing covariance of regression residuals
  ## 2015-03-03 AP correcting error in computing result.new[["cov.bhat"]] for full.cov.bhat == FALSE
  ## 2015-06-30 AP changes for improving efficiency
  ## 2015-07-17 AP TtT passed to function, new name for function
  ## 2015-08-19 AP changes for computing covariances under long-tailed distribution of epsilon
  
  family = match.arg( family )
  
  ## store flags for controlling output
  
  ccov.bhat               <- cov.bhat
  ffull.cov.bhat          <- full.cov.bhat
  ccov.betahat            <- cov.betahat
  ccov.bhat.betahat       <- cov.bhat.betahat
  ccov.delta.bhat         <- cov.delta.bhat
  ffull.cov.delta.bhat    <- full.cov.delta.bhat
  ccov.delta.bhat.betahat <- cov.delta.bhat.betahat
  
  ## setting flags for computing required items
  
  cov.ystar        <- FALSE
  cov.bhat.b       <- FALSE
  cov.bhat.e       <- FALSE
  cov.betahat.b    <- FALSE
  cov.betahat.e    <- FALSE
  
  if( aux.cov.pred.target ){
    cov.bhat.b <- cov.betahat.b <- TRUE
  }
  
  if( cov.ehat.p.bhat ){
    cov.betahat <- cov.betahat.b <- cov.betahat.e <- TRUE
  }
  
  if( cov.ehat ){
    cov.delta.bhat.betahat <- cov.delta.bhat <- cov.betahat <- cov.bhat.e <- cov.betahat.e <- TRUE
    if( full.cov.ehat ) full.cov.delta.bhat <- TRUE
  }
  
  if( cov.delta.bhat.betahat ){
    cov.betahat.b <- cov.bhat.betahat <- TRUE
  }
  
  if( cov.delta.bhat ){
    cov.bhat <- cov.bhat.b <- TRUE
    if( full.cov.delta.bhat ) full.cov.bhat <- TRUE
  }
  
  if( any( cov.bhat, cov.betahat, cov.bhat.betahat, cov.delta.bhat, cov.delta.bhat.betahat, cov.ehat ) ){
    cov.ystar <- TRUE
  }
  
  ## compute required auxiliary items 
  
  switch(
    family,
    gaussian = {
      var.psi     <- expectations["var.gauss.psi"] 
      exp.dpsi    <- expectations["exp.gauss.dpsi"]
      var.eps     <- nugget
      cov.psi.eps <- expectations["exp.gauss.dpsi"]
    },
    long.tailed = {
      var.psi     <- expectations["var.f0.psi"]
      exp.dpsi    <- expectations["var.f0.psi"]
      var.eps     <- nugget * expectations["var.f0.eps"]
      cov.psi.eps <- 1.
    }
  )
  
  V <- eta * nugget * Valphaxi.objects[["Valphaxi"]] 
  VTtT <- t( TtT * V )
  
  result.new <- list( error = FALSE )
  
  #   ## auxiliary items for checking computation for Gaussian case
  #   
  #   cov.bhat = TRUE
  #   full.cov.bhat = TRUE
  #   cov.betahat = TRUE 
  #   cov.bhat.betahat = TRUE
  #   cov.delta.bhat = TRUE
  #   full.cov.delta.bhat = TRUE
  #   cov.delta.bhat.betahat = TRUE
  #   cov.ehat = TRUE
  #   full.cov.ehat = TRUE
  #   cov.ehat.p.bhat = TRUE
  #   full.cov.ehat.p.bhat = TRUE
  #   aux.cov.pred.target = TRUE
  #   
  #   Gammat <- VTtT
  #   diag( Gammat ) <- diag( Gammat ) + nugget
  #   Gammati <- solve( Gammat )
  #   VGammati <- VTtT %*% Gammati
  #   tmp1 <- solve( crossprod( XX, Gammati ) %*% XX )
  #   tmp2 <- XX[TT,] %*% tmp1 %*% t(XX[TT,])
  
  ## compute S_alphaxi
  
  aux <- Valphaxi.inverse.Palphaxi / ( exp.dpsi * eta )
  diag( aux ) <- diag( aux ) + TtT
  aux <- try( chol( aux ), silent = TRUE )
  if( identical( class( aux ), "try-error" ) ){
    result.new[["error"]] <- TRUE
    return( result.new )            
  }
  Salphaxi <- chol2inv( aux )
  
  ## factors to compute bhat and betahat from zhat
  
  if( any( c( cov.ystar, cov.bhat.b, cov.bhat.e, cov.bhat, cov.bhat.betahat ) ) ){
    PaSa <- pmm( Palphaxi, Salphaxi, control.pmm )
  }
  
  if( any( c( cov.betahat.b, cov.betahat.e, cov.betahat, cov.bhat.betahat) ) ){
    AaSa <- Aalphaxi %*% Salphaxi
  }
  
  ## covariance of huberized observations
  
  if( cov.ystar ){
    cov.ystar <- TtT * VTtT
    diag( cov.ystar ) <- diag( cov.ystar ) + (var.psi * nugget / exp.dpsi^2) * TtT
    PaSa.cov.ystar <- pmm( PaSa, cov.ystar, control.pmm )
  }
  
  ## covariance of bhat and betahat with B and epsilon
  
  if( cov.bhat.b )    cov.bhat.b      <- pmm( PaSa, t( VTtT ), control.pmm )
  if( cov.bhat.e )    cov.bhat.e      <- (nugget / exp.dpsi * cov.psi.eps * PaSa)[, TT]
  if( cov.betahat.b ){
    cov.betahat.b <- tcrossprod( AaSa, VTtT )
    #     cov.betahat.b <- AaSa %*% t( VTtT )
    TX.cov.betahat.bT <- (XX %*% cov.betahat.b)[TT,TT]
  }
  if( cov.betahat.e ){
    cov.betahat.e <- (nugget / exp.dpsi * cov.psi.eps * AaSa)[, TT]
    TX.cov.betahat.e <- (XX %*% cov.betahat.e)[TT,]
  }
  
  ## compute now the requested covariances ...
  
  ## ... of bhat (debugging status ok)

  if( cov.bhat ){
    t.cov.bhat <- if( full.cov.bhat )
    {
      aux <- pmm( PaSa.cov.ystar, t(PaSa ), control.pmm )
      attr( aux, "struc" ) <- "sym"
      aux
    } else {
      aux <- rowSums( PaSa.cov.ystar * PaSa )
      names( aux ) <- rownames( XX )
      aux
    }
    if( ccov.bhat ) result.new[["cov.bhat"]] <- if( ffull.cov.bhat ){
      t.cov.bhat
    } else {
      f.diag( t.cov.bhat )
    }
  }
  
  #   print( 
  #     summary( 
  #       c( 
  #         t.cov.betahat - 
  #         (VGammati %*% VTtT - VGammati %*% tmp2 %*% t( VGammati ))
  #       )
  #     )
  #   )
  
  ## ... of betahat (debugging status ok)
  
  if( cov.betahat ){
    t.cov.betahat <- tcrossprod( tcrossprod( AaSa, cov.ystar ), AaSa )
    attr( t.cov.betahat, "struc" ) <- "sym"
    if( ccov.betahat ) result.new[["cov.betahat"]] <- t.cov.betahat
  }
  
  #   print( summary( c( 
  #         t.cov.betahat - 
  #         tmp1
  #       )
  #     )
  #   )
  
  ##  ... of bhat and betahat (debugging status ok)
    
  if( cov.bhat.betahat ){
    t.cov.bhat.betahat <- tcrossprod( PaSa.cov.ystar, AaSa )
    if( ccov.bhat.betahat ) result.new[["cov.bhat.betahat"]] <- t.cov.bhat.betahat
  }
  
  #   print( summary( c( t.cov.bhat.betahat ) ) )
  
  ## ... of (b - bhat) (debugging status ok)
  
  if( cov.delta.bhat ){
    t.cov.delta.bhat <- if( full.cov.delta.bhat ){
      aux <- V + t.cov.bhat - cov.bhat.b - t( cov.bhat.b )
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( rownames( XX ), rownames( XX ) )
      aux
    } else {
      #       aux <- diag( V ) - 2 * diag( cov.bhat.b ) + (if( full.cov.bhat ){
      #         diag( t.cov.betahat )
      #       } else {
      #         t.cov.betahat
      #       })
      aux <- diag( V ) - 2 * diag( cov.bhat.b ) + f.diag( t.cov.bhat )
      names( aux ) <- rownames( XX )
      aux
    }
    if( ccov.delta.bhat ) result.new[["cov.delta.bhat"]] <- if( ffull.cov.delta.bhat ){
      t.cov.delta.bhat
    } else {
      f.diag( t.cov.delta.bhat )
    }
  }
  
  #   print( 
  #     summary( 
  #       c( 
  #         t.cov.delta.bhat - 
  #         (VTtT - VGammati %*% VTtT + VGammati %*% tmp2 %*% t( VGammati ))
  #       )
  #     )
  #   )
  
  ## ... of (b - bhat) and betahat (debugging status ok)
  
  if( cov.delta.bhat.betahat ){
    t.cov.delta.bhat.betahat <- t( cov.betahat.b ) - t.cov.bhat.betahat
    dimnames( t.cov.delta.bhat.betahat ) <- dimnames( XX )
    if( ccov.delta.bhat.betahat ){
      result.new[["cov.delta.bhat.betahat"]] <- t.cov.delta.bhat.betahat
    }
  }
  
  #   print( 
  #     summary( 
  #       c( 
  #         t.cov.delta.bhat.betahat - 
  #         (VGammati %*% XX %*% tmp1)
  #       )
  #     )
  #   )
  
  ## ... of ehat (debugging status ok)
 
  if( cov.ehat ){
    aux1 <- tcrossprod( t.cov.delta.bhat.betahat, XX )[TT,TT]
    result.new[["cov.ehat"]] <- if( full.cov.ehat )
    {
      aux <- t.cov.delta.bhat[TT,TT] + 
        tcrossprod( tcrossprod( XX, t.cov.betahat ), XX )[TT,TT] -
        aux1 - t(aux1) - cov.bhat.e[TT,] - t(cov.bhat.e)[,TT] - 
        TX.cov.betahat.e - t(TX.cov.betahat.e)
      diag( aux ) <- diag( aux ) + var.eps
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux   
    } else {
      #       aux <- (if( full.cov.delta.bhat ){
      #         diag( t.cov.delta.bhat )[TT] 
      #       } else {
      #         t.cov.delta.bhat[TT]
      #       }) + rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] -
      #         2 * diag( aux1 ) - 2 * diag( cov.bhat.e[TT,] ) - 2 * diag( TX.cov.betahat.e ) + 
      #         nugget
      aux <- f.diag( t.cov.delta.bhat )[TT] +
        rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] -
        2 * diag( aux1 ) - 2 * diag( cov.bhat.e[TT,] ) - 2 * diag( TX.cov.betahat.e ) + 
        var.eps
        #         t.cov.delta.bhat[TT]
        #       }) + rowSums( XX * (XX %*% t.cov.betahat) )[TT] -
        #         2 * diag( aux1 ) - 2 * diag( cov.bhat.e[TT,] ) - 2 * diag( TX.cov.betahat.e ) + 
        #         nugget
        names( aux ) <- names.yy
      aux
    }
  }
  
  #   tmp3 <- -VGammati
  #   diag(tmp3) <- diag(tmp3) + 1
  #   print( 
  #     summary( 
  #       c( 
  #         result.new[["cov.ehat"]] - 
  #         (tmp3 %*% (Gammat - tmp2 ) %*% tmp3)
  #       )
  #     )
  #   )
  
  
  ## ... of ehat + bhat (debugging status ok)
  
  if( cov.ehat.p.bhat ){
    result.new[["cov.ehat.p.bhat"]] <- if( full.cov.ehat.p.bhat )
    {
      aux <- tcrossprod( tcrossprod( XX, t.cov.betahat ), XX )[TT,TT] - 
        TX.cov.betahat.bT - t(TX.cov.betahat.bT) -
        TX.cov.betahat.e - t(TX.cov.betahat.e) + V[TT,TT]
      diag( aux ) <- diag( aux ) + var.eps
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux   
    } else {
      aux <- rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] - 
        2 * diag( TX.cov.betahat.bT ) - 
        2 * diag( TX.cov.betahat.e ) + diag( V )[TT] + var.eps
#       aux <- rowSums( XX * (XX %*% t.cov.betahat) )[TT] - 
#         2 * diag( TX.cov.betahat.bT ) - 
#         2 * diag( TX.cov.betahat.e ) + diag( V )[TT] + nugget
      names( aux ) <- names.yy
      aux
    }
  }
    
  #   print( 
  #     summary( 
  #       c( 
  #         result.new[["cov.ehat.p.bhat"]] - 
  #         (Gammat - tmp2)
  #       )
  #     )
  #   )
  
  ## ...  auxiliary item to compute covariance of kriging predictions
  ## and observations
  
  if( aux.cov.pred.target ){
    result.new[["cov.pred.target"]] <- pmm(
      rbind( cov.bhat.b, cov.betahat.b ),
      Valphaxi.objects[["Valphaxi.inverse"]] / eta / nugget,
      control.pmm
    )
  }
  
  return( result.new )
  
}

##   ##############################################################################

update.zhat <- 
  function( 
    XX, yy, res, TT, 
    nugget, eta, reparam,
    Valphaxi.inverse.Palphaxi,
    psi.function, tuning.psi, 
    verbose
  )
{
  
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-29 AP solving linear system of equation by cholesky decomposition
  ## 2015-12-02 AP correcting error in try(chol(m))
  
  ## function computes (1) updated IRWLS estimates zhat of linearized
  ## normal equations, (2) the associated rweights,
  ## (3) the unstandardized residuals (= estimated epsilons); the results
  ## are returned as a list
  
  ## compute rweights (cf. p. 7, proposal HRK of 26 July 2010)

  ## 2013-04-23 AP new names for robustness weights
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram
  
  if( reparam ){
    
    Wi <- rep( 1., length( res ) )
    
  } else {
    
    std.res <- res / sqrt( nugget )
    
    ##  construct left-hand side matrix M and right-hand side vector of
    ##  linear equation system
    
    Wi <- ifelse( 
      abs( std.res ) < sqrt( .Machine[["double.eps"]] ),
      1.,
      psi.function( std.res, tuning.psi ) / std.res
    )
    
  }
  
  ##  aggregate rweights for replicated observations
  
  if( sum( duplicated( TT ) ) > 0 ){
    
    TtWiT  <- as.vector( tapply( Wi, factor( TT ), sum ) )
    TtWiyy <- as.vector( tapply( Wi * yy, factor( TT ), sum ) )
    
  } else {
    
    TtWiT <- Wi
    TtWiyy <- Wi * yy
    
  }
  
  ##  construct left-hand side matrix M and right-hand side vector b of
  ##  linearized system of equations
    
  M <- Valphaxi.inverse.Palphaxi / eta
  diag( M ) <- diag( M ) + TtWiT
  
  b <- TtWiyy
  
  ##  solve linear system
  
  result <- list( error = TRUE )
  
#   r.solve <- try( solve( M, b ), silent = TRUE ) 
#   
#   if( !identical( class( r.solve ), "try-error" ) ) {
    
  t.chol <- try( chol( M ), silent = TRUE )
  if( !identical( class( t.chol ), "try-error" ) ){
    
    r.solve <- forwardsolve( t( t.chol ), b )
    r.solve <- backsolve( t.chol, r.solve )
    
    
    ##  collect output
    
    result[["error"]]      <- FALSE
    result[["zhat"]]      <- r.solve
    result[["residuals"]]  <- yy - result[["zhat"]][TT]
    result[["rweights"]]   <- Wi
    
  }
  
  return( result )
  
}


##    ##############################################################################

estimating.equations.B <- function( 
  res, TT, zhat, 
  nugget, eta, reparam,
  Valphaxi.inverse.Palphaxi, 
  psi.function, tuning.psi
){
  
  ## auxiliary function to compute estimating equations for zhat

  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram

  if( reparam ){
    
    Ttpsi <- res
    if( sum( duplicated( TT ) > 0 ) ){
      Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    }
    
  } else {
    
    Ttpsi <- sqrt( nugget ) * psi.function( res / sqrt( nugget ), tuning.psi )
    if( sum( duplicated( TT ) > 0 ) ){
      Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    }
    
  }
      
  Ttpsi - drop( Valphaxi.inverse.Palphaxi %*% zhat ) / eta

}

##    ##############################################################################

estimate.zhat <- 
  function(
    compute.zhat,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, 
    maxit, ftol,
    nugget, eta, reparam,
    Valphaxi.inverse,
    control.pmm,
    verbose
  )
{
  
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-06-30 AP new method to determine convergence
  ## 2015-07-17 AP computing residual sums of squares (UStar)
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram
   
  ## function computes (1) estimates zhat, bhat, betahat by
  ## solving robustified estimating equations by IRWLS,
  ## (2) the weights of the IRWLS, (3) the unstandardized residuals
  ## (= estimated epsilons); the results are returned as a list
  
  ##  compute projection matrix Palphaxi and related items
  
  result <- list( error = FALSE )
  
  aux <- crossprod( XX, Valphaxi.inverse )
  
  if( col.rank.XX[["deficient"]] ){
    s <- svd( aux %*% XX )
    s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    Palphaxi <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )                # Moore-Penrose inverse
  } else {
    t.chol <- try( chol( aux %*% XX ), silent = TRUE )
    if( !identical( class( t.chol ), "try-error" ) ){
      Palphaxi <- chol2inv( t.chol )    
    } else {
      result[["error"]] <- TRUE
      return( result )    
    }
  }
  
  result[["Aalphaxi"]]             <- Palphaxi %*% aux
  dimnames( result[["Aalphaxi"]] ) <- dimnames( t(XX) )
  
  result[["Palphaxi"]]             <- -XX %*% result[["Aalphaxi"]]
  diag( result[["Palphaxi"]] )     <- diag( result[["Palphaxi"]] ) + 1.
  rownames( result[["Palphaxi"]] ) <- rownames( XX )
  colnames( result[["Palphaxi"]] ) <- rownames( XX )
  
  result[["Valphaxi.inverse.Palphaxi"]] <- pmm( 
    Valphaxi.inverse, result[["Palphaxi"]], control.pmm 
  )
  rownames( result[["Valphaxi.inverse.Palphaxi"]] )      <- rownames( XX )
  colnames( result[["Valphaxi.inverse.Palphaxi"]] )      <- rownames( XX )
  attr( result[["Valphaxi.inverse.Palphaxi"]], "struc" ) <- "sym"
  
  if( compute.zhat ){
      
    res <- yy - zhat[TT]
    
    eeq.old <- estimating.equations.B(     
      res, TT, zhat, 
      nugget, eta, reparam,
      result[["Valphaxi.inverse.Palphaxi"]], 
      psi.function, tuning.psi
    )
    eeq.old.l2 <- sum( eeq.old^2 ) / 2.
    
    if( !is.finite( eeq.old.l2 ) ) {
      result[["error"]] <- TRUE
      return( result )
    }
    
    converged <- FALSE
    
    if( verbose > 2 ) cat(
      "\n  IRWLS\n",
      "      it     Fnorm.old     Fnorm.new   largest |f|\n", sep = ""
    )
    
    ##  IRWLS
    
    for( i in 1:maxit ){
      
      ##  compute new estimates 
      
      new <- update.zhat(
        XX, yy, res, TT, 
        nugget, eta, reparam,
        result[["Valphaxi.inverse.Palphaxi"]],
        psi.function, tuning.psi, 
        verbose
      )
      
      if( new[["error"]] ) {
        result[["error"]] <- TRUE
        return( result )
      }
      
      ##  evaluate estimating equations for xi and compute its l2 norm
      
      eeq.new <- estimating.equations.B(       
        new[["residuals"]], TT, new[["zhat"]], 
        nugget, eta, reparam,
        result[["Valphaxi.inverse.Palphaxi"]], 
        psi.function, tuning.psi
      )
      eeq.new.l2 <- sum( eeq.new^2 ) / 2.
      
      if( !is.finite( eeq.new.l2 ) ) {
        result[["error"]] <- TRUE
        return( result )
      }
      
      if( verbose > 2 ) cat( 
        format( i, width = 8 ),
        format( 
          signif( 
            c( eeq.old.l2, eeq.new.l2, max(abs(eeq.new) ) ), digits = 7 
          ), scientific = TRUE, width = 14 
        ), "\n", sep = ""
      )
      
      ##  check for convergence
      
      if( max( abs( eeq.new ) ) <= ftol && i >= 2L ){
        converged <- TRUE
        break
      }
      
      ##  update zhat, residuals and eeq.old.l2
      
      eeq.old.l2 <- eeq.new.l2
      zhat      <- new[["zhat"]]
      res        <- new[["residuals"]]
      
    }
    
    ## compute scaled residuals sum of squares 
    
    ##  collect output
    
    result[["zhat"]]            <- new[["zhat"]]
    names( result[["zhat"]] )   <- rownames( XX )
    
    result[["residuals"]]        <- new[["residuals"]]
    result[["rweights"]]         <- new[["rweights"]]
    result[["converged"]]        <- converged
    result[["nit"]]              <- i
    
  } else {
    
    result[["zhat"]]            <- zhat
    names( result[["zhat"]] )   <- rownames( XX )
    
    result[["residuals"]]        <- yy - zhat[TT]
    
    if( reparam ){
      result[["rweights"]]       <- rep( 1., length( result[["residuals"]] ) )
    } else {
      result[["rweights"]]       <- ifelse( 
        abs( std.res <- result[["residuals"]] / sqrt( nugget ) ) < sqrt( .Machine[["double.eps"]] ),
        1.,
        psi.function( std.res, tuning.psi ) / std.res
      )
    }
    
    result[["converged"]]        <- NA
    result[["nit"]]              <- NA_integer_
    
  }
  
  result[["bhat"]]             <- drop( result[["Palphaxi"]] %*% result[["zhat"]] )
  names( result[["bhat"]] )    <- rownames( XX )
  
  result[["betahat"]]          <- drop( result[["Aalphaxi"]] %*% result[["zhat"]] )
  names( result[["betahat"]] ) <- colnames( XX )
  
  result[["Valphaxi.inverse.bhat"]] <- drop( Valphaxi.inverse %*% result[["bhat"]] )
  
  result[["RSS"]] <- f.aux.RSS(
    res = result[["residuals"]],
    TT = TT, TtT = TtT, 
    bhat = result[["bhat"]],
    Valphaxi.inverse.bhat = result[["Valphaxi.inverse.bhat"]],
    eta = eta
  )
  
  return( result )
  
}


##    ##############################################################################

f.aux.gcr <- 
  function( 
    lag.vectors, variogram.model, param, xi, aniso, 
    irf.models = control.georob()[["irf.models"]],
    control.pmm, verbose
  )
{
  
  ##  Function computes the generalized correlation (matrix) for the lag
  ##  distances in lag.vectors.  The result is a generalized correlation matrix
  ##  that is positive definite.
  
  ##  cf. HRK's notes of 2011-06-17 on "Robust Kriging im intrinsischen
  ##  Fall"
  
  ##  2011-12-27 ap
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-03-05 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-23 AP Valpha (correlation matrix without spatial nugget) no longer stored
  
  result <- list( error = TRUE )
  
  ## matrix for coordinate transformation
  
  A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]
  
  ## prepare model
  
  model.list <- list( variogram.model )
  model.list <- c( model.list, as.list( param[-1] ) )
  model.list <- list( "$", var = 1., A = A, model.list )

  ##  negative semivariance matrix
  
  ## functions of version 3 of RandomFields
  
  RFoptions(newAniso=FALSE)
  
  ## auxiliary function to compute generalized correlations in parallel
  
  f.aux <- function(i, s, e, lag.vectors, model.list ){
    result <- try(
      -RFvariogram(
        x = lag.vectors[s[i]:e[i], ], model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
      ),
      silent = TRUE
    )
    if( !(identical( class( result ), "try-error" ) || any( is.na( result ) )) ){
      result
    } else {
      "RFvariogram.error"
    }
  }
  
  ## definition of junks to be evaluated in parallel
  
  k <- control.pmm[["f"]] * control.pmm[["ncores"]]
  n <- NROW(lag.vectors)
  dn <- floor( n / k )
  s <- ( (0:(k-1)) * dn ) + 1
  e <- (1:k) * dn
  e[k] <- n
  
  ## compute generalized correlations in parallel
  
  if( control.pmm[["ncores"]] > 1L ){
    
    if( identical( .Platform[["OS.type"]], "windows") ){
      
      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = control.pmm[["ncores"]] )
        #         junk <- sfLibrary( RandomFields, verbose = FALSE )
        junk <- sfLibrary( georob, verbose = FALSE )
      }
      
        
      Valpha <- sfLapply( 
        1:k, f.aux, s = s, e = e, lag.vectors = lag.vectors, model.list = model.list 
      )
      
      if( control.pmm[["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }
      
    } else {
      
      Valpha <- mclapply( 
        1:k, f.aux, s = s, e = e, lag.vectors = lag.vectors, model.list = model.list,  
        mc.cores = control.pmm[["ncores"]] 
      )
      
    }
    
    not.ok <- any( sapply( Valpha, function( x ) identical( x, "RFvariogram.error" ) ) )
    
  } else {
    
    Valpha <- try(
      -RFvariogram(
        x = lag.vectors, model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
      ),
      silent = TRUE
    )
    
    not.ok <- identical( class( Valpha ), "try-error" ) || any( is.na( Valpha ) )
    
  }
  
  if( !not.ok ){
    
    Valpha <- unlist( Valpha )
        
    ## convert semivariance vectors to symmetric matrices
    
    Valpha <- list( 
      diag = rep( 0., 0.5 * ( 1 + sqrt( 1 + 8 * length( Valpha ) ) ) ),
      tri = Valpha
    )
    attr( Valpha, "struc" ) <- "sym"
    
    Valphaxi <- Valpha
    Valphaxi[["tri"]] <- (1. - xi) * ( Valphaxi[["tri"]] + 1. ) - 1.
    
    Valphaxi  <- expand( Valphaxi )
    #     Valpha  <- expand( Valpha )
    
    ##  compute additive constant for positive definiteness and
    
    if( variogram.model %in% irf.models ){
      gcr.constant <- max( -Valphaxi ) * 2.                    
    } else {
      gcr.constant <- 1.
    }
    
    ##  collect results
    
    result[["error"]]        <- FALSE
    result[["gcr.constant"]] <- gcr.constant
    #     result[["Valpha"]]       <- Valpha + gcr.constant   # correlation matrix that includes spatial nugget
    result[["Valphaxi"]]     <- Valphaxi + gcr.constant   # correlation matrix that includes spatial nugget
    
  } else {
    
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 3 ) cat(
      "\n an error occurred when computing the negative semivariance matrix\n"
    )
    
  }
  
  return( result )
    
}




##    ##############################################################################

likelihood.calculations <- 
  function(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pmm, 
    verbose
  )
{
  
  ## 2011-12-10 AP modified for replicated observations
  ## 2012-02-03 AP modified for geometrically anisotropic variograms
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-02 AP modification computing ilcf
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2012-11-27 AP changes in check allowed parameter range
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP new transformation of rotation angles
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-16 AP extended variogram parameter transformations, elimination of unused variables
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-07-17 AP TtT passed to function
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-23 AP changes for avoiding computation of Valphaxi object if not needed
  ## 2015-12-02 AP reparametrized variogram parameters renamed

  ##  function transforms (1) the variogram parameters back to their
  ##  original scale; computes (2) the correlation matrix, its inverse
  ##  and its inverse lower cholesky factor; (3) computes betahat,
  ##  bhat and further associates items; and (4) computes the
  ##  matrices A and the cholesky factor of the matrix Q
  
  d2r <- pi / 180. 
  
  ## load lik.item object
  
  lik.item <- get( "lik.item", pos = as.environment( envir ) )
  
  ##  transform variogram parameters back to original scale
  
  param <- c( adjustable.param, fixed.param )[param.name]
  
  param <- sapply(
    param.name,
    function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
    bwd.tf = bwd.tf,
    param.tf = param.tf,
    param = param
  )
  names( param ) <- param.name
  
  
  ##  transform anisotropy parameters back to original scale
  
  aniso <- c( adjustable.param, fixed.param )[aniso.name]
  
  aniso <- sapply(
    aniso.name,
    function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
    bwd.tf = bwd.tf,
    param.tf = param.tf,
    param = aniso
  )
  names( aniso ) <- aniso.name
    
  ##  check whether the current variogram parameters and the variogram
  ##  parameters that were used in the previous call to
  ##  likelihood.calculations are the same
  
  same.param <- isTRUE( all.equal( 
      c( param, aniso ), 
      c( lik.item[["param"]], lik.item[["aniso"]][["aniso"]] )
    ))
  
  if( same.param && !is.null( lik.item[["zhat"]] ) ){
    
    if( verbose > 4 ) cat(
      "\n     likelihood.calculations: exit without computing any objects\n"    
    )
    return( lik.item )
  }
  
  ## check whether variogram parameters are within reasonable bounds and
  ## return an error otherwise
  
  if( length( c( param, aniso ) ) && any( c( param, aniso ) > safe.param ) ){
    if( verbose > 1 ){
      t.param <- param
      if( !lik.item[["aniso"]][["isotropic"]] ) t.param <- c( 
        t.param, aniso / c( rep( 1., 2), rep( d2r, 3 ) ) 
      )
      cat( "\n\n                      ",
        format( names( t.param ), width = 14, justify = "right" ), 
        "\n", sep = ""
      )
      cat( "  Variogram parameters", 
        format( 
          signif( t.param, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n" , sep = ""
      )
    }
    return( lik.item )  
  }
  
  ## check whether extra variogram parameters are within allowed bounds and
  ## return an error otherwise
  
  ep <- param.names( model = variogram.model )
  param.bounds <- param.bounds( variogram.model, NCOL( lag.vectors ) )
  ep.param <- param[ep]
  
  if( !is.null( param.bounds ) ) t.bla <- sapply(
    1:length( ep.param ),
    function( i, param, bounds ){
      if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) cat(
        "value of parameter '", names( param[i] ), "' outside of allowed range", sep = "" 
      )
      return( lik.item )
    }, 
    param = ep.param,
    bounds = param.bounds
  )
  
  ##  update variogram and parameters and compute eta
  
  lik.item[["param"]] <- param
  
  if( reparam ){
    lik.item[["eta"]] <- unname( lik.item[["param"]]["nugget"] )
    lik.item[["xi"]]  <- unname( lik.item[["param"]]["snugget"] )
  } else {
    t.param <- f.reparam.fwd( lik.item[["param"]] )
    lik.item[["eta"]] <- unname( t.param["nugget"] )
    lik.item[["xi"]]  <- unname( t.param["snugget"] )
  }
  
  ##  update anisotropy parameters and the coordinate transformation
  ##  matrices
  
  lik.item[["aniso"]][["aniso"]] <- aniso
  lik.item[["aniso"]][["sincos"]] <- list(
    co = unname( cos( aniso["omega"] ) ),
    so = unname( sin( aniso["omega"] ) ),
    cp = unname( cos( aniso["phi"] ) ),
    sp = unname( sin( aniso["phi"] ) ),
    cz = unname( cos( aniso["zeta"] ) ),
    sz = unname( sin( aniso["zeta"] ) )
  )
  
  n <- NCOL( lag.vectors)
  
  if( n <= 3 ){
    
    lik.item[["aniso"]][["rotmat"]] <- with( 
      lik.item[["aniso"]][["sincos"]],
      rbind(
        c(             sp*so,             sp*co,       cp ),
        c( -cz*co + sz*cp*so,  co*sz*cp + cz*so,   -sp*sz ),
        c( -co*sz - cz*cp*so, -cz*co*cp + sz*so,    cz*sp )
      )[ 1:n, 1:n, drop = FALSE ]
    )
    
    
    lik.item[["aniso"]][["sclmat"]] <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1:n ]
    
  } else {  # only isotropic case for n > 3
    
    lik.item[["aniso"]][["rotmat"]] <- diag( n )
    lik.item[["aniso"]][["sclmat"]] <- rep( 1., n )
    
  }
  
  ## print updated variogram parameters
  
  if( verbose > 1 ) {
    
    t.param <- lik.item[["param"]]
    
    if( reparam ){
      tmp <- names( t.param )
      tmp[tmp %in% "snugget"] <- "xi"
      tmp[tmp %in% "nugget"]  <- "eta"
      names( t.param ) <- tmp  
    }
    
    if( !lik.item[["aniso"]][["isotropic"]] ) t.param <- c( 
      t.param, lik.item[["aniso"]][["aniso"]] / c( rep( 1., 2), rep( d2r, 3 ) )
    )
    cat( "\n\n                      ",
      format( names( t.param ), width = 14, justify = "right" ), 
      "\n", sep = ""
    )
    cat( "  Variogram parameters", 
      format( 
        signif( t.param, digits = 7 ), 
        scientific = TRUE, width = 14
      ), "\n" , sep = ""
    )
  }
  
  
  ##  compute updates of required likelihood items if the variogram
  ##  parameters differ and are all within allowed bounds
  
  if( !same.param || is.null( lik.item[["Valphaxi"]] ) ){
  
    if( verbose > 4 ) cat(
      "\n     likelihood.calculations: computing 'Valphaxi' object\n"    
    )
    
    lik.item[["Valphaxi"]][["error"]] <- TRUE
    
    ##  calculate generalized correlation matrix, its inverse and its
    ##  inverse cholesky factor
    
    t.Valphaxi <- f.aux.Valphaxi(
      lag.vectors = lag.vectors, variogram.model = variogram.model,
      param = lik.item[["param"]][!names(lik.item[["param"]]) %in% c( "variance", "snugget", "nugget")],
      xi = lik.item[["xi"]],
      aniso = lik.item[["aniso"]], control.pmm = control.pmm,
      verbose = verbose
    )
    
    if( !t.Valphaxi[["error"]] ){
      lik.item[["Valphaxi"]] <- t.Valphaxi
    } else {
      return( lik.item )
    }
    
  }
  
  ##  estimate fixed and random effects (zhat, betahat, bhat, residuals )
  ##  and estimate of signal variance for Gaussian (RE)ML
  
  if( verbose > 4 ) cat(
    "\n     likelihood.calculations: computing 'zhat' object\n"    
  )
  
  ##  either take initial guess of betahat and bhat for the current
  ##  irwls iteration from initial.object or from previous iteration
  
  if( 
    !irwls.initial && !is.null( lik.item[["zhat"]][["zhat"]] )
  ){
    zhat <- lik.item[["zhat"]][["zhat"]]
  }
  
  lik.item[["zhat"]] <- estimate.zhat( 
    compute.zhat,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, 
    irwls.maxiter, irwls.ftol,
    lik.item[["param"]]["nugget"], lik.item[["eta"]], reparam,
    lik.item[["Valphaxi"]][["Valphaxi.inverse"]],
    control.pmm,
    verbose
  )
  
  if( lik.item[["zhat"]][["error"]] ) return( lik.item )     ##  an error occurred
  
  ##  compute Q matrix and its Cholesky factor (required for Gaussian
  ##  (RE)ML estimation)
  
  if( tuning.psi >= tuning.psi.nr ) {
    
    ## compute matrix Q
    
    if( verbose > 4 ) cat(
      "\n     likelihood.calculations: computing 'Q' object\n"    
    )
    
    t.Q <- f.aux.Qstar( 
      TT = TT, TtT = TtT, 
      XX = XX, col.rank.XX = col.rank.XX, min.condnum = min.condnum,
      Vi = lik.item[["Valphaxi"]][["Valphaxi.inverse"]], 
      eta = lik.item[["eta"]], 
      ml.method = ml.method, control.pmm = control.pmm
    )
    
    if( !t.Q[["error"]] ){
      lik.item[["Q"]] <- t.Q
    } else {
      return( lik.item )
    }

  }
  
  ##  store updated lik.item object
  
  assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
#   print( str( lik.item ) ); stop()
  
  return( lik.item )
    
}


##   ##############################################################################

partial.derivatives.variogram <- 
  function(
    x, variogram.model, param, d.param, aniso, verbose 
  )
{
  
  ##  Function to compute partial derivatives of generalized
  ##  correlation matrix with respect to scale and extra parameters
  
  ##  Arguments:
  ##  x             lag vectors for all pairs of distinct locations
  ##  variogram.model         Covariance Model as in Variogram{RandomFields}
  ##  param         Vector with variogram parameters
  ##  d.param       String, Parameter for which to determine the derivative
  
  ##  Value:
  ##  Vector or Matrix with partial derivative of Valphaxi for scale and extra parameters
  ##                named a, b, c, ... as in Variogram{RandomFields}
  
  ##  References:
  ##  help(Variogram)
  ##  Chiles and Delfiner, Section 2.5
  
  ##  06 Apr 2011  C.Schwierz
  ##  2011-07-17 ap
  ##  2012-01-24 ap RMcauchytbm and RMlgd models added
  ##  2012-01-25 ap extra model parameter with same names as in Variogram{RandomFields}
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ##  2013-06-12 AP substituting [["x"]] for $x in all lists
  ##  2014-05-15 AP changes for version 3 of RandomFields
  ##  2015-07-17 AP new name of function, scaling with 1-xi eliminated
  
  aniso.name <- names( aniso[["aniso"]] )
  alpha <- unname( param["scale"] )
  n = NCOL( x )
  aux <- aniso[["rotmat"]] %*% t(x)
  
  ## scaled lag distance
  
  hs <- sqrt( colSums( ( aniso[["sclmat"]] * aux )^2 ) ) / alpha
  
  ## partial derivatives of scaled lag distance with respect to
  ## anisotropy parameters
  
  dhs.daniso <- switch(
    
    d.param,
    
    f1 = {
      colSums(
        ( c( 0., -1. / aniso[["aniso"]]["f1"]^2, 0. )[1:n] * aniso[["sclmat"]] ) * aux^2 
      )
    },
    
    f2 = { 
      colSums(
        ( c( 0., 0., -1. / aniso[["aniso"]]["f2"]^2 )[1:n] * aniso[["sclmat"]] ) * aux^2 
      )
    },
    omega = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(             sp*co,            -sp*so, 0. ),
          c(  co*sz*cp + cz*so,  cz*co - sz*cp*so, 0. ),
          c( -cz*co*cp + sz*so,  co*sz + cz*cp*so, 0. )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    phi = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(     cp*so,     cp*co,    -sp ),
          c( -sz*sp*so, -co*sz*sp, -cp*sz ),
          c(  cz*sp*so,  cz*co*sp,  cz*cp )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    zeta = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(                0.,               0.,     0. ),
          c(  co*sz + cz*cp*so, cz*co*cp - sz*so, -cz*sp ),
          c( -cz*co + sz*cp*so, co*sz*cp + cz*so, -sp*sz )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    NA
  ) / ( hs * alpha^2 )
  
  ##  partial derivative of scaled lag distance with respect to scale
  ##  parameter
  
  dhs.dscale <- -hs / alpha
  
  ##  compute derivative of generalized correlation matrix with
  ##  respect to scale and extra parameters
  
  result <- switch(
    variogram.model,
    
    RMbessel = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 2^A * besselJ( hs, 1+A ) * gamma( 1+A ) ) / hs^A
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^A * besselJ( hs, 1+A ) * gamma(1 + A) ) / hs^A,
      #   0.
      # )
      
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector( 
            attr( 
              numericDeriv( 
                expr = quote( 
                  2^nu * gamma( nu+1 ) * besselJ( hs, nu ) / hs^nu 
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            ) 
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMbessel
    
    RMcauchy = {
      
      A <- unname( param["gamma"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * A * hs * ( 1+hs^2 )^(-1-A)        
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        gamma = {
          -( 1 + hs^2 )^(-A) * log( 1 + hs^2 )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcauchy
    
#     RMcauchytbm = {
#       
#       A <- unname( param["alpha"] )
#       B <- unname( param["beta"] )
#       C <- unname( param["gamma"] )
#       
#       ## derivative with of generalized covariance with respect to
#       ## scaled lag distance
#       
#       dgc.dhs <- -( 
#         B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - (-B + C) * hs^A ) 
#       ) / C
#       # dgc.dhs <- ifelse(
#       #   hs > 0.,
#       #   -( 
#       #     B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - B * hs^A + C * hs^A) 
#       #   ) / C,
#       #   if( A > 1. ){
#       #     0.
#       #   } else if( identical( A, 1. ) ){
#       #     -B * (1+C) / C
#       #   } else {
#       #     -Inf
#       #   }
#       # )
#       
#       
#       switch(
#         d.param,
#         scale = dgc.dhs * dhs.dscale,
#         # scale = {
#         #   ( B * hs^A * (1+hs^A)^(-2-B/A) * (A + C + (-B+C) * hs^A ) ) / ( C * scale )
#         # },
#         alpha = {
#           ( B * (1+hs^A)^(-2 - B/A) * (
#               -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
#               ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
#             ) 
#           ) / (A^2 * C )
#         },
#         # alpha = {
#         #   ifelse(
#         #     hs > 0.,
#         #     ( B * (1+hs^A)^(-2 - B/A) * (
#         #         -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
#         #         ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
#         #       ) 
#         #     ) / (A^2 * C ),
#         #     0.
#         #   )
#         # },
#         beta = {
#           ( -( A * hs^A) - (C + (-B+C) * hs^A ) * log( 1+hs^A ) ) / 
#           ( A*C * (1+hs^A)^( (A+B)/A ) )
#         },
#         gamma = {
#           ( B * hs^A ) / ( C^2 * (1+hs^A)^( (A+B)/A) )
#         },
#         dgc.dhs * dhs.daniso
#       )
#     }, ##  end case RMcauchytbm
    
    RMcircular = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( -4 * sqrt( 1-hs[sel]^2 ) ) / pi
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcircular
    
    RMcubic = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- hs[sel] * ( -14. + 26.25*hs[sel] - 17.5*hs[sel]^3 + 5.25*hs[sel]^5 )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcubic
    
    RMdagum = {
      
      A <- unname( param["beta"] )
      B <- unname( param["gamma"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -B / ( hs * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( B / ( hs * ( 1+ hs^(-A) )^(B/A) * (1 + hs^A ) ) ),
      #   -Inf
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     B / ( ( 1 + hs^(-A) )^(B/A) * (1 + hs^A) * scale ),
        #     0.
        #   )
        # },
        beta = {
          -( B * ( A * log(hs) + (1+hs^A) * log( 1+hs^(-A) ) ) ) /
          ( A^2 * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) )
        },
        # beta = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * ( A * log(hs) + (1+hs^A) * log( 1+hs^(-A) ) ) ) /
        #     ( A^2 * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) ),
        #     0.
        #   )
        # },
        gamma = {
          log( 1 + hs^(-A) ) / ( A * (1 + hs^(-A) )^(B/A) )
        },
        # gamma = {
        #   ifelse(
        #     hs > 0.,
        #     log( 1 + hs^(-A) ) / ( A * (1 + hs^(-A) )^(B/A) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdagum
    
    RMdampedcos = {
      
      A <- unname( param["lambda"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( A * cos(hs) + sin(hs) ) / exp( A*hs ) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        lambda = {
          -exp( -A * hs ) * hs * cos( hs )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdampedcos
    
    RMdewijsian = {
      
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -A / ( hs + hs^(1-A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A / ( hs + hs^(1-A) ) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * hs^A )/( scale + hs^A * scale )
        # },
        alpha = {
          -( ( hs^A * log( hs ) ) / ( 1 + hs^A ) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( ( hs^A * log( hs ) ) / ( 1 + hs^A ) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdewijsian
    
    
    RMexp = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -exp( -hs )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case exponential
    
    RMfbm = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -A * hs^(-1+A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * hs^(-1+A) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   A * hs^A / scale
        # },
        alpha = {
          -hs^A * log( hs )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( hs^A * log( hs ) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMfbm
    
    RMgauss = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * hs / exp( hs^2 )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgauss
    
    RMgenfbm = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["delta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( A * B * hs^(-1+A) * (1+hs^A)^(-1+B))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * B * hs^(-1+A) * (1+hs^A)^(-1+B)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * B * hs^A * (1+hs^A)^(-1+B) ) / scale
        # },
        alpha = {
          -( B * hs^A * (1+hs^A)^(-1+B) * log(hs) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * hs^A * (1+hs^A)^(-1+B) * log(hs) ),
        #     0.
        #   )
        # },
        delta = {
          -( (1 + hs^A )^B * log( 1 + hs^A ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgenfbm
    
    RMgencauchy = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( B * hs^(-1+A)) / (1+hs^A)^((A+B)/A))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( B * hs^(-1+A)) / (1+hs^A)^((A+B)/A)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( B * hs^A ) / ( (  1+hs^A)^((A+B)/A) * scale )
        # },
        alpha = {
          B * ( 1 + hs^A )^(-(A+B)/A) * (
            -A * hs^A * log( hs ) +
            ( 1 + hs^A ) * log( 1 + hs^A )
          ) / A^2
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     B * ( 1 + hs^A )^(-(A+B)/A) * (
        #       -A * hs^A * log( hs ) +
        #       ( 1 + hs^A ) * log( 1 + hs^A )
        #     ) / A^2,
        #     0.
        #   )
        # },
        beta = {
          -( log( 1+hs^A ) / ( A * (1+hs^A)^(B/A) ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgencauchy
    
    #     gengneiting = { Version 2 of package RandomFields
    #       
    #       
    #       A <- unname( param["n"] )
    #       B <- unname( param["alpha"] )
    #       
    #       ## derivative with of generalized covariance with respect to
    #       ## scaled lag distance
    #       
    #       dgc.dhs <- rep( 0., length( hs ) )
    #       sel <- hs < 1.
    #       dgc.dhs[sel] <- if( identical( A, 1 ) ){
    #         -( (1+B) * (2+B) * (1-hs[sel])^B * hs[sel] )
    #       } else if( identical( A, 2 ) ){
    #         -( (3+B) * (4+B) * (1-hs[sel])^(1+B) * hs[sel] * ( 1 + hs[sel] + B*hs[sel]) ) / 3.
    #       } else if( identical( A, 3 ) ){
    #         -( 
    #           (5+B) * (6+B) * (1-hs[sel])^(2+B) * hs[sel] * ( 3 + 3 * (2+B) * hs[sel] + (1+B) * (3+B) * hs[sel]^2 ) 
    #         ) / 15.
    #       } else {
    #         stop( "gengneiting model undefined for 'n' != 1:3" )
    #       }
    #       
    #       result <- rep( 0., length( hs ) )
    #       
    #       switch(
    #         d.param,
    #         scale = dgc.dhs * dhs.dscale,
    #         alpha = {
    #           result[sel] <- if( identical( A, 1 ) ){
    #             (1-hs[sel])^(1+B) * ( hs[sel] + (1 + hs[sel] + B*hs[sel]) * log( 1-hs[sel]) )
    #             
    #           } else if( identical( A, 2 ) ){
    #             (
    #               (1-hs[sel])^(2+B) * ( 
    #                 hs[sel] * ( 3 + 2 * (2+B) *hs[sel] ) + 
    #                 ( 3 + 3 * ( 2+B) * hs[sel] + ( 1+B) * (3+B) * hs[sel]^2 ) * log( 1-hs[sel] )
    #               )
    #             ) / 3.
    #           } else if( identical( A, 3 ) ){
    #             ( 
    #               (1-hs[sel])^(3+B) * ( 
    #                 hs[sel] * ( 15 + hs[sel] * ( 36 + 23*hs[sel] + 3 * B * ( 4 + (6+B)*hs[sel] ) ) ) + 
    #                 ( 15 + 15 * (3+B) * hs[sel] + ( 45 + 6 * B * (6+B) ) * hs[sel]^2 + (1+B) * (3+B) * (5+B) * hs[sel]^3 ) * 
    #                 log( 1-hs[sel]) 
    #               ) 
    #             ) / 15.
    #           }
    #           result
    #         },
    #         dgc.dhs * dhs.daniso
    #       )
    #       
    #       
    #     }, ##  end case Gengneiting
    
    RMgengneiting = {
      
      
      A <- unname( param["kappa"] )
      B <- unname( param["mu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- if( identical( A, 1 ) ){
        (2.5+B) * (1-hs[sel])^(2.5+B) - (2.5+B) * (1-hs[sel])^(1.5+B) * (1+(2.5+B) * hs[sel])
      } else if( identical( A, 2 ) ){
        (1 - hs[sel])^(4.5+B) * (4.5+B + 2/3 * (3.5+B) * (5.5+B) * hs[sel] ) - 
        (4.5+B) * (1 - hs[sel])^(3.5+B) * (
          1 + hs[sel]*(4.5 + B + 6.416666666666666*hs[sel] + B/3. * (9.+B) * hs[sel] ) 
        )
      } else if( identical( A, 3 ) ){
        (1 - hs[sel])^(6.5+B) * (6.5 + B + 0.8 * (5.275255128608411+B) * (7.724744871391589+B) * hs[sel] + 
          0.2 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^2) - 
        (6.5+B) * (1 - hs[sel])^(5.5+B) * (1 + (6.5+B) * hs[sel] + 0.4 * (5.275255128608411+B) * 
          (7.724744871391589+B) * hs[sel]^2 + 0.2/3 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^3
        )
      } else {
        stop( "RMgengneiting model undefined for 'n' != 1:3" )
      }
      
      result <- rep( 0., length( hs ) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        mu = {
          result[sel] <- if( identical( A, 1 ) ){
            (1 - hs[sel])^(2.5+B) * (hs[sel] + (1 + (2.5+B) * hs[sel]) * log(1 - hs[sel]))
          } else if( identical( A, 2 ) ){
            (1 - hs[sel])^(4.5+B) * (hs[sel] + 2/3 * (4.5+B) * hs[sel]^2 + (1 + hs[sel] * (
                  4.5 + B + 6.416666666666666*hs[sel] +  B/3. * (9.+B) * hs[sel]) ) * log(1 - hs[sel])
            )
          } else if( identical( A, 3 ) ){
            (1 - hs[sel])^(6.5 + B)*
            (hs[sel] + (5.2 + 0.8*B)*hs[sel]^2 + 
              0.2*(5.345299461620754 + B)*
              (7.654700538379246 + B)*hs[sel]^3 + 
              (1 + hs[sel]*(6.5 + 1.*B + 
                  0.4*(5.275255128608411 + B)*
                  (7.724744871391589 + B)*hs[sel] + 
                  0.06666666666666667*(4.5 + B)*

                  (6.5 + B)*(8.5 + B)*hs[sel]^2))*
              log(1 - hs[sel]))          
          }
          result
        },
        dgc.dhs * dhs.daniso
      )
      
      
    }, ##  end case Gengneiting
    

    RMgneiting = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < (1. / 0.301187465825)
      dgc.dhs[sel] <- (1. - 0.301187465825*hs[sel])^7 * (
        -1.9957055705418814*hs[sel] -  4.207570523270417*hs[sel]^2 - 2.896611435848653*hs[sel]^3
      )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgneiting
    
    RMlgd = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( A * B * hs^(-1-B) ) / (A+B)
      sel <- hs <= 1.
      dgc.dhs[sel] <- -( A * B * hs[sel]^(-1+A) ) / (A+B)
      
      # dgc.dhs <- ifelse(
      #   hs > 0.
      #   ifelse(
      #     hs <= 1.,
      #     -( A * B * hs^(-1+A) ) / (A+B),
      #     -( A * B * hs^(-1-B) ) / (A+B)
      #   ),
      #   if( identical( A, 1. ) ){
      #     -B / ( B + 1 )
      #   } else if( A < 1. ){
      #     -Inf
      #   } 
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse( 
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       A * B * hs^A,
        #       A * B / hs^B
        #     ) / ( (A+B) * scale ),
        #     0.
        #   )
        # },
        alpha = {
          result <- B / ( (A+B)^2 * hs^B )
          sel <- hs <= 1.
          result[sel] <- -( B * hs[sel]^A * ( -1 + (A+B ) * log( hs[sel] ) ) ) / (A+B)^2
          result
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       -( B * hs^A * ( -1 + (A+B ) * log( hs ) ) ) / (A+B)^2,
        #       B / ( (A+B)^2 * hs^B )              
        #     ),
        #     0.
        #   )
        # },
        beta = {
          result <- -A * ( 1 + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
          sel <- hs <= 1.
          result[sel] <- -A * hs[sel]^A / (A+B)^2
          result
        },
        # beta = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       -A * hs^A / (A+B)^2,
        #       -A * ( 1 + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
        #     ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMlgd
    
    RMmatern = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 
        2^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2*A)*hs, -1+A )
      ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -(
      #     ( 2^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2) * sqrt(A) * hs , -1+A )
      #     ) / gamma(A)
      #   ),
      #   if( A < 0.5 ){
      #     -Inf
      #   } else if( identical( A, 0.5 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse( 
        #     hs > 0.,
        #     ( 2^(1.5 - A/2.) * ( sqrt(A) * hs )^(1+A) * 
        #       besselK( sqrt(2*A) * hs, A-1) 
        #     ) / (scale * gamma(A) ),
        #     0.
        #   )
        # },
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector( 
            attr( 
              numericDeriv( 
                expr = quote( 
                  2^(1.-nu) / gamma(nu) * 
                  ( sqrt( 2*nu ) * hs )^nu * besselK( sqrt( 2*nu ) * hs, nu )
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            ) 
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMmatern
    
    RMpenta = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( 11 * (-1+hs[sel])^5 * hs[sel] * (2+hs[sel]) * ( 4 + hs[sel] * ( 18 + 5 * hs[sel] * (3+hs[sel]) ) ) ) / 6.
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMpenta
    
    RMaskey = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -(A * (1-hs[sel])^(-1+A))         
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          result <- rep( 0., length( hs ) )
          sel <- hs < 1.
          result[sel] <- ( 1 - hs[sel] )^A * log( 1 - hs[sel] )
          result
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMaskey
    
    RMqexp = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- 2 * (-A + exp(hs) ) / ( (-2+A ) * exp(2*hs) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          ( 2 * exp( -2*hs ) * ( -1 + exp( hs ) ) ) / (-2+A)^2
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMqexp
    
    
    RMspheric = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -1.5 + 1.5 * hs[sel]^2          
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMspheric
    
    RMstable = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( A * hs^(-1+A) ) / exp(hs^A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( A * hs^(-1+A) ) / exp(hs^A) ),
      #   if( A > 1. ){
      #     0.
      #   } else if( identical( A, 1. ) ){
      #     -1.
      #   } else {
      #     -Inf            
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * exp( -hs^A ) * hs^A ) / scale
        # },
        alpha = {
          -exp( -hs^A ) * hs^A * log( hs )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -exp( -hs^A ) * hs^A * log( hs ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMstable
    
    RMwave = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- ( hs * cos(hs) - sin(hs) ) / hs^2
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   ( hs * cos(hs) - sin(hs) ) / hs^2,
      #   0.
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case wave
    
    RMwhittle = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 2^(1-A) * hs^A * besselK( hs, -1+A ) ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^(1-A) * hs^A * besselK( hs, -1+A ) ) / gamma(A),
      #   if( A < 0.5 ){
      #     -Inf
      #   } else if( identical( A, 0.5 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     ( 
        #       2^(1-A) * h * hs^A * besselK( hs, -1+A ) 
        #     ) / ( scale^2 * gamma(A) ),
        #     0.
        #   )
        #   
        # },
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector( 
            attr( 
              numericDeriv( 
                expr = quote( 
                  2^(1.-nu) / gamma(nu) * hs^nu * besselK( hs, nu )
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            ) 
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case whittle
    
    stop(
      paste( 
        variogram.model, 
        "model: derivatives for scale, extra variogram and anisotropy parameters undefined" 
      )
    )
    
  ) ##  end switch cov.model
    
  ##  convert to matrix
  
  result <- list(
    diag = rep( 0., 0.5 * ( 1 + sqrt( 1 + 8 * length( result ) ) ) ),
    tri = result
  )
  attr( result, "struc" ) <- "sym"
  result <- expand( result )
  
  return( result )
  
}

##   ##############################################################################

estimating.equations.theta <- 
  function(
    adjustable.param,
    envir,
    variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, 
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    force.gradient,
    expectations,
    error.family.estimation,
    control.pmm, 
    verbose
  )
{
  
  ## function evaluates the robustified estimating equations of
  ## variogram parameters derived from the Gaussian log-likelihood
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP TtT passed to function
  ## 2015-07-17 AP new function interface, improved efficiency
  ## 2015-07-27 AP changes to further improve efficiency
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2015-08-19 AP control about error families for computing covariances added  

  ##  get lik.item
  
  lik.item <- likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param, reparam = FALSE,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE, 
    control.pmm = control.pmm,
    verbose = verbose
  )
  
#   print( str( lik.item ) ); stop()
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    t.result <- rep( Inf, length( adjustable.param ) )
    names( t.result ) <- names( adjustable.param )
    return( t.result )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    t.result <- rep( Inf, length( adjustable.param ) )
    names( t.result ) <- names( adjustable.param )
    return( t.result )
  }
  
  ##  check whether estimating equations should be computed for fixed parameters
  
  if( length( adjustable.param ) == 0 && force.gradient ){
    adjustable.param <- fixed.param
  }
  
  ##  evaluate estimating equations
  
  if( length( adjustable.param ) > 0 ){
    
    ##  compute Cov[bhat]
    
    r.cov <- covariances.fixed.random.effects(
      Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
      Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
      Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
      Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
      rweights = lik.item[["zhat"]][["rweights"]],
      XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
      nugget = lik.item[["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations, family = error.family.estimation,
      cov.bhat = TRUE, full.cov.bhat = TRUE,
      cov.betahat = FALSE,
      cov.bhat.betahat = FALSE,
      cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
      cov.delta.bhat.betahat = FALSE,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = FALSE,
      control.pmm = control.pmm,
      verbose = verbose
    )
    
    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0 ) cat(
        "\nan error occurred when computing the covariances of fixed and random effects\n"
      )
      t.result <- rep( Inf, length( adjustable.param ) )
      names( t.result ) <- names( adjustable.param )
      return( t.result )
    }
    
    ## compute further required items
    
    ## Valphaxi^-1 Cov[bhat]
        
    Valphaxii.cov.bhat <- pmm( 
      lik.item[["Valphaxi"]][["Valphaxi.inverse"]], r.cov[["cov.bhat"]], 
      control.pmm
    )
    
    ##  computation of estimating equations for all elements of adjustable.param
    
    t.eeq <- sapply(
      names( adjustable.param ),
      f.aux.eeq,
      param = lik.item[["param"]], aniso = lik.item[["aniso"]],
      xi = lik.item[["xi"]],
      Valphaxii = lik.item[["Valphaxi"]][["Valphaxi.inverse"]], 
      Valphaxii.cov.bhat = Valphaxii.cov.bhat,
      bh = lik.item[["zhat"]][["bhat"]],
      bhVaxi = lik.item[["zhat"]][["Valphaxi.inverse.bhat"]],
      r.cov = r.cov, lik.item = lik.item,
      TtT = TtT, 
      lag.vectors = lag.vectors, variogram.model = variogram.model, 
      control.pmm = control.pmm, verbose = verbose
    )
    
    eeq.exp <- t.eeq["eeq.exp", ]
    eeq.emp <- t.eeq["eeq.emp", ]
    
    names( eeq.exp ) <- names( adjustable.param )
    names( eeq.emp ) <- names( adjustable.param )
    
    if( verbose > 1 ) {
      cat( "\n                      ",
        format( names( eeq.emp), width = 14, justify = "right" ), 
        "\n", sep =""
      )
      cat( "  EEQ                :", 
        format( 
          signif( eeq.emp / eeq.exp - 1, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n", sep = "" 
      )
      if( verbose > 2 ){
        cat( "      empirical terms:", 
          format( 
            signif( eeq.emp, digits = 7 ), 
            scientific = TRUE, width = 14
          ), "\n", sep = "" 
        )
        cat( "      expected  terms:", 
          format( 
            signif( eeq.exp, digits = 7 ), 
            scientific = TRUE, width = 14
          ), "\n", sep = ""
        )
      }
      cat("\n")
    }
    
    ##  store terms in lik.item object
    
    lik.item[["eeq"]] <- list(
      eeq.emp = eeq.emp,
      eeq.exp = eeq.exp
    )
    
    assign( "lik.item", lik.item, pos = as.environment( envir ) )
    
    return( eeq.emp / eeq.exp - 1. )
        
  } else {
    
    ##  all parameters are fixed
    
    return( NA_real_ )
    
  }
  
}


##   ##############################################################################

negative.loglikelihood <- 
  function(
    adjustable.param,
    envir,
    variogram.model, fixed.param, 
    param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, 
    tuning.psi, tuning.psi.nr, ml.method, reparam,
    irwls.initial, irwls.maxiter, irwls.ftol,
    control.pmm, 
    verbose,
    ...
  )
{
  
  ## function computes to negative (un)restricted loglikelihood
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-03 AP changes for estimating zhat
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-27 AP correcting error in likelihood for original parametrization (reparam == FALSE)
  #     sel <- !c( param.name, aniso.name ) %in% names( fixed.param )
  #     names( adjustable.param ) <- c( param.name, aniso.name )[sel]
  
  ##  compute required items (param, eta, Valphaxi.inverse, Valphaxi.ilcf, 
  ##  betahat, bhat, residuals, etc.)
  
  lik.item <- likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE, 
    control.pmm = control.pmm,
    verbose
  )
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    return( NA )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( NA )
  }
  
  ##  check whether Q matrix not positive definite
  
  if( lik.item[["Q"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\nan error occurred when determinants required for",
      "Gaussian log-likelihood were computed\n"
    )
    return( NA )
  }
  
  ##  compute negative (restricted || profile) loglikelihood 
  
  t.q <- t.qmp <- NROW(XX)
  if( identical( ml.method, "REML" ) ){
    t.qmp <- t.q - col.rank.XX[["rank"]]  
  }
  
  if( reparam ){
    
    ## (restricted) profile loglikelihood
    
    r.neg.loglik <- 0.5 * (
      t.qmp * ( 
        1. + log(2*pi) - log(t.qmp ) + 
        log( lik.item[["zhat"]][["RSS"]] )
      ) - 
      t.q * log( lik.item[["eta"]] ) + 
      lik.item[["Q"]][["log.det.Qstar"]] +
      lik.item[["Valphaxi"]][["log.det.Valphaxi"]]
    )
    
  } else {
  
    ## (restricted) loglikelihood
    
    t.param <- f.reparam.fwd( lik.item[["param"]] )
    
    r.neg.loglik <- 0.5 * (
      t.qmp * ( 
        log(2*pi) + log( t.param["variance"] )
      ) - 
      t.q * log( lik.item[["eta"]] ) + 
      lik.item[["Valphaxi"]][["log.det.Valphaxi"]] +
      lik.item[["Q"]][["log.det.Qstar"]] +
      lik.item[["zhat"]][["RSS"]] / t.param["variance"]
    )
    
  }
  
  attributes( r.neg.loglik ) <- NULL
  
  if( verbose > 1 ) cat(
    "\n  Negative. restrict. loglikelihood:", 
    format( 
      signif( r.neg.loglik, digits = 7 ), 
      scientific = TRUE, width = 14
    ), "\n", sep = ""
  )
  
  return( r.neg.loglik )
  
}


##   ##############################################################################

gradient.negative.loglikelihood <- 
  function(
    adjustable.param,
    envir,
    variogram.model, fixed.param,
    param.name, aniso.name,
    param.tf, deriv.fwd.tf, bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, 
    tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.ftol,
    force.gradient,
    control.pmm,
    verbose
  )
{
  
  ##  function computes gradient of Laplace approximation of negative 
  ##  restricted log-likelihood with respect to covariance parameters
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP correction of values returned on error
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2015-12-02 AP reparametrized variogram parameters renamed
  ##  get lik.item
  
  lik.item <- likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat, 
    psi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pmm = control.pmm,
    verbose
  )
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether Q matrix not positive definite
  
  if( lik.item[["Q"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0 ) cat(
      "\nan error occurred when determinants required for ",
      "Gaussian log-likelihood were computed\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether gradient should be computed for fixed parameters
  
  if( length( adjustable.param ) == 0 && force.gradient ){
    adjustable.param <- fixed.param
  }
  
  ##  evaluate gradient
  
  if( length( adjustable.param ) > 0 ){
    
    ##  compute auxiliary items
    
    n <- nrow( XX )
    
    if( reparam ){
      
      Qsi       <- lik.item[["Q"]][["Qstar.inverse"]]
      Valphaxii <- Valphaxii <- lik.item[["Valphaxi"]][["Valphaxi.inverse"]]
      bh        <- lik.item[["zhat"]][["bhat"]]
      bhVaxi    <- lik.item[["zhat"]][["Valphaxi.inverse.bhat"]]
      if( !identical( names( adjustable.param ), "nugget" ) ){
        Qst11Vai  <- pmm( Qsi[1:n, 1:n], Valphaxii, control.pmm )
      } else {
        Qst11Vai  <- NULL
      }
      
    } else {
      
      t.param <- f.reparam.fwd( lik.item[["param"]] )
      
      Qi     <- lik.item[["Q"]][["Qstar.inverse"]] * t.param["variance"]
      Vi     <- lik.item[["Valphaxi"]][["Valphaxi.inverse"]] / t.param["variance"]
      bh     <- lik.item[["zhat"]][["bhat"]]
      bhVi   <- lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] / t.param["variance"]
      if( !identical( names( adjustable.param ), "nugget" ) ){
        Qt11Vi <- pmm( Qi[1:n, 1:n], Vi, control.pmm )
      } else {
        Qt11Vi <- NULL
      }

    }
    
    ##  compute gradient for elements of adjustable.param
          
    if( reparam ){
      
      ## gradient for reparametrized variogram parameters
      
      r.gradient <- sapply(
        names( adjustable.param ),
        f.aux.gradient.npll,
        param = lik.item[["param"]], aniso = lik.item[["aniso"]],
        eta = lik.item[["eta"]], xi = lik.item[["xi"]],
        TT = TT, TtT = TtT, XX = XX, col.rank.XX = col.rank.XX,
        res = lik.item[["zhat"]][["residuals"]],
        Ustar = lik.item[["zhat"]][["RSS"]],
        Qsi = Qsi, Qst11Vai = Qst11Vai, 
        Valphaxii = Valphaxii,
        bh = bh, bhVaxi = bhVaxi,
        lag.vectors = lag.vectors, variogram.model = variogram.model,
        param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf, 
        ml.method = ml.method, 
        control.pmm = control.pmm, verbose = verbose
      )
      
    } else {
      
      ## gradient for reparametrized variogram parameters
            
      r.gradient <- sapply(
        names( adjustable.param ),
        f.aux.gradient.nll,
        param = lik.item[["param"]], aniso = lik.item[["aniso"]],
        TT = TT, TtT = TtT, XX = XX, 
        res = lik.item[["zhat"]][["residuals"]],
        Qi = Qi, Vi = Vi, Qt11Vi = Qt11Vi, 
        bh = bh, bhVi = bhVi,
        lag.vectors = lag.vectors, variogram.model = variogram.model,
        param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf, 
        ml.method = ml.method, 
        control.pmm = control.pmm, verbose = verbose
      )
      
    }
        
    names( r.gradient ) <- names( adjustable.param )
        
    ##  rearrange elements of gradient and change sign (for negative
    ##  log-likelihood)
    
    r.gradient <- -r.gradient[names( adjustable.param )]
    
    if( verbose > 1 ){
      
      t.gradient <- r.gradient
      
      if( reparam ){
        tmp <- names( t.gradient )
        tmp[tmp %in% "snugget"] <- "xi"
        tmp[tmp %in% "nugget"]  <- "eta"
        names( t.gradient ) <- tmp  
      }
      
      cat( "\n                      ",
        format( names( t.gradient ), width = 14, justify = "right" ), 
        "\n", sep = ""
      )
      cat( "  Gradient           :", 
        format( 
          signif( t.gradient, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n" , sep = ""
      )
    }
    
    return( r.gradient )
    
  } else {
    
    ##  all parameters are fixed
    
    return( NA_real_ )
    
  }
}


##  ##   ##############################################################################
##      
##      f.compute.df <- function( Valphaxi, XX, param ){
##          
##          ##  function computes three estimates of the degrees of freedom of
##          ##  the smoothing universal kriging predictor, cf.  Hastie &
##          ##  Tibshirani, 1990, Generalized additive models, pp.52
##          
##          ##  2011-07-05
##          ##  Andreas Papritz
##          
##          sigma <- param["variance"] * Valphaxi
##          diag( sigma ) <- diag( sigma ) + param["nugget"]
##          
##          ##  compute inverse lower cholesky factor of covariance matrix of
##          ##  data
##          
##          ilcf <- t( backsolve( chol( sigma ), diag( nrow( Valphaxi ) ), k = nrow( Valphaxi ) ) )
##          
##          ##  compute hat matrix
##          
##          q <- qr.Q( qr( xtilde <- ilcf %*% XX ) )
##          s <- -tcrossprod( q )
##          
##          diag( s ) <- diag( s ) + 1
##          s <- -param["nugget"] * t( ilcf ) %*% s %*% ilcf
##          diag( s ) <- diag( s ) + 1
##          
##          ##  compute degrees of freedom
##          
##          df.1 <- sum( diag( s ) )
##          df.3 <- sum( s^2 )
##          df.2 <- 2 * df.1 - df.3
##          
##          return( 
##              c( 
##                  df.SSt    = t.df.2 <- sum( s^2 ), 
##                  df.S      = t.df.1 <- sum( diag( s ) ), 
##                  df.2SmSSt = 2 * t.df.1 - t.df.2
##              ) 
##          )
##              
##      }

##   ##############################################################################

georob.fit <- 
  function(
    ## root.finding,
    initial.objects,
    variogram.model, param, fit.param,
    aniso, fit.aniso,
    param.tf, 
    fwd.tf, 
    deriv.fwd.tf, 
    bwd.tf,
    georob.object,
    safe.param,
    tuning.psi, 
    error.family.estimation, error.family.cov.effects, error.family.cov.residuals,
    cov.bhat, full.cov.bhat,
    cov.betahat, 
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    min.condnum, col.rank.XX,
    psi.func,
    tuning.psi.nr,
    ml.method,
    maximizer,
    reparam,
    irwls.initial,
    irwls.maxiter, 
    irwls.ftol, 
    force.gradient,
    zero.dist,
    control.nleqslv,
    control.optim, 
    control.nlminb,
    hessian,
    control.pmm,
    verbose
  )
{
  
  ## 2011-06-24 ap
  ## 2011-06-24 cs
  ## 2011-06-29 ap, cs
  ## 2011-07-22 ap
  ## 2011-07-28 ap
  ## 2011-08-12 ap
  ## 2011-10-14 ap
  ## 2011-12-19 ap
  ## 2011-12-22 ap
  ## 2011-12-23 AP modified for estimating variogram model with spatial
  ##               nugget (micro-scale variation)
  ## 2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2012-02-20 AP replacement of ifelse
  ## 2012-02-27 AP rescaled rho-, psi-function etc.
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-21 AP arguments lower, upper passed to optim
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2012-11-27 AP changes in check allowed parameter range
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP changes in stored items of Valphaxi object
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-03 AP new transformation of rotation angles
  ## 2013-07-09 AP catching errors occuring when fitting anisotropic
  ##               variograms with default anisotropy parameters
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-02-18 AP correcting error when fitting models with offset
  ## 2014-05-28 AP change in check for initial variogram parameter values
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-10 AP changes for reparametrization of variogram  
  ## 2015-03-16 AP elimination of unused variables, own function for psi function
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram, 
  ##               nlminb added as maximizer of loglikelihood
  ## 2015-07-23 AP changes for avoiding computation of Valphaxi object if not needed, 
  ##               rearrangement of output item (Valpha.objects, zhat.objects)
  ## 2015-08-19 AP variances of eps and psi(eps/sigma) for long-tailed error distribution; 
  ##               computing covariances of residuals under long-tailed error model,
  ##               control about error families for computing covariances added
  ## 2015-12-02 AP catching error in computation of covariances
  ## 2016-01-26 AP refined check of initial values of variogram parameters
  
  ##  ToDos:
  
  ##  main body of georob.fit
  
  d2r <- pi / 180.
  
  ##  define rho-function and derivatives (suppress temporarily warnings issued by gamma())
  
  old.op <- options( warn = -1 )
  rho.psi.etc <- f.psi.function( x = psi.func, tp = tuning.psi )
  options( old.op )
  
  ##  set number of IRWLS iterations for estimating bhat and betahat to
  ##  1 for non-robust REML case
  
  if( tuning.psi >= tuning.psi.nr ){
    irwls.maxiter <- 1        
  }
  
  ##  copy items of initial.objects to local environment
  
  XX          <- initial.objects[["x"]]
  yy          <- initial.objects[["y"]]
  betahat     <- coefficients( initial.objects[["initial.fit"]] )
  bhat        <- initial.objects[["bhat"]]
  coordinates <- initial.objects[["locations.objects"]][["coordinates"]]
  
  ##  check for multiple observations at same location and generate
  ##  designmatrix of replicated observations
  
  dist0 <- as.matrix( dist( coordinates ) ) <= zero.dist
  first.dist0 <- unname( apply( dist0, 1, function( x ) ( (1:length(x))[x])[1] ) )
  
  
  TT <- matrix( 0, nrow = length( yy ), ncol = length( yy ) )
  TT[ cbind( 1:nrow(TT), first.dist0 ) ] <- 1
  rep.obs <- (1:ncol(TT))[ apply( TT, 2, function( x ) all( x == 0 ) ) ]
  if( length( rep.obs ) > 0 )  TT <- TT[, -rep.obs]
  
  ## check whether explanatory variables are the identical for the replicated
  ## observations and issue an error if not
  
  apply( 
    TT, 
    2, 
    function( i, XX ){
      XX <- XX[as.logical(i), , drop = FALSE]
      apply( 
        XX, 
        2, 
        function( x ){
          if( length(x) > 1 && any( x[-1] != x[1] ) ) stop(
            "explanatory variables differ for some replicated observations" 
          )
        }
      )    
    },
    XX = XX
  )
  
  ## store row indices of replicated observations only
  
  TT <- drop( TT %*% 1:ncol( TT ) )
  TtT <- as.vector( table( TT ) )
  
  ##  omit elements corresponding to replicated observations in XX, bhat
  ##  and coordinates
  
  if( length( rep.obs ) > 0 ) {
    XX          <- XX[ -rep.obs, , drop = FALSE]
    bhat      <- bhat[ -rep.obs ]
    coordinates <- coordinates[ -rep.obs, , drop = FALSE]
    if( verbose > 0 ) cat( "\n", length(rep.obs), "replicated observations at", 
      length( unique( TT[rep.obs] ) ), "sampling locations\n" 
    )
  }
  
  ## compute lag vectors for all pairs of coordinates (or restore from
  ## georob.object if available and coordinates are the same)
  
  if( 
    !is.null( georob.object ) && 
    isTRUE( all.equal( georob.object[["locations.objects"]][["coordinates"]], coordinates ) )
  ){
    lag.vectors <- georob.object[["locations.objects"]][["lag.vectors"]]
  } else {
    indices.pairs <- combn( NROW( coordinates ), 2 )
    lag.vectors <- coordinates[ indices.pairs[2,], ] - coordinates[ indices.pairs[1,], ]
  }
  
  ## set snugget to zero if snugget has not been specified or if there are
  ## no replicated observations
  
  if( !"snugget" %in% names( param ) | sum( duplicated( TT ) ) == 0 ){
    param["snugget"] <- 0.
    fit.param["snugget"] <- FALSE
  }
  
  ##  check whether fitting of chosen variogram model is implemented and
  ##  return names of extra parameters (if any)
  
  ep <- param.names( model = variogram.model )
  
  ## check names of initial variogram parameters and flags for fitting
  
  param.name <- c( "variance", "snugget", "nugget", "scale", ep )
  
  if( !all( names( param ) %in% param.name ) ) stop( 
    "error in names of initial values of variogram parameters" 
  )
  
  if( !all( param.name  %in% names( param ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    paste( param.name[ !param.name %in% names( param ) ], collapse= ", "), "'"
  )
  
  if( !all( names( fit.param ) %in% param.name ) ) stop( 
    "error in names of control flags for fitting variogram parameters" 
  )
  
  if( length( param ) != length( fit.param ) || 
    !all( names( fit.param ) %in% names( param ) )
  ) stop( 
    "names of variogram parameters and control flags for fitting do not match" 
  )
  
  if( !all( is.numeric( param ) ) ) stop(
    "initial values of variogram parameters must be of mode 'numeric'"
  )
  if( !all( is.logical( fit.param ) ) ) stop(
    "fitting control flags of variogram parameters must be of mode 'logical'"
  )
  
  ##  rearrange initial variogram parameters
  
  param <- param[param.name]
  
  ## check whether intitial values of variogram parameters are valid
  
  if( param["variance"] < 0. ) stop("initial value of 'variance' must be > 0" )
  if( param["snugget"] < 0. )  stop("initial value of 'snugget' must be > 0" )
  if( param["nugget"] <= 0. ) stop("initial value of 'nugget' must be >= 0" )
  if( param["scale"] <= 0. ) stop("initial value of 'scale' must be >= 0" )
  
  param.bounds <- param.bounds( variogram.model, NCOL( coordinates ) )
  ep.param <- param[ep]
  
  if( !is.null( param.bounds ) ) t.bla <- sapply(
    1:length( ep.param ),
    function( i, param, bounds ){
      if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) stop(
        "initial value of parameter '", names( param[i] ), "' outside of allowed range" 
      )
    }, 
    param = ep.param,
    bounds = param.bounds
  )
  
  ##  rearrange and check flags controlling variogram parameter fitting 
  
  fit.param <- fit.param[param.name]
  
  if( 
    variogram.model %in% (t.models <- c( "RMfbm" ) ) && 
    ( 
      sum( duplicated( TT ) > 0 ) && all( 
        fit.param[c( "variance", "snugget", "scale" ) ] 
      ) ||
      sum( duplicated( TT ) == 0 ) && all( 
        fit.param[c( "variance", "scale" ) ] 
      ) 
    )
  ) stop( 
    "'variance', 'scale' (and 'snugget') cannot be fitted simultaneously for variograms ",
    paste( t.models, collapse = " or "), "; \n  'scale' parameter must be fixed"
  )
  
  ## re-parametrize variogram parameters for Gaussian (RE)ML estimation
  
  t.param <- f.reparam.fwd( param )
  
  reparam <- reparam && tuning.psi >= tuning.psi.nr && 
    sum( fit.param[c("variance", "snugget", "nugget")] ) > 1
    
  if( reparam ){
    param      <- t.param
    param      <- param[!names(param) %in% "variance"]
    param.name <- param.name[!param.name %in% "variance"]
    fit.param  <- fit.param[!names(fit.param) %in% "variance"]
  }
  
  ##  preparation for variogram parameter transformations
  
  all.param.tf <- param.tf
  
  t.sel <- match( param.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some variogram parameters" )
  } else {
    param.tf <- all.param.tf[t.sel]
  }
  param.tf <- sapply(
    param.tf,
    function( x ) if( length(x) > 1L ) x[variogram.model] else x
  )
  names( param.tf ) <- param.name 
  
  ##  transform initial variogram parameters
  
  transformed.param <- sapply(
    param.name,
    function( x, fwd.tf, param.tf, param ) fwd.tf[[param.tf[x]]]( param[x] ),
    fwd.tf = fwd.tf,
    param.tf = param.tf,
    param = param
  )
  
  names( transformed.param ) <- param.name 
  
  ## check names of initial anisotropy parameters and flags for fitting
  
  aniso.name <- c( "f1", "f2", "omega", "phi", "zeta" )
  
  if( !all( names( aniso ) %in% aniso.name ) ) stop( 
    "error in names of initial values of anisotropy parameters" 
  )
  
  if( !all( aniso.name  %in% names( aniso ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    aniso.name[ !aniso.name %in% names( aniso ) ], "'"
  )
  
  if( !all( names( fit.aniso ) %in% aniso.name ) ) stop( 
    "error in names of control flags for fitting  anisotropy parameters"
  )
  
  if( length( aniso ) != length( fit.aniso ) || 
    !all( names( fit.aniso ) %in% names( aniso ) )
  ) stop( 
    "names of anisotropy parameters and control flags for fitting do not match" 
  )
  
  if( !all( is.numeric( aniso ) ) ) stop(
    "initial values of anisotropy parameters must be of mode 'numeric'"
  )
  if( !all( is.logical( fit.aniso ) ) ) stop(
    "fitting control flags of anisotropy parameters must be of mode 'logical'"
  )
  
  ##  rearrange initial anisotropy parameters
  
  aniso <- aniso[aniso.name]
  
  ## check whether intitial values of anisotropy parameters are valid
  
  if( aniso["f1"] < 0. ||  aniso["f1"] > 1. ) stop(
    "initial value of parameter 'f1' must be in [0, 1]" 
  )
  if( aniso["f2"] < 0. ||  aniso["f1"] > 1. ) stop(
    "initial value of parameter 'f2' must be in [0, 1]" 
  )
  if( aniso["omega"] < 0. ||  aniso["omega"] > 180. ) stop(
    "initial value of parameter 'omega' must be in [0, 180]" 
  )
  if( aniso["phi"] < 0. ||  aniso["phi"] > 180. ) stop(
    "initial value of parameter 'phi' must be in [0, 180]" 
  )
  if( aniso["zeta"] < -90. ||  aniso["zeta"] > 90. ) stop(
    "initial value of parameter 'zeta' must be in [-90, 90]" 
  )
  
  ## adjust default initial values of anisotropy parameters if these are
  ## fitted
  
  if( fit.aniso["omega"] && identical( aniso["f1"], 1. ) ) aniso["f1"] <- aniso["f1"] - sqrt( .Machine$double.eps )
  if( fit.aniso["phi"] ){
    if( identical( aniso["f1"], 1. ) ) aniso["f1"] <- aniso["f1"] - 0.0001
    if( identical( aniso["f2"], 1. ) ) aniso["f2"] <- aniso["f2"] - 0.0001
  }
  if( fit.aniso["zeta"] && identical( aniso["f2"], 1. ) ) aniso["f2"] <- aniso["f2"] - 0.0001
  
  ##  rearrange and check flags controlling anisotropy parameter fitting 
  
  fit.aniso <- fit.aniso[aniso.name]
  
  ##  preparation for anisotropy parameter transformations
  
  t.sel <- match( aniso.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some anisotropy parameters" )
  } else {
    aniso.tf <- all.param.tf[t.sel]
  }  
  aniso.tf <- sapply(
    aniso.tf,
    function( x ) if( length(x) > 1L ) x[variogram.model] else x
  )  
  names( aniso.tf ) <- aniso.name 
  
  ##  convert angles to radian
  
  aniso[c("omega", "phi", "zeta" )] <- aniso[c("omega", "phi", "zeta" )] * d2r
  
  ##  transform initial anisotropy parameters
  
  transformed.aniso <- sapply(
    aniso.name,
    function( x, fwd.tf, aniso.tf, aniso ){
      fwd.tf[[aniso.tf[x]]]( aniso[x] )
    },
    fwd.tf = fwd.tf,
    aniso.tf = aniso.tf,
    aniso = aniso
  )
  names( transformed.aniso ) <- aniso.name 
  
  ##  create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  likelihood.calculations
  
  envir <- new.env()
  lik.item <- list()
  
  ##  initialize values of variogram parameters stored in the environment
  
  lik.item[["param"]]   <- param
  if( reparam ){
    lik.item[["eta"]] <- unname( lik.item[["param"]]["nugget"] )
    lik.item[["xi"]]  <- unname( lik.item[["param"]]["snugget"] )
  } else {
    lik.item[["eta"]] <- unname( t.param[["nugget"]] )
    lik.item[["xi"]]  <- unname( t.param[["snugget"]] )
  }
  
  lik.item[["aniso"]]   <- list( 
    isotropic = initial.objects[["isotropic"]], 
    aniso = aniso
  )
  names( lik.item[["aniso"]][["aniso"]] ) <- aniso.name
  
  lik.item[["aniso"]][["sincos"]] <- list(
    co = unname( cos( aniso["omega"] ) ),
    so = unname( sin( aniso["omega"] ) ),
    cp = unname( cos( aniso["phi"] ) ),
    sp = unname( sin( aniso["phi"] ) ),
    cz = unname( cos( aniso["zeta"] ) ),
    sz = unname( sin( aniso["zeta"] ) )
  )
  
  n <- NCOL( lag.vectors)
  
  if( n <= 3 ){
    
    lik.item[["aniso"]][["rotmat"]] <- with( 
      lik.item[["aniso"]][["sincos"]],
      rbind(
        c(             sp*so,             sp*co,       cp ),
        c( -cz*co + sz*cp*so,  co*sz*cp + cz*so,   -sp*sz ),
        c( -co*sz - cz*cp*so, -cz*co*cp + sz*so,    cz*sp )
      )[ 1:n, 1:n, drop = FALSE ]
    )
    
    
    lik.item[["aniso"]][["sclmat"]] <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1:n ]
    
  } else {  # only isotropic case for n > 3
    
    lik.item[["aniso"]][["rotmat"]] <- diag( n )
    lik.item[["aniso"]][["sclmat"]] <- rep( 1., n )
    
  }
  
  ## restore Valphaxi object from georob.object if available and variogram
  ## parameters are the same
  
  if( 
    !is.null( georob.object ) && 
    isTRUE( all.equal( georob.object[["param"]], f.reparam.bkw( t.param ) ) ) &&
    isTRUE( all.equal( 
        georob.object[["aniso"]][["aniso"]] * c( 1, 1, rep( d2r, 3 ) ),
        lik.item[["aniso"]][["aniso"]] ) 
    ) &&
    isTRUE( all.equal( 
        NROW( XX ), 
        length( georob.object[["Valphaxi.objects"]][["Valphaxi"]][["diag"]] ) )
    )
   ){
    lik.item[["Valphaxi"]]  <- expand( georob.object[["Valphaxi.objects"]] )
  } 
  
  assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
  ##  compute various expectations of psi, and its derivative etc.
  
  expectations <- numeric()
  
  ##  ... E[ psi'(x) ]  with respect to nominal Gaussian model
  
  if( is.null( rho.psi.etc[["exp.gauss.dpsi"]] ) ){
    
    t.exp <- integrate( 
      function( x, dpsi.function, tuning.psi ) {
        dnorm( x ) * dpsi.function( x, tuning.psi )
      }, 
      lower = -Inf, upper = Inf, 
      dpsi.function = rho.psi.etc[["dpsi.function"]], 
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["exp.gauss.dpsi"] <- t.exp[["value"]]
    
  } else {
    
    expectations["exp.gauss.dpsi"] <- rho.psi.etc[["exp.gauss.dpsi"]]( tuning.psi )
    
  }

  ##  ... E[ psi(x)^2 ]  with respect to nominal Gaussian model
  
  if( is.null( rho.psi.etc[["var.gauss.psi"]] ) ){
    
    t.exp <- integrate( 
      function( x, psi.function, tuning.psi ) {
        dnorm( x ) * ( psi.function( x, tuning.psi ) )^2
      }, 
      lower = -Inf, upper = Inf, 
      psi.function = rho.psi.etc[["psi.function"]],
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["var.gauss.psi"] <- t.exp[["value"]]
    
  } else {
    
    expectations["var.gauss.psi"] <- rho.psi.etc[["var.gauss.psi"]]( tuning.psi )
    
  }
  
  ## ...  E[ x^2 ] with respect to assumed longtailed distribution 
  ## f0 \propto 1/sigma exp( - rho(x/sigma) )
  
  if( is.null( rho.psi.etc[["var.f0.eps"]] ) ){

    t.exp <- integrate( 
      function( x, f0, tuning.psi ) {
        f0( x, tuning.psi, sigma = 1. ) * x^2
      }, 
      lower = -Inf, upper = Inf, 
      f0 = rho.psi.etc[["f0"]],
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["var.f0.eps"] <- t.exp[["value"]]

  } else {
    
    expectations["var.f0.eps"] <- rho.psi.etc[["var.f0.eps"]](
      tuning.psi, sigma = 1.  
    )
    
  }

  
  ## ...  E[ psi(x)^2 ] ( = E[ psi'(x) ]) with respect to assumed longtailed distribution 
  ## f0 \propto 1/sigma exp( - rho(x/sigma) )
  
  expectations["var.f0.psi"] <- rho.psi.etc[["var.f0.psi"]]( tuning.psi )
  
  if( verbose > 2 ) cat( 
    "\n expectation of psi'(epsilon/sigma) under nominal Gaussian model  :", 
    signif( expectations["exp.gauss.dpsi"] ), "\n", 
    "variance of psi(epsilon/sigma) under nominal Gaussian model      :", 
    signif( expectations["var.gauss.psi"] ), "\n", 
    "variance of epsilon under long-tailed model                      :", 
    signif( expectations["var.f0.eps"] ), "\n", 
    "variance of psi(epsilon/sigma) under long-tailed model           :", 
    signif( expectations["var.f0.psi"] ), "\n" 
  )
  
  
  ## zhat
  
  sel <- !is.na (betahat )
  zhat <- drop( XX[, sel, drop=FALSE] %*% betahat[sel] + bhat )
  names( zhat ) <- rownames( XX )
  
  r.hessian <- NULL
  
  if( tuning.psi < tuning.psi.nr ) {
    
    ## robust REML estimation
    
    hessian <- FALSE
      
    if( any( c( fit.param, fit.aniso ) ) ){
      
      ##  find roots of estimating equations
            
      r.root <- nleqslv(
        x = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        fn = estimating.equations.theta,
        method = control.nleqslv[["method"]],
        global = control.nleqslv[["global"]],
        xscalm = control.nleqslv[["xscalm"]],
        control = control.nleqslv[["control"]],
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        expectations = expectations,
        error.family.estimation = error.family.estimation,
        control.pmm = control.pmm,
        verbose = verbose
      ) 
      
      #       r.param <- r.root[["x"]] names( r.param ) <- names(
      #       transformed.param[ fit.param ] )
      
      r.gradient <- r.root[["fvec"]]
      names( r.gradient ) <- c(
        names( transformed.param[ fit.param ] ),
        names( transformed.aniso[ fit.aniso ] )
      )
      
      r.converged <- r.root[["termcd"]] == 1
      r.convergence.code <- r.root[["termcd"]] 
      
      r.counts <- c( nfcnt = r.root[["nfcnt"]], njcnt = r.root[["njcnt"]] )
            
    } else {
      
      ##  all variogram parameters are fixed
      
      ##  evaluate estimating equations
      
      r.gradient <- estimating.equations.theta(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        expectations = expectations,
        error.family.estimation = error.family.estimation,
        control.pmm = control.pmm,
        verbose = verbose
      )
      
      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )
      
    }
    
    r.opt.neg.loglik <- NA_real_
    
  } else {
    
    if( any( c( fit.param, fit.aniso ) ) ){
      
      ##  Gaussian REML estimation
      
      error.family.cov.effects <- error.family.cov.residuals <- "gaussian"
      
      if( identical( maximizer, "optim" ) ){
      
        r.opt.neg.restricted.loglik <- optim(
          par = c( 
            transformed.param[ fit.param ],
            transformed.aniso[ fit.aniso ] 
          ),
          fn = negative.loglikelihood,
          gr = gradient.negative.loglikelihood,
          method = control.optim[["method"]], 
          lower = control.optim[["lower"]],
          upper = control.optim[["upper"]],
          control = control.optim[["control"]],
          hessian = hessian && !reparam,
          envir = envir,        
          variogram.model = variogram.model,
          fixed.param = c( 
            transformed.param[ !fit.param ],
            transformed.aniso[ !fit.aniso ] 
          ),
          param.name = param.name, 
          aniso.name = aniso.name,
          param.tf = c( param.tf, aniso.tf ),
          deriv.fwd.tf = deriv.fwd.tf,
          bwd.tf = bwd.tf,
          safe.param = safe.param,
          reparam = reparam,
          lag.vectors = lag.vectors,
          XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
          yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
          psi.function = rho.psi.etc[["psi.function"]], 
          tuning.psi = tuning.psi,
          tuning.psi.nr = tuning.psi.nr,
          ml.method = ml.method,
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.ftol = irwls.ftol,
          control.pmm = control.pmm,
          verbose = verbose,
          force.gradient = force.gradient
        )
        
        r.opt.neg.loglik <- r.opt.neg.restricted.loglik[["value"]]     
        r.converged <- r.opt.neg.restricted.loglik[["convergence"]] == 0
        r.convergence.code <- r.opt.neg.restricted.loglik[["convergence"]]      
        r.counts <- r.opt.neg.restricted.loglik[["counts"]]
        
      } else {
        
        r.opt.neg.restricted.loglik <- nlminb(
          start = c( 
            transformed.param[ fit.param ],
            transformed.aniso[ fit.aniso ] 
          ),
          objective = negative.loglikelihood,
          gradient = gradient.negative.loglikelihood,
          control = control.nlminb[["control"]],
          lower = control.nlminb[["lower"]],
          upper = control.nlminb[["upper"]],
          envir = envir,        
          variogram.model = variogram.model,
          fixed.param = c( 
            transformed.param[ !fit.param ],
            transformed.aniso[ !fit.aniso ] 
          ),
          param.name = param.name, 
          aniso.name = aniso.name,
          param.tf = c( param.tf, aniso.tf ),
          deriv.fwd.tf = deriv.fwd.tf,
          bwd.tf = bwd.tf,
          safe.param = safe.param,
          reparam = reparam,
          lag.vectors = lag.vectors,
          XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
          yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
          psi.function = rho.psi.etc[["psi.function"]], 
          tuning.psi = tuning.psi,
          tuning.psi.nr = tuning.psi.nr,
          ml.method = ml.method,
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.ftol = irwls.ftol,
          control.pmm = control.pmm,
          verbose = verbose,
          force.gradient = force.gradient
        )
        
        r.opt.neg.loglik <- r.opt.neg.restricted.loglik[["objective"]]     
        r.converged <- r.opt.neg.restricted.loglik[["convergence"]] == 0
        r.convergence.code <- r.opt.neg.restricted.loglik[["convergence"]]      
        r.counts <- r.opt.neg.restricted.loglik[["evaluations"]]
        
      }
      
      r.gradient <- gradient.negative.loglikelihood(
        adjustable.param = r.opt.neg.restricted.loglik[["par"]],
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ],
          transformed.aniso[ !fit.aniso ] 
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        control.pmm = control.pmm,
        verbose = verbose
      )
      
    } else {
      
      ##  all variogram parameters are fixed
      
      hessian <- FALSE
      
      ##  compute negative restricted loglikelihood and gradient
      
      r.opt.neg.loglik <- negative.loglikelihood(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        control.pmm = control.pmm,
        verbose = verbose
      )
      
      r.gradient <- gradient.negative.loglikelihood(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        control.pmm = control.pmm,
        verbose = verbose
      )
      
      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )
      
    }
    
  }
  
  ##  get the other fitted items
    
  lik.item <- get( "lik.item", pos = as.environment( envir ) )
  
  ##  compute signal variance and re-compute original parameters
  
  if( reparam ){
    
    t.qmp <- NROW(XX)
    if( identical( ml.method, "REML" ) ){
      t.qmp <- t.qmp - col.rank.XX[["rank"]]  
    }
    t.sigma02 <- lik.item[["zhat"]][["RSS"]] / t.qmp
    
    lik.item[["param"]] <- f.reparam.bkw( c(variance = t.sigma02, lik.item$param) )
    fit.param <- c( variance = TRUE, fit.param )
    param.name <- c( "variance", param.name )
  
  }
  
  ## extract or recompute Hessian for Gaussian (RE)ML
  
  if( hessian ){
    
    if( reparam || identical( maximizer, "nlminb" ) ){
      
      t.sel <- match( param.name, names( all.param.tf ) )
      
      param.tf <- sapply(
         all.param.tf[t.sel],
        function( x ) if( length(x) > 1L ) x[variogram.model] else x
      )
      names( param.tf ) <- param.name 
      
      ##  transform initial variogram parameters
      
      transformed.param <- sapply(
        param.name,
        function( x, fwd.tf, param.tf, param ) fwd.tf[[param.tf[x]]]( param[x] ),
        fwd.tf = fwd.tf,
        param.tf = param.tf,
        param = lik.item[["param"]]
      )
      
      names( transformed.param ) <- param.name 
      
      r.hessian <- optimHess(
        par = c( 
          transformed.param[ fit.param ],
          transformed.aniso[ fit.aniso ] 
        ),
        fn = negative.loglikelihood,
        gr = gradient.negative.loglikelihood,
        control = control.optim[["control"]],
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ],
          transformed.aniso[ !fit.aniso ] 
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = c( param.tf, aniso.tf ),
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = FALSE,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        control.pmm = control.pmm,
        verbose = 0,
        force.gradient = force.gradient
      )    
      
    } else {

      r.hessian <- r.opt.neg.restricted.loglik[["hessian"]]
    
    }
  
  }
  
  
  ##  compute the covariances of fixed and random effects under nominal
  ##  Gaussian model
  
  r.cov <- list()
  
  if( any( c( 
        cov.bhat, cov.betahat, cov.bhat.betahat, 
        cov.delta.bhat, cov.delta.bhat.betahat, 
        aux.cov.pred.target
      ) 
    ) 
  ){
    
    ##  compute the covariances of fixed and random effects under nominal
    ##  Gaussian model
    
    r.cov <- covariances.fixed.random.effects(
      Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
      Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
      Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
      Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
      rweights = lik.item[["zhat"]][["rweights"]],
      XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
      nugget = lik.item[["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations, family = error.family.cov.effects,
      cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
      cov.betahat = cov.betahat, 
      cov.bhat.betahat = cov.bhat.betahat,
      cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
      cov.delta.bhat.betahat = cov.delta.bhat.betahat,
      cov.ehat = FALSE, full.cov.ehat = FALSE, 
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = aux.cov.pred.target,
      control.pmm = control.pmm,
      verbose = verbose
    )
    
    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0 ) cat(
        "\nan error occurred when computing the covariances of fixed and random effects\n"
      )
    } else {
      r.cov <- r.cov[!names(r.cov) %in% "error"]
    }
    
  }
  
  ##  compute covariances of residuals under assumed long-tailed error distribution 
  
  if( any( c( cov.ehat, cov.ehat.p.bhat ) ) ){
    
    ##  compute the covariances of fixed and random effects under nominal
    ##  Gaussian model
    
    r.cov <- c(
      r.cov,
      covariances.fixed.random.effects(
        Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
        Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
        Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
        Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
        rweights = lik.item[["zhat"]][["rweights"]],
        XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
        nugget = lik.item[["param"]]["nugget"],
        eta = lik.item[["eta"]],
        expectations = expectations, family = error.family.cov.residuals,
        cov.bhat = FALSE, full.cov.bhat = FALSE,
        cov.betahat = FALSE, 
        cov.bhat.betahat = FALSE,
        cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
        cov.delta.bhat.betahat = FALSE,
        cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat, 
        cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
        aux.cov.pred.target = FALSE,
        control.pmm = control.pmm,
        verbose = verbose
      )
    )
    
    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0 ) cat(
        "\nan error occurred when computing the covariances of residuals\n"
      )
    } else {
      r.cov <- r.cov[!names(r.cov) %in% "error"]
    }
    
  }
  
  ## stop SNOW and snowfall clusters
  
  #   f.stop.cluster()
  
  if( length( lik.item[["defaultCluster"]] ) > 0 ){
    cl <- lik.item[["defaultCluster"]]
    
    junk <- parLapply( cl, 1:length(cl), function( i ) sfStop() )
    junk <- stopCluster( cl )
    sfStop()
    options( error = NULL )  
  } 
  
  if( sfIsRunning() ){
    sfStop()
    options( error = NULL )  
  }
  
  if( file.exists( "SOCKcluster.RData" ) ){
    file.remove( "SOCKcluster.RData" )
  } 
  
  
  attr( r.gradient, "eeq.emp" )    <- lik.item[["eeq"]][["eeq.emp"]]
  attr( r.gradient, "eeq.exp" )    <- lik.item[["eeq"]][["eeq.exp"]]
  
  ##      ##  compute residual degrees of freedom 
  ##      
  ##      r.df <- f.compute.df( 
  ##          Valphaxi = lik.item[["Valphaxi"]][["Valphaxi"]],
  ##          XX = XX, 
  ##          param = lik.item[["param"]]
  ##      )
  
  result.list <- list(
    loglik = -r.opt.neg.loglik,
    variogram.model = variogram.model,
    param = lik.item[["param"]],
    aniso = lik.item[["aniso"]],
    gradient = r.gradient,
    tuning.psi = tuning.psi,
    coefficients = lik.item[["zhat"]][["betahat"]],
    fitted.values = drop( XX %*% lik.item[["zhat"]][["betahat"]] )[TT],
    bhat = lik.item[["zhat"]][["bhat"]],
    residuals = lik.item[["zhat"]][["residuals"]],
    rweights = lik.item[["zhat"]][["rweights"]],
    converged = r.converged,
    convergence.code = r.convergence.code,
    iter = r.counts,
    Tmat = TT
  )
  names( result.list[["fitted.values"]] ) <- names( result.list[["residuals"]] )
  names( result.list[["rweights"]] )      <- names( result.list[["residuals"]] )
  
  if( any( c( 
        cov.bhat, cov.betahat, cov.bhat.betahat, 
        cov.delta.bhat, cov.delta.bhat.betahat, 
        cov.ehat, cov.ehat.p.bhat, aux.cov.pred.target
      ) 
    ) 
  ){
    
    result.list[["cov"]] <- compress( r.cov )
    
  }
  
  ## map angles to halfcircle
  
  result.list[["aniso"]][["aniso"]] <- result.list[["aniso"]][["aniso"]] /
    c( rep( 1., 2 ), rep( d2r, 3 ) )
  
  if( !result.list[["aniso"]][["isotropic"]] ){
    
    if( result.list[["aniso"]][["aniso"]]["omega"] < 0. ){
      result.list[["aniso"]][["aniso"]]["omega"] <- 
      result.list[["aniso"]][["aniso"]]["omega"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["omega"] > 180. ){
      result.list[["aniso"]][["aniso"]]["omega"] <- 
      result.list[["aniso"]][["aniso"]]["omega"] - 180.
    }
    if( result.list[["aniso"]][["aniso"]]["phi"] < 0. ){
      result.list[["aniso"]][["aniso"]]["phi"] <- 
      result.list[["aniso"]][["aniso"]]["phi"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["phi"] > 180. ){
      result.list[["aniso"]][["aniso"]]["phi"] <- 
      result.list[["aniso"]][["aniso"]]["phi"] - 180.
    }
    if( result.list[["aniso"]][["aniso"]]["zeta"] < 90. ){
      result.list[["aniso"]][["aniso"]]["zeta"] <- 
      result.list[["aniso"]][["aniso"]]["zeta"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["zeta"] > 90. ){
      result.list[["aniso"]][["aniso"]]["zeta"] <- 
      result.list[["aniso"]][["aniso"]]["zeta"] - 180.
    }
  
  }
  
  result.list[["expectations"]]       <- expectations
  result.list[["Valphaxi.objects"]]   <- compress( lik.item[["Valphaxi"]] )
  result.list[["zhat.objects"]]       <- compress( 
    lik.item[["zhat"]][c( "Aalphaxi", "Palphaxi", "Valphaxi.inverse.Palphaxi" )]
  )
  
  result.list[["locations.objects"]] <- initial.objects[["locations.objects"]]
  result.list[["locations.objects"]][["lag.vectors"]] <- lag.vectors
  
  result.list[["initial.objects"]] <- list(
    coefficients = initial.objects[["betahat"]],
    bhat = initial.objects[["bhat"]],
    param = param,
    fit.param = fit.param,
    aniso = aniso,
    fit.aniso = fit.aniso
  )
  if( !is.null( r.hessian ) ){
    result.list[["hessian"]] <- r.hessian
  }
  ##      result.list[["df.model"]] <- r.df
  
  
  return(result.list)
  
}

#  ##############################################################################

getCall.georob <- 
  function( object )
{
  
  ## Function replaces the name of a formula object in the call component
  ## of a georob object by the formula itself (needed for update.default to
  ## work correctly)

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists

  if( is.null( call <- getElement( object, "call" ) ) )  stop(
    "need an object with call component"
  )
  call[["formula"]] <- update.formula( formula(object), formula( object ) )
  
  return( call )
  
}


################################################################################

f.aux.gamma <- 
  function( 
    lag.vectors, variogram.model, param, aniso
  )
{

  ## auxiliary function to compute semivariances for an anisotropic model
  
  ## arguments:
  
  ## param                                        vector with variogram parameters in standard order
  ## aniso                                        list with component rotmat and sclmat for coordinate
  ##                                                                  transformation in 3d
  ## lag.vectors    
  
  ## 2012-04-13 A. Papritz
  ## 2012-05-23 ap correction in model.list for models with more than 4 parameters
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-03-05 AP changes for version 3 of RandomFields
  ## 2015-07-17 new function name

  
  ## matrix for coordinate transformation
  
  A <- with(
    aniso,
    sclmat * rotmat / param["scale"]
  )
  
  ## set up variogram model object
  
  model.list <- list( "+",
    list( "$", 
      var = param["variance"], 
      A = A, 
      if( length( param[-(1:4)] ) > 0 ){
         c( list( variogram.model ) , as.list(param[-(1:4)]) )
      } else {
        list( variogram.model )
      }
    ),
    list( "$", 
      var = sum( param[ c("nugget", "snugget") ] ), 
      list( "nugget" )
    )
  )
  
  ##  semivariance 
  
  ## functions of version 3 of RandomFields
  
  RFoptions(newAniso=FALSE)
  
  r.gamma <- try(
    RFvariogram(
      x = lag.vectors, model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
    ),
    silent = TRUE
  )

  ## functions of version 3 of RandomFields

  ##   RFoldstyle()
  ##   r.gamma <- try(
  ##     Variogram( lag.vectors, model = model.list ),
  ##     silent = TRUE
  ##   )
  
  return( r.gamma )
  
}

################################################################################

f.aux.eeq <- function( 
  x, 
  nme.var.snug,
  param, aniso, xi,
  Valphaxii, Valphaxii.cov.bhat,
  bh, bhVaxi,
  r.cov, lik.item,
  TtT, 
  lag.vectors, variogram.model, control.pmm, verbose
){
  
  ##  auxiliary function to compute robustified estimating equations
  ##  (called by estimating.equations.theta)
  
  ## 2014-07-29 A. Papritz
  ## 2015-03-03 AP changes to optimize computing effort
  ## 2015-07-17 AP new function interface and improved efficiency of computation
  ## 2015-07-27 AP changes to further improve efficiency
  
  switch(
    x,
    nugget = {
      
      ## nugget
      
      #       eeq.exp <- sum( diag( 
      #           ( 1/TtT * lik.item[["Valphaxi"]][["Valphaxi.inverse"]] ) %*% r.cov[["cov.bhat"]] %*% lik.item[["Valphaxi"]][["Valphaxi.inverse"]]
      #         ) 
      #       )
      #       eeq.emp <- sum( 
      #         ( lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )^2 / TtT
      #       )
      
      eeq.exp <- sum( 1/TtT * Valphaxii * Valphaxii.cov.bhat )
      eeq.emp <- sum( bhVaxi^2 / TtT )
      
    },
    snugget = {
      
      ## snaugget
      
      #       eeq.exp <- sum( 
      #         diag( 
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% r.cov[["cov.bhat"]]
      #         )
      #       )
      #       eeq.emp <- sum( lik.item[["zhat"]][["Valphaxi.inverse.bhat"]]^2 )
      
      eeq.exp <- sum( Valphaxii * Valphaxii.cov.bhat )
      eeq.emp <- sum( bhVaxi^2 )
            
    },
    variance = {
      
      ## variance
      
      #       eeq.exp <- sum( 
      #         diag(
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% lik.item[["Valphaxi"]][["Valpha"]] %*% 
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% r.cov[["cov.bhat"]]
      #         )
      #       )
      #       eeq.emp <- sum( 
      #         lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] * drop( lik.item[["Valphaxi"]][["Valpha"]] %*% lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )
      #       )
            
      if( xi > 0. ){
        
        aux <- -xi * Valphaxii
        diag( aux ) <- diag( aux ) + 1.
        eeq.exp <- sum( aux * Valphaxii.cov.bhat ) / ( 1. - xi )
        eeq.emp <- sum( bhVaxi * ( bh - xi * bhVaxi ) ) / ( 1. - xi )
      
      } else {
        
        eeq.exp <- sum( diag( Valphaxii.cov.bhat ) )
        eeq.emp <- sum( bh * bhVaxi )
        
      }
      
    },
    {
      
      ## scale and extra parameters 
      
      dVa <- partial.derivatives.variogram(
        x = lag.vectors, variogram.model = variogram.model, param = param, 
        d.param = x,
        aniso = aniso,
        verbose = verbose
      )
      #       aux <- pmm( dVa, lik.item[["Valphaxi"]][["Valphaxi.inverse"]], control.pmm )
      #       eeq.exp <- sum(
      #         pmm( lik.item[["Valphaxi"]][["Valphaxi.inverse"]], aux, control.pmm ) * r.cov[["cov.bhat"]]
      #       )
      #       eeq.emp <- sum( 
      #         lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] * drop( dVa %*% lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )
      #       )
      
      aux <- pmm( dVa, Valphaxii, control.pmm )
      eeq.exp <- sum( aux * Valphaxii.cov.bhat ) * param["variance"]
      eeq.emp <- sum( bhVaxi * drop( dVa %*% bhVaxi ) ) * param["variance"]
    }
  )
  
  c( eeq.exp = eeq.exp, eeq.emp = eeq.emp )
  
}

################################################################################

f.aux.gradient.nll <- function( 
  x, 
  nme.var.snug,
  param, aniso, 
  TT, TtT, XX, 
  res, 
  Qi, Vi, Qt11Vi,
  bh, bhVi, 
  lag.vectors, variogram.model,
  param.tf, deriv.fwd.tf, 
  ml.method,
  control.pmm, verbose
){
  
  ##  auxiliary function to compute gradient of (restricted) log-likelihood
  ##  (called by gradient.negative.loglikelihood)
  ##  original parametrization of variogram
  
  ## 2014-07-29 A. Papritz
  ## 2015-07-17 AP new parametrization of loglikelihood
  ## 2015-07-17 AP new function interface, improve efficiency
  ## 2015-07-27 AP changes to further improve efficiency
  
  t.q <- NROW(Vi)
  
  switch(
    x,
    nugget = {
      
      ## compute partial derivative of (restricted) log-likelihood with
      ## respect to nugget
      
      ## partial derivate of U with respect to nugget
      
      Ttres <- as.vector( tapply( res, factor( TT ), sum ) )
      dU <- -sum( Ttres^2 / TtT ) / param["nugget"]^2
      
      ## partial derivative of log(det(Q)) with respect to nugget
      
      if( identical( ml.method, "REML" ) ){
        TtTX <- TtT * XX
        dlogdetQ <- -sum( 
          Qi * rbind( 
            cbind( diag( TtT ),                TtTX   ),
            cbind(     t(TtTX), crossprod( XX, TtTX ) )
          )
        ) / param["nugget"]^2
      } else {
        dlogdetQ <- -sum( TtT * diag(Qi) ) / param["nugget"]^2
      }
      
      ## partial derivative of loglik with respect to transformed nugget
      
      result <- c( nugget = ( -0.5 * ( t.q / param["nugget"] + dlogdetQ + dU ) ) / 
        deriv.fwd.tf[[param.tf["nugget"]]]( param["nugget"] )
      )
            
    },
    snugget = {
      
      ##  compute partial derivative of (restricted) log-likelihood with
      ##  respect to spatial nugget
      
      ## partial derivate of U with respect to spatial nugget
      
      dU <- -sum( bhVi^2 )
      
      ## partial derivative of V with respect to spatial nugget
      
      dlogdetV <- sum( diag( Vi ) )
      
      ## partial derivative of log( det( Q ) ) with respect to spatial nugget
      
      dlogdetQ <- -sum( Qt11Vi * Vi )
      
      ## partial derivative of loglik with respect to transformed spatial
      ## nugget
      
      result <- c(
        snugget = ( -0.5 * ( dlogdetV + dlogdetQ + dU ) ) / 
        deriv.fwd.tf[[param.tf["snugget"]]]( param["snugget"] )
      )
      
    },
    variance = {
    
      ##  compute partial derivative of restricted log-likelihood with
      ##  respect to variance
      
      if( param["snugget"] > 0. ){
        
        aux <- -param["snugget"] * Vi
        diag(aux) <- diag( aux ) + 1.
      
        ## partial derivate of U with respect to variance
        
        dU <- -sum( bhVi * ( bh - param["snugget"] * bhVi ) )
        
        # partial derivative of V with respect to variance
        
        dlogdetV <- sum( diag( aux ) )
        
        ## partial derivative of log(det Q)) with respect to variance

        dlogdetQ <- -sum( Qt11Vi * aux )
        
      } else {
      
        ## partial derivate of U with respect to variance
        
        dU <- -sum( bhVi * bh )
        
        # partial derivative of V with respect to variance
        
        dlogdetV <- NROW( Vi )
        
        ## partial derivative of log(det Q)) with respect to variance

        dlogdetQ <- -sum( diag( Qt11Vi ) )
      }
                  
      ## partial derivative of loglik with respect to transformed variance
      
      result <- c(
        variance = ( -0.5 * ( dlogdetV + dlogdetQ + dU ) / param["variance"] ) / 
        deriv.fwd.tf[[param.tf["variance"]]]( param["variance"] )
      )
      
    },
    {
      
      ## compute partial derivative of (restricted) log-likelihood with
      ## respect to scale and extra parameters
      
      
      ## partial derivative of V_0(alpha) with respect to scale and extra
      ## parameters

      dVa <- partial.derivatives.variogram(
        x = lag.vectors, variogram.model = variogram.model, param = param, 
        d.param = x,
        aniso = aniso,
        verbose = verbose
      )
      
      dVaVi <- pmm( dVa, Vi, control.pmm ) 
      
      ##  derivate of U
      
      dU <- -sum( bhVi * drop( dVa %*% bhVi ) )

      ## partial derivatie of V
      
      dlogdetV <- sum( diag( dVaVi ) )
      
      ##  derivative of log(det(Q))
      
      dlogdetQ <- -sum( Qt11Vi * t( dVaVi ) )
      
      ## partial derivative of loglik
        
      result <- ( -0.5 * ( dlogdetV + dlogdetQ + dU ) * param["variance"] ) / 
        deriv.fwd.tf[[param.tf[x]]]( 
          c( param, aniso[["aniso"]] )[x] 
        )
      names( result ) <- x
            
    }
  )
  
  return( result )
  
}

################################################################################

f.aux.gradient.npll <- function( 
  x, 
  nme.var.snug,
  param, aniso,
  eta, xi,
  TT, TtT, XX, col.rank.XX,
  res, Ustar,
  Qsi, Qst11Vai, 
  Valphaxii, 
  bh, bhVaxi,
  lag.vectors, variogram.model,
  param.tf, deriv.fwd.tf, 
  ml.method,
  control.pmm, verbose
){
  
  ##  auxiliary function to compute gradient of (restricted) profile
  ##  log-likelihood (called by gradient.negative.loglikelihood)
  ## reparametrized variogram
  
  ## 2015-07-17 A. Papritz
  ## 2015-07-27 AP changes to improve efficiency
  
  t.q <- t.qmp <- NROW(XX)
  if( identical( ml.method, "REML" ) ){
    t.qmp <- t.q - col.rank.XX[["rank"]]  
  }
  
  switch(
    x,
    nugget = {
      
      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to eta
      
      ## partial derivative of log(Ustar) with respect to eta

      Tt.res <- as.vector( tapply( res, factor( TT ), sum ) )
      dU <- sum( Tt.res^2 / TtT ) / Ustar
      
      ## partial derivative of log(det(Qstar)) with respect to eta

      if( identical( ml.method, "REML" ) ){
        TtTX <- TtT * XX
        dlogdetQ <- sum( 
          Qsi * rbind( 
            cbind( diag( TtT ),                TtTX   ),
            cbind(     t(TtTX), crossprod( XX, TtTX ) )
          )
        )
      } else {
        dlogdetQ <- sum( TtT * diag(Qsi) )
      }
      
      ## partial derivative of profile loglik with respect to transformed eta
      
      result <- c( nugget = -0.5 * ( t.qmp * dU - t.q / eta + dlogdetQ ) / 
        deriv.fwd.tf[[param.tf["nugget"]]]( param["nugget"] )
      )
            
    },
    snugget = {
      
      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to xi
      
      ## partial derivative of log(Ustar) with respect to xi

      dU <- -sum( bhVaxi * ( bhVaxi - bh) ) / Ustar
      
      ## partial derivative of log(det(Valphaxi)) with respect to xi
      
      dlogdetV <- sum( diag( Valphaxii ) ) - NROW( Valphaxii )
      
      ## partial derivative of log(det(Qstar)) with respect to xi
      
      dlogdetQ <- -( sum( Qst11Vai * Valphaxii ) - sum( diag( Qst11Vai ) ) )
      
      ## partial derivative of profile loglik with respect to transformed xi
      
      result <- c( 
        snugget = -0.5 / ( 1. -xi ) * ( t.qmp * dU + dlogdetV + dlogdetQ ) / 
        deriv.fwd.tf[[param.tf["snugget"]]]( param["snugget"] )
      )
      
    },
    {
      
      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to scale and extra parameters
      
      
      ## partial derivative of V_0(alpha) with respect to scale and extra
      ## parameters
      
      dVa <- partial.derivatives.variogram(
        x = lag.vectors, variogram.model = variogram.model, param = param, 
        d.param = x,
        aniso = aniso,
        verbose = verbose
      )
      
      dVaVai <- pmm( dVa, Valphaxii, control.pmm )  ### !!!dVaVai not symmetric!!!!
      
      ## partial derivative of log(Ustar) (up to factor (1-xi))

      dU <- -sum( bhVaxi * drop( dVa %*% bhVaxi ) ) / Ustar
      
      ## partial derivative of log(det(Valphaxi)) (up to factor (1-xi))
      
      dlogdetV <- sum( diag( dVaVai ) )
      
      ## partial derivative of log(det(Qstar)) (up to factor (1-xi))
      
      dlogdetQ <- -sum( Qst11Vai * t( dVaVai ) )
      
      ## partial derivative of profile loglik
      
      result <- -0.5 * ( 1. - xi ) * ( t.qmp * dU + dlogdetV + dlogdetQ ) / 
        deriv.fwd.tf[[param.tf[x]]]( 
          c( param, aniso[["aniso"]] )[x] 
        )
      names( result ) <- x
            
    }
  )
  
  return( result )
  
}

################################################################################

f.aux.Qstar <- function( 
  TT, TtT, 
  XX, col.rank.XX, min.condnum, 
  Vi, 
  eta, 
  ml.method, control.pmm 
){
  
  ## auxiliary function to compute matrix Qstar used for Gaussian
  ## log-likelihood (called by likelihood.calculations)
  
  ## 2014-07-29 A. Papritz
  ## 2015-07-17 AP new function interface and new name
  
  result <- list( error = TRUE, log.det.Qstar = NULL, Qstar.inverse = NULL )
  
  TtTX <- eta * TtT * XX 
  
  ##  compute matrix Qstar
  
  Qstar <-  Vi
  diag( Qstar ) <- diag( Qstar ) + eta *TtT
  
  if( identical( ml.method, "REML" ) ){
    Qstar <- rbind( 
      cbind( Qstar,       TtTX                 ),
      cbind( t(TtTX), crossprod( XX, TtTX) )
    )
  }
  
  if( col.rank.XX[["deficient"]] && ml.method == "REML" ){
    
    ## compute log(pseudo.det(Qstar)) and (Moore-Penrose) pseudo inverse of Qstar by svd
    
    result[["error"]] <- FALSE
    s <- svd( Qstar )
    sel <- s[["d"]] / max( s[["d"]] ) > min.condnum
    #     result[["log.det.Qstar"]] <- sum( log( s[["d"]][s[["d"]] / max( s[["d"]] ) > min.condnum] ) )
    #     s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    #         result[["Qstar.inverse"]] <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )
    result[["log.det.Qstar"]] <- sum( log( s[["d"]][sel] ) )
    s[["d"]] <- ifelse( sel,  1. / s[["d"]], 0. )
    result[["Qstar.inverse"]] <- pmm(
      s[["v"]], s[["d"]] * t( s[["u"]] ), control.pmm
    )
    
  } else {
    
    ##  compute log(det(Qstar)) and inverse of Qstar by cholesky decomposition
    
    t.chol <- try( chol( Qstar ), silent = TRUE )
    
    if( !identical( class( t.chol ), "try-error" ) ) {
      
      result[["error"]] <- FALSE
      result[["log.det.Qstar"]] <- 2 * sum( log( diag( t.chol) ) )
      result[["Qstar.inverse"]] <- chol2inv( t.chol )
      
    }
    
  }
  
  result
  
}


################################################################################

f.aux.Valphaxi <- function(
  lag.vectors, variogram.model, param, xi, aniso, control.pmm, verbose
){
  
  ## auxiliary function to compute generalized correlation matrix and
  ## related items (called by likelihood.calculations)
  
  ## 2014-07-29 A. Papritz
  ## 2015-07-23 AP Valpha (correlation matrix without spatial nugget) no longer stored
  
  result <- list( 
    error = TRUE, gcr.constant = NULL, Valphaxi = NULL, 
    Valphaxi.inverse = NULL
  )
  
  cormat <- f.aux.gcr(
    lag.vectors = lag.vectors, variogram.model = variogram.model, 
    param = param, xi = xi, aniso = aniso, 
    control.pmm = control.pmm,
    verbose = verbose
  )
  if( cormat[["error"]] ) return( result )
  
  t.vchol <- try( chol( cormat[["Valphaxi"]] ), silent = TRUE )

  if( !identical( class( t.vchol ), "try-error" ) ) {

    result[["error"]]            <- FALSE
    result[["gcr.constant"]]     <- cormat[["gcr.constant"]]
    result[["Valphaxi"]]         <- cormat[["Valphaxi"]]
    #     result[["Valpha"]]           <- cormat[["Valpha"]]
    result[["Valphaxi.inverse"]] <- chol2inv( t.vchol )
    result[["log.det.Valphaxi"]] <- 2 * sum( log( diag( t.vchol) ) )
    
    attr( result[["Valphaxi"]], "struc" )         <- "sym"
    attr( result[["Valphaxi.inverse"]], "struc" ) <- "sym"
  }
  result
}


################################################################################

f.stop.cluster <- function( cl = NULL ){
  
  ## function to stop snow and snowfall clusters 
  ## 2014-07-31 A. Papritz
  
  if( sfIsRunning() ){
    sfStop()
  }
  
  if( file.exists( "SOCKcluster.RData" ) ){
    if( is.null( cl ) ) load( "SOCKcluster.RData" )
    file.remove( "SOCKcluster.RData" )
  } 
  
  ## stop cluster started by child processes in recursive paralellized
  ## computations
  
  if( !is.null( cl ) ){
    junk <- parLapply( cl, 1:length(cl), function( i ) sfStop() )
    junk <- stopCluster( cl )
  }
  options( error = NULL )  

}



################################################################################

f.psi.function <- function( x = c( "logistic", "t.dist", "huber" ), tp ){
  
  ## define psi-function and derivatives
  ## 2015-03-16 A. Papritz 
  ## 2015-08-19 AP pdf and variances of eps and psi(eps/sigma) for long-tailed error distribution, 
  ##               new parametrisation for t-dist family
  
  x <- match.arg( x )
  
  switch(
    x,
    logistic = list(
      #       rho.function = function( x, tuning.psi ) {
      #         tuning.psi * (-x + tuning.psi * 
      #           ( -log(2.) + log( 1. + exp( 2. * x / tuning.psi ) ) )
      #         )
      #       }, 
      psi.function = function( x, tuning.psi ) {
        tuning.psi * tanh( x / tuning.psi )
      }, 
      dpsi.function = function( x, tuning.psi ) {
        1. / (cosh( x / tuning.psi ))^2
      },
      f0 = if( is.finite( gamma((1. + tp^2)/2.) ) ){
        function( x, tuning.psi, sigma = 1. ) {              # pdf f0 propto 1/sigma exp(-rho(eps/sigma))
          ifelse( 
            is.finite( exp((2.*x)/(tuning.psi*sigma)) ),
            ( exp( (tuning.psi * ( x - tuning.psi * sigma * log( 
                      (1. + exp((2.*x)/(tuning.psi*sigma)))/2.
                    )))/sigma ) * gamma((1. + tuning.psi^2)/2.)) / 
            ( tuning.psi * sqrt(pi) * sigma * gamma(tuning.psi^2/2.) ),
            0.
          )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = NULL,                        # variance of eps under f0 (= expectation of psi'(eps/sigma) under f0)
      var.f0.psi = function( tuning.psi){       # variance of psi(eps/sigma) under f0
        tuning.psi^2 / ( (1. + tuning.psi^2) )
      },
      exp.gauss.dpsi = NULL,                    # expectation of psi'(eps/sigma) under N(0, sigma^2)
      var.gauss.psi = NULL                      # variance of psi(eps/sigma) under N(0, sigma^2)
    ),
    t.dist = list(
      #       rho.function = function( x, tuning.psi ){
      #         0.5 * tuning.psi^2 * log( (x^2 + tuning.psi^2) / tuning.psi^2 )
      #       },            
      psi.function = function( x, tuning.psi ){
        tuning.psi^2 * x / ( x^2 + tuning.psi^2 )
      },
      dpsi.function = function( x, tuning.psi ) {
        tuning.psi^2 * ( tuning.psi^2 - x^2 ) / ( x^2 + tuning.psi^2 )^2
      },
      f0 = if( is.finite( gamma((-1.+tp^2)/2.) ) ){
        function( x, tuning.psi, sigma = 1. ) {
          ( (tuning.psi^2)^(tuning.psi^2/2.) * gamma( tuning.psi^2/2.) ) / 
          ( tuning.psi*sqrt(pi) * (tuning.psi^2 + (x/sigma)^2 )^(tuning.psi^2/2.) * sigma * 
            gamma((-1.+tuning.psi^2)/2.) )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = if( tp > sqrt(3.) ){
        function( tuning.psi, sigma = 1. ){
          ( tuning.psi*sigma )^2 / ( -3. + tuning.psi^2 )
        } 
      } else NULL,
      var.f0.psi = function( tuning.psi ){
        1. - 3./(2. + tuning.psi^2)
      },
      exp.gauss.dpsi = if( is.finite( exp(tp^2/2.) ) ){
        function( tuning.psi ){
          tuning.psi^2 - tuning.psi^3 * exp(tuning.psi^2/2.) * sqrt(pi/2.) * 
          2. * (1.-pnorm(tuning.psi))
        }
      } else {
        function( tuning.psi ){ 1. }
      },
      var.gauss.psi = if( is.finite( exp(tp^2/2.) ) ){
        function( tuning.psi ){
          ( tuning.psi^3 * (-2.*tuning.psi + (1 + tuning.psi^2) * exp(tuning.psi^2/2.) * 
              sqrt(2.*pi) * 2. * (1.-pnorm(tuning.psi)) ) ) / 4.
        }
      } else {
        function( tuning.psi ){ 1. }
      }      
      ),
    huber = list(
      #       rho.function = function( x, tuning.psi ) {
      #         ifelse( 
      #           abs(x) <= tuning.psi, 
      #           0.5 * x^2, 
      #           tuning.psi * abs(x) - 0.5 * tuning.psi^2 
      #         )
      #       },
      psi.function = function( x, tuning.psi ) {
        ifelse( 
          abs(x) <= tuning.psi, 
          x, 
          sign(x) * tuning.psi
        )
      },
      dpsi.function = function( x, tuning.psi ) {
        ifelse( abs(x) <= tuning.psi, 1., 0. )
      },
      f0 = if( is.finite( (tp*exp(0.5*tp^2) ) ) ){
        function( x, tuning.psi, sigma = 1. ){
          1. / ( exp(
              ifelse( 
                abs(x/sigma) <= tuning.psi, 
                0.5 * (x/sigma)^2, 
                tuning.psi * abs(x/sigma) - 0.5 * tuning.psi^2 
              )) 
            * ( (2.*sigma) / (tuning.psi*exp(0.5*tuning.psi^2) ) +
              sqrt(2*pi) * sigma * (2.*pnorm(tuning.psi)-1.))
          )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = function( tuning.psi, sigma = 1. ){
        sigma^2 * ( 1. + 
          ( sqrt(2./pi) * ( 2. + tuning.psi^2 ) ) / 
          ( tuning.psi^2 * ( sqrt(2./pi) + tuning.psi *exp(0.5*tuning.psi^2) * (2.*pnorm(tuning.psi) - 1.) ) )
        )      
      },
      var.f0.psi = function( tuning.psi ){
        1. + 1. / (
          -1. - 0.5*sqrt(2*pi) * tuning.psi * exp( 0.5*tuning.psi^2 ) * (2.*pnorm(tuning.psi) - 1.) 
        )
      },
      exp.gauss.dpsi = function( tuning.psi ){
        (2.*pnorm(tuning.psi) - 1.)
      },
      var.gauss.psi = function( tuning.psi ){
        (2.*pnorm(tuning.psi) - 1.) +  tuning.psi*(-(sqrt(2./pi)/exp(tuning.psi^2/2.)) +  
          tuning.psi * 2. * (1.-pnorm(tuning.psi)) )
      }
    )
  )
  
}

################################################################################

## functions to (back)transform variogram parameters for alternative
## parametrization of variogram for Gaussian (RE)ML

## 2015-07-17 A. Papritz

f.reparam.fwd <- function( param ){
  param["variance"] <- tmp <- sum( param[c( "snugget", "variance" ) ] ) #sigma_0^2
  param["nugget"]   <- tmp / param["nugget"]    # eta
  param["snugget"]  <- param["snugget"] / tmp   # xi
  param  
}

f.reparam.bkw <- function( param ){
  param["nugget"]   <- 1. / param["nugget"] * param["variance"]   # tau^2
  param["snugget"]  <- param["snugget"] * param[ "variance" ]     # sigma_1^2
  param["variance"] <- param["variance"] - param["snugget"]       # sigma_2^2
  param  
}


################################################################################

f.aux.RSS <- function( res, TT, TtT, bhat, Valphaxi.inverse.bhat, eta ){

  ## function computes residual sums of squares required for
  ## likelihood computations
  
  ## 2015-07-17 A. Papritz
  
  Ttres <- res
  
  if( sum( duplicated( TT ) > 0 ) ){
    Ttres <- as.vector( tapply( Ttres, factor( TT ), sum ) )
  }
  
  eta * sum( Ttres^2 / TtT ) + sum( bhat * Valphaxi.inverse.bhat )

}


##  ###########################################################################

## auxiliary function to extract diagonal of square matrix and to return
## x unchanged otherwise

f.diag <- function( x ){
  switch( 
    class( x ),
    matrix = diag( x ),
    x
  )
}


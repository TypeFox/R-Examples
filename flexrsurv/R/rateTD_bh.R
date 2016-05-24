# idem rateTD.R but baseline hasard is gamma0(t), not exp(gamma0(t))
rateTD_bh_beta0alphabeta<- function(T, iT, gamma0, Zbeta0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # spline bases for baseline hazard
#  YT0 <- bs(T, knots=Knots_t0, intercept=Intercept_t0, degree=degree_t0, Boundary.knots =  Boundary.knots_t0)
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
    if(!is.null(Zbeta0)) {
      if(!is.null(Zalphabeta)) {
        (YT0 %*% gamma0 )* 
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
                fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
      else {
        (YT0 %*% gamma0) *
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
      }
    }
    else if(!is.null(Zalphabeta)) {
        (YT0 %*% gamma0)*
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
}

rateTD_bh_alphabeta<- function(T, iT, gamma0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
 
  # spline bases for baseline hazard
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
    YT  <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value 
  (YT0 %*% gamma0) * exp( YT %*% Zalphabeta[iT,])
  
}

rateTD_bh<- function(T, iT, gamma0, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE){
  # compute the contribution of the baseline hazard rate to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
 
  # spline bases for baseline hazard
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # returned value

  YT0 %*% gamma0
  
}



######################################################################
# same as rateTD* but witout the baseline


# idem rateTD.R but baseline hasard is gamma0(t), not exp(gamma0(t))
ratioTD_bh_beta0alphabeta<- function(T, iT, Zbeta0, Zalphabeta, 
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # spline bases for baseline hazard
  # spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
    if(!is.null(Zbeta0)) {
      if(!is.null(Zalphabeta)) {
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
                fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
      else {
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
      }
    }
    else if(!is.null(Zalphabeta)) {
            exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
}

ratioTD_bh_alphabeta<- function(T, iT, Zalphabeta, 
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
 
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value 
  exp( YT %*% Zalphabeta[iT,])
  
}

ratioTD_bh<- function(T, iT){
  # compute the contribution of the baseline hazard rate to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
 
  # returned value

  return(rep(1, length(T)))
  
}






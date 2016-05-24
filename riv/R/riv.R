riv <- function(Y, Xend, Xex=NULL, Zinst, dummies=NULL,
                method = c('S-est', 'SD-est', 
                           'MCD-est', 'classical')) {

  method <- match.arg(method)
  
  if (is.null(dummies)) {
    riv_noDummies(Y, Xend, Xex, Zinst,
                  intercept=TRUE,
                  method=method)
  } else {
    riv_withDummies(Y, Xend, Xex, Zinst, dummies,
                    method = method)
  }
}

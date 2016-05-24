riv_noDummies <- function(Y, Xend, Xex = NULL, Zinst, intercept = TRUE,
                          method = c('S-est', 'SD-est',
                                     'MCD-est', 'classical')) {
  Zinst <- as.matrix(Zinst)
  if (is.null(colnames(Zinst)))
      colnames(Zinst) <- paste('Zinst',
                               seq(ncol(Zinst)),
                               sep='')
  
  Xend <- as.matrix(Xend)
  if (is.null(colnames(Xend)))
      colnames(Xend) <- paste('Xend',
                              seq(ncol(Xend)),
                              sep='')
  
  kend <- ncol(Xend)
  k <- ncol(Zinst)
  if (kend > k)
    stop('the number of instruments must be equal or bigger to the number of endogenous variables')
  
  if (!is.null(Xex)) {
    Xex <- as.matrix(Xex)
    if (is.null(colnames(Xex)))
      colnames(Xex) <- paste('Xex',
                             seq(ncol(Xex)),
                             sep='')
  }
  
  if (any(is.na(c(Y, Xend, Xex, Zinst))))
    stop('missing values are not allowed')
  
  method <- match.arg(method)

  switch(method,
         `S-est`     = riv_sest(Y, Xend, Xex, Zinst, intercept),
         `SD-est`    = riv_sdest(Y, Xend, Xex, Zinst, intercept),
         `MCD-est`   = riv_mcdest(Y, Xend, Xex, Zinst, intercept),
         `classical` = riv_classical(Y, Xend, Xex, Zinst, intercept))
}

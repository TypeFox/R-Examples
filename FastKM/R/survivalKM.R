#----------------------------------------------------------------------#
# For survival trait type, iterate over coxKM methods until a stable   #
# solution is obtained.                                                #
# Return number of eigenvectors kept and the p.values.                 #
#----------------------------------------------------------------------#
survivalKM <- function(A, 
                       cKernel,  
                       y,  
                       delta,  
                       aEV,  
                       bEV, ...) {

  #------------------------------------------------------------------#
  # Load coxKM and survival packages.                                #
  #------------------------------------------------------------------#
  hasCoxKM <- requireNamespace("coxKM", quietly=TRUE)
  if( !hasCoxKM ) {
    stop(paste("R package coxKM is required for survival traits. ",
               "See documentation for further information.", sep=""))
  }

  hasSurvival <- requireNamespace("survival", quietly=TRUE)
  if( !hasSurvival ) {
    stop("R package survival is required for survival traits. ")
  }

  #------------------------------------------------------------------#
  # Match user inputs to coxph with formal arguments of coxph.       #
  #------------------------------------------------------------------#
  userArgs <- list(...)

  coxFormalNames <- names(formals(survival::coxph))
  coxFormals <- list()

  for( i in names(userArgs) ) {

    if( i %in% coxFormalNames ) coxFormals[[i]] <- userArgs[[i]]

  }

  dt <- data.frame(A, y, delta)
  coxFormals[["data"]] <- dt

  #-----------------------------------------------------------------#
  # Match user inputs to coxKM with formal arguments of coxKM.      #
  #-----------------------------------------------------------------#
  sFormalNames <- names(formals(coxKM::coxKM))
  sFormals <- list()

  for( i in names(userArgs) ) {

    if( i %in% sFormalNames ) sFormals[[i]] <- userArgs[[i]]

  }

  sFormals[["Z"]] <- NULL
  sFormals[["U"]] <- y
  sFormals[["Delta"]] <- delta
  sFormals[["kernel"]] <- cKernel
  sFormals[["is_check_genotype"]] <- FALSE
  sFormals[["weights"]] <- NULL

  nZa <- as.integer(round(length(aEV),0L))
  nZb <- as.integer(round(length(bEV),0L))
  nX <- ncol(A) - nZa - nZb

  obj <- survival::Surv(y,delta)

  cat("Running coxKM.\n")

  while(TRUE) {

    #--------------------------------------------------------------#
    # Attempt to calculate model parameters.                       #
    #--------------------------------------------------------------#
    nms <- paste("obj~",
                 paste(colnames(dt)[1:ncol(A)],collapse="+"),
                 sep="")
    coxFormals[["formula"]] <- stats::as.formula(nms)

    gamma <- tryCatch(expr =  do.call(what = survival::coxph,
                                      args = coxFormals),
                      error = function(e){NULL},
                      warning = function(w){NULL})

    if( is(gamma, "NULL") ) {
      #----------------------------------------------------------#
      # If tryCatch returned NULL, remove a column from A.       #
      #----------------------------------------------------------#
      res <- shrinkA(A = A, 
                     aEV = aEV, 
                     bEV = bEV, 
                     nZa = nZa, 
                     nZb = nZb, 
                     nX = nX)

      if( is(res,"NULL") ) {
        #------------------------------------------------------#
        # If res is NULL, no stable solution. Set p.value to   #
        # NULL and exit.                                       #
        #------------------------------------------------------#
        p.value <- NULL
        break
      } else {
        A = res$A
        nZa = res$nZa
        nZb = res$nZb
        rm(res)
        next
      } 
    }

    #--------------------------------------------------------------#
    # Attempt to calculate kernel association test.                #
    #--------------------------------------------------------------#
    sFormals[["X"]] <- A
    sFormals[["gamma"]] <- gamma$coef
    sFormals[["npert.check"]] <- FALSE

    SKATResult <- tryCatch(expr = do.call(what = coxKM::coxKM, 
                                          args = sFormals),
                           error = function(e){NULL},
                           warning = function(w){NULL})

    if( is(SKATResult, "NULL") ) {
      #----------------------------------------------------------#
      # If tryCatch returned NULL, remove a column from A        #
      #----------------------------------------------------------#
      res <- shrinkA(A = A, 
                     aEV = aEV, 
                     bEV = bEV, 
                     nZa = nZa, 
                     nZb = nZb, 
                     nX = nX)

      if( is(res,"NULL") ) {
        #------------------------------------------------------#
        # If res is NULL, no stable solution.  Set p.value to  #
        # NULL and exit.                                       #
        #------------------------------------------------------#
        p.value <- NULL
        break
      } else {
        A = res$A
        nZa = res$nZa
        nZb = res$nZb
        rm(res)
        next
      } 
    }

    p.value <- SKATResult$p.value

    break
  }

  return(list( "pv" = p.value,
              "nZa" = nZa,
              "nZb" = nZb))

}

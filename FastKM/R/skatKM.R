#----------------------------------------------------------------------#
# For continuous or dichotomous trait types, iterate over SKAT methods #
# until a stable  solution is obtaind.                                 #
# Return number of eigenvectors kept and the p.values.                 #
#----------------------------------------------------------------------#
skatKM <- function(A, 
                   cKernel,  
                   y,  
                   aEV,  
                   bEV, 
                   opt, ...) {

  #------------------------------------------------------------------#
  # Load SKAT package.                                               #
  #------------------------------------------------------------------#
  hasSKAT <- requireNamespace("SKAT", quietly = TRUE)
  if( !hasSKAT ) {
    stop("R package SKAT is required for continuous trait analyses.",
         call. = FALSE)
  }

  #------------------------------------------------------------------#
  # Match user provided inputs to SKAT_Null_Model with formal        #
  # arguments of SKAT_Null_Model.                                    #
  #------------------------------------------------------------------#
  userArgs <- list(...)

  snmFormals <- as.list(formals(SKAT::SKAT_Null_Model))

  for( i in names(userArgs) ) {

    if( i %in% names(snmFormals) ) snmFormals[[i]] <- userArgs[[i]]

  }

  #------------------------------------------------------------------#
  # Set formal argument "out_type" indicating type of trait          #
  #------------------------------------------------------------------#
  snmFormals[["out_type"]] <- opt

  #------------------------------------------------------------------#
  # Match user provided inputs to SKAT with formal arguments of SKAT.#
  #------------------------------------------------------------------#
  sFormals <- as.list(formals(SKAT::SKAT))

  for( i in names(userArgs) ) {

    if( i %in% names(sFormals) ) sFormals[[i]] <- userArgs[[i]]

  }

  #------------------------------------------------------------------#
  # Set formal argument "kernel" to interaction kernel and           #
  # is_check_genotype to FALSE.                                      #
  #------------------------------------------------------------------#
  sFormals[["kernel"]] <- cKernel
  sFormals[["is_check_genotype"]] <- FALSE
  sFormals[["is_dosage"]] <- FALSE

  nZa <- as.integer(round(length(aEV),0L))
  nZb <- as.integer(round(length(bEV),0L))
  nX <- ncol(A) - nZa - nZb

  cat("Running SKAT.\n")

  while(TRUE) {

    #--------------------------------------------------------------#
    # Attempt to calculate model parameters.                       #
    #--------------------------------------------------------------#
    snmFormals[["formula"]] <- y ~ A
    snmFormals[["data"]] <- data.frame(y,A)

    SKATNull <- tryCatch(expr = do.call(what = SKAT::SKAT_Null_Model,
                                        args = snmFormals),
                         error = function(e){NULL},
                         warning = function(w){NULL})

    if( is(SKATNull, "NULL") ) {
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
        #--------------------------------------------------------#
        # If res is NULL, no stable solution. Set p.value to     #
        # NULL and exit.                                         #
        #--------------------------------------------------------#
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
    sFormals[["Z"]] <- A
    sFormals[["obj"]] <- SKATNull

    SKATResult <- tryCatch(expr = do.call(what = SKAT::SKAT,
                                          args = sFormals),
                           error = function(e){NULL},
                           warning = function(w){NULL})

    if( is(SKATResult, "NULL") ) {
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
        #--------------------------------------------------------#
        # If res is NULL, no stable solution. Set p.value to     #
        # NULL and exit.                                         #
        #--------------------------------------------------------#
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

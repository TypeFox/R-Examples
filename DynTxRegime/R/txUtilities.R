txProcess <- function(txVar, data, fSet){

  txVec <- data[,txVar]

  if( any(is.na(txVec)) ) {
    if( is(fSet, "NULL") ) {
      msg <- "Treatment variable cannot be NA unless fSet is defined."
      UserError("input", msg)
    }
  }

  if( is(txVec, "factor") ) {
    opt <- levels(txVec)
  } else if( is(txVec, "numeric") ) {
    txVec <- as.factor(round(txVec,0L))
    opt <- levels(txVec)
    txVec <- levels(txVec)[txVec]
    if( !isTRUE(all.equal(txVec, data[,txVar])) ) {
      msg <- "Treatment must be of type factor or integer."
      UserError("input, msg")
    }
  } else  if( is(txVec, "integer") ) {
    txVec <- as.factor(txVec)
    opt <- levels(txVec)
    txVec <- levels(txVec)[txVec]
  } else {
    msg <- "Treatment must be of type factor or integer."
    UserError("input, msg")
  }

  if( length(opt) <= 1L ) {
    msg <- paste("Number of treatment options is ", length(opt), 
                 ". Verify txName input.",sep="")
    UserError("input", msg)
  }

  txSubsets <- feasibility(superSet = opt, 
                           fSet = fSet, 
                           txName = txVar,
                           data = data)

  result <- new("TxInfo",
                txName     = txVar,
                superSet   = opt,
                subsets    = txSubsets$subsets,
                ptsSubset  = txSubsets$ptsSubset,
                subsetRule = fSet)


}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# feasibility - Determines subset of treatments available to each patient      #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# superSet : super set of tx options given as vector of integers               #
#                                                                              #
# fSet     : a function defining the feasibility rule (NULL is acceptable)     #
#                                                                              #
# txName   : the column header of data for the treatment variable              #
#                                                                              #
# data     : data.frame; covariates and treatment histories                    #
#                                                                              #
#==============================================================================#
#                                                                              #
#= Returns a list:                                                            =#
#=  subsets   : a list of vectors.                                            =#
#=              Each vector is a treatment subset.                            =#
#=              The name of each element is the subset nickname               =#
#=  ptsSubset : a character vector.                                           =#
#=              The ith element is the nickname of the tx subset available to =#
#=              the ith patient                                               =#
#                                                                              #
#==============================================================================#
feasibility <- function(superSet, 
                        fSet, 
                        txName, 
                        data){

  txVec <- data[,txName]

  levels_in_data <- superSet

  if( !is(fSet, "NULL") ){
    #----------------------------------------------------------------------#
    # Determine what formals are indicated in user specified function.     #
    #----------------------------------------------------------------------#
    formal_arguments <- names(formals(fSet))

    if( length(formal_arguments) == 0L || is(formal_arguments,"NULL")) {
      msg <- "Formal arguments of fSet could not be identified."
      UserError("input", msg)
    }

    #----------------------------------------------------------------------#
    # Identify if formals match the column names of the data set           #
    #----------------------------------------------------------------------#
    matched_colnames <- match(formal_arguments, colnames(data))

    if( any(is.na(matched_colnames)) ){
      #------------------------------------------------------------------#
      # If any formals are not matched to the column headers of data,    #
      # assume that formal argument is data itself. However, if the      #
      # length of matched_colnames is greater than 1, this is not true;  #
      # throw error.                                                     #
      #------------------------------------------------------------------#
      if( length(formal_arguments) > 1L ) {
        msg <- paste("fSet formal arguments ", 
                     paste(formal_arguments[is.na(matched_colnames)],
                           collapse = ", "), 
                     " could not be found in dataset.", sep="")
        UserError("input", msg)
      }

      if( !(formal_arguments %in% c("data")) ) {
        msg <- paste("fSet formal arguments ", 
                     paste(formal_arguments, collapse = ", "), 
                     " could not be found in dataset.", sep="")
        UserError("input", msg)
      }

      #------------------------------------------------------------------#
      # General function that passes the covariates of a single patient  #
      # to fSet to obtain the list of feasible txs for the individual    #
      #------------------------------------------------------------------#
      func1 <- function(x, formal_arguments, fSet){
        argList <- list()
        argList[[ formal_arguments ]] <- x

        txvec <- do.call(what = fSet, args = argList)

        txvec <- verifyTxvec(txvec, levels_in_data)

        return( txvec )
      }

      #------------------------------------------------------------------#
      # For each patient, call fSet to obtain feasible tx set.           #
      #------------------------------------------------------------------#
      ptsSubset <- vector(mode="character", length=nrow(data))
      subsets <- list()

      for( i in 1L:nrow(data) ) {

        tt <- func1(data[i,,drop=FALSE], formal_arguments, fSet)

        subsets[[ tt[[1]] ]] <- tt[[2]]

        ptsSubset[i] <- tt[[1]]

      }

    } else {
      #------------------------------------------------------------------#
      # If all formals are matched to the column names of data           #
      #------------------------------------------------------------------#
      #------------------------------------------------------------------#
      # General function that passes the covariates of a single patient  #
      # to fSet to obtain a vector of feasible txs for the individual    #
      #------------------------------------------------------------------#
      func2 <- function(x, nms, cls, fSet){
        argList <- list()
        argList[nms] <- x[cls]
        txvec <- do.call(what = fSet, args = argList)

        txvec <- verifyTxvec(txvec, levels_in_data)

        return( txvec )
      }

      #------------------------------------------------------------------#
      # For each patient, call fSet to obtain feasible tx set.           #
      #------------------------------------------------------------------#
      ptsSubset <- vector(mode="character", length=nrow(data))
      subsets <- list()

      for( i in 1L:nrow(data) ) {

        tt <- func2(x = data[i,,drop=FALSE], 
                    nms = formal_arguments,  
                    cls = matched_colnames,  
                    fSet = fSet)

        subsets[[ tt[[1]] ]] <- tt[[2]]

        ptsSubset[i] <- tt[[1]]

      }

    }
  } else {
    #----------------------------------------------------------------------#
    # If no feasibility rules are defined, use superSet for each patient   #
    #----------------------------------------------------------------------#
    ptsSubset <- rep("SS", nrow(data))
    subsets <- list("SS" = superSet)
  }

  return( list(  "subsets" = subsets,
               "ptsSubset" = ptsSubset) )

}


verifyTxvec <- function(txvec, levels_in_data) {

  if( !is(txvec, "list") ) {
    msg <- "fSet must return a list."
    UserError("input", msg)
  }

  if( as.integer(round(length(txvec),0L)) != 2L ) {
    msg <- "fSet must return a list of length 2."
    UserError("input", msg)
  }

  if( !is(txvec[[1]], "character") ) {
    msg <- "First element of list returned by fSet must be a character."
    UserError("input", msg)
  }

  if( as.integer(round(length(txvec[[2]]),0L)) == 0L ) {
    msg <- paste("At least one treatment option must be available ",
                 "to every patient --- verify fSet.", sep="")
    UserError("input", msg)
  }

  if( is(txvec[[2]], "factor") ) {
    txvec[[2]] <- levels(txvec[[2]])[txvec[[2]]]
  } else if( is(txvec[[2]], "numeric") ) {
    txvec[[2]] <- as.factor(round(txvec[[2]],0L))
    txvec[[2]] <- levels(txvec[[2]])[txvec[[2]]]
  } else  if( is(txvec[[2]], "integer") ) {
    txvec[[2]] <- as.factor(txvec[[2]])
    txvec[[2]] <- levels(txvec[[2]])[txvec[[2]]]
  } else {
    msg <- "fSet defined treatment must be of type factor or integer"
    UserError("input, msg")
  }

  if( !all(txvec[[2]] %in% c(NA,levels_in_data)) ) {
    msg <- paste("fSet returned treatment options not present in ",
                 "the original data.", paste(txvec[[2]],collapse=","),
                 "{",paste(levels_in_data,collapse=","),"}",sep="")
    UserError("input", msg)
  }

  return(txvec)

}

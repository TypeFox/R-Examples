remoteCalc <- function(input1Name, input2Name = NULL, FUN, outputName, input1Pos = '.GlobalEnv', input2Pos = '.GlobalEnv', outputPos = '.GlobalEnv'){
  FUN <- match.fun(FUN)
  status <- mpi.remote.exec(localCalc, input1Name, input2Name, FUN, outputName,
                  input1Pos, input2Pos, outputPos, ret = TRUE)
  if("try-error" %in% sapply(status, class))
    stop("remoteCalc: error on slaves:\n", status)  
  invisible(NULL)
}

localCalc <- function(input1Name, input2Name, FUN, outputName, input1Pos, input2Pos, outputPos){
  # function that allows generic operations to be performed on objects on the slave processes; e.g., subtracting two vectors or matrices 
  # this function assumes dimensions of the two inputs are consistent - e.g. that both are distributed (and therefore broken up into pieces that lie within a dProblem or both are complete
  if(is.null(input2Name))
    status <- try(assign(outputName, FUN(get(input1Name, pos = eval(as.name(input1Pos)))), pos = eval(as.name(outputPos)))) else
  status <- try(assign(outputName, FUN(get(input1Name, pos = eval(as.name(input1Pos))), get(input2Name, pos = eval(as.name(input2Pos)))),
                       pos = eval(as.name(outputPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}



remoteForwardsolve <- function(cholName, inputName, outputName, cholPos = '.GlobalEnv', inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2 = NULL, h1 = 1, h2 = NULL){
  status <- mpi.remote.exec(localForwardsolve, cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class))
    stop("remoteForwardsolve: error on slaves:\n", status)  
  invisible(NULL)
}

localForwardsolve <- function(cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1, h2){
  status <- try( {
    output <- alloc(inputName, inputPos) # this creates a copy of the input
    if(is.null(n2)){
      .Call("forwardsolve_wrapper", output, get(cholName, eval(as.name(cholPos))), as.integer(n1), as.integer(h1), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    } else{
      .Call("forwardsolve_matrix_wrapper", output, get(cholName, eval(as.name(cholPos))), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    }
    assign(outputName, output, pos = eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

remoteBacksolve <- function(cholName, inputName, outputName, cholPos = '.GlobalEnv', inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2 = NULL, h1 = 1, h2 = NULL){
  if(!is.null(n2)) {
    stop("remoteBacksolve: backsolve into a matrix not currently implemented.")
  } else {
    status <- mpi.remote.exec(localBacksolve, cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
    if("try-error" %in% sapply(status, class))
      stop("remoteBacksolve: error on slaves:\n", status)  
    invisible(NULL)
  }
}

localBacksolve <- function(cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1, h2){
  status <- try( {
    output <- alloc(inputName, inputPos) # this creates a copy of the input
    if(is.null(n2)){
      .Call("backsolve_wrapper", output, get(cholName, eval(as.name(cholPos))), as.integer(n1), as.integer(h1), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    } else{
      stop("localBacksolve: backsolve into a matrix not currently implemented.")
    }
    assign(outputName, output, pos = eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}


remoteCalcChol <- function(matName, cholName, matPos = '.GlobalEnv', cholPos = '.GlobalEnv', n, h = 1){
  out <- mpi.remote.exec(localCalcChol, matName, cholName, matPos, cholPos, n, h, ret = TRUE)
  if("try-error" %in% sapply(out, class)) 
    stop("remoteCalcChol: error on slaves:\n", out)
  if(min(out) != 0 || max(out) != 0)
    stop("remoteCalcChol: Problem with Cholesky decomposition: either input matrix is not positive definite or there is a bug in the linear algebra implementation.")
  invisible(NULL)
}


localCalcChol <- function(matName, cholName, matPos, cholPos, n, h){
  status <- try( {
    L <- alloc(matName, matPos) # creates a copy of the input matrix
    out <- .Call("cholesky_wrapper", as.double(L), as.integer(n), as.integer(h), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    assign(cholName, L, eval(as.name(cholPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(out)
}

remoteMultChol <- function(cholName, inputName, outputName, cholPos = '.GlobalEnv', inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2 = NULL, h1 = 1, h2 = NULL){
  status <- mpi.remote.exec(localMultChol, cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteMultChol: error on slaves:\n", status)
  invisible(NULL)
}
  
localMultChol <-  function(cholName, inputName, outputName, cholPos, inputPos, outputPos, n1, n2, h1 = 1, h2){
  status <- try( {
    output <- alloc(inputName, inputPos) # creates a copy of the input vector
    if(is.null(n2)){
      .Call("mult_chol_vector_wrapper", output, get(cholName, eval(as.name(cholPos))), get(inputName, eval(as.name(inputPos))), as.integer(n1), as.integer(h1), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    } else {
      .Call("mult_chol_matrix_wrapper", output, get(cholName, eval(as.name(cholPos))), get(inputName, eval(as.name(inputPos))), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    }    
    assign(outputName, output, eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}
 
remoteCrossProdMatVec <- function(matName, inputName, outputName, matPos = '.GlobalEnv', inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2, h1 = 1, h2 = 1){
  status <- mpi.remote.exec(localCrossProdMatVec, matName, inputName, outputName, matPos, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteCrossProdMatVec: error on slaves:\n", status)
  invisible(NULL)
}

localCrossProdMatVec <- function(matName, inputName, outputName, matPos, inputPos, outputPos, n1, n2, h1, h2){
  status <- try( {
    output <- alloc(getDistributedVectorLength(n2, h2))
    .Call("mult_cross_wrapper", output, get(matName, eval(as.name(matPos))), get(inputName, eval(as.name(inputPos))), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    assign(outputName, output, eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

remoteCrossProdMatSelf <- function(inputName, outputName, inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2, h1 = 1, h2 = 1){
  status <- mpi.remote.exec(localCrossProdMatSelf, inputName, outputName, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteCrossProdMatSelf: error on slaves:\n", status)
  invisible(NULL)
}

localCrossProdMatSelf <- function(inputName, outputName, inputPos, outputPos, n1, n2, h1 = 1, h2 = 1){
  status <- try( {
    output <- alloc(getDistributedTriangularMatrixLength(n2, h2))
    .Call("cross_prod_self_wrapper", output, get(inputName, eval(as.name(inputPos))), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    assign(outputName, output, eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

remoteCrossProdMatSelfDiag <- function(inputName, outputName, inputPos = '.GlobalEnv', outputPos = '.GlobalEnv', n1, n2, h1 = 1, h2 = 1){
  status <- mpi.remote.exec(localCrossProdMatSelfDiag, inputName, outputName, inputPos, outputPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteCrossProdSelfDiag: error on slaves:\n", status)
  invisible(NULL)
}

localCrossProdMatSelfDiag <- function(inputName, outputName, inputPos, outputPos, n1, n2, h1, h2){
  status <- try( {
    output <- alloc(getDistributedVectorLength(n2, h2))
    .Call("cross_prod_self_diag_wrapper", output, get(inputName, eval(as.name(inputPos))), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    assign(outputName, output, eval(as.name(outputPos)))
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

remoteConstructRnormVector <- function(objName, objPos = '.GlobalEnv', n, h = 1){
  status <- mpi.remote.exec(localConstructRnormMatrix, objName, objPos, n1 = n, n2 = NULL, h1 = h, h2 = NULL, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteConstructRnormMatrix: error on slaves:\n", status)
  invisible(NULL)
}

remoteConstructRnormMatrix <- function(objName, objPos = '.GlobalEnv', n1, n2, h1 = 1, h2 = 1){
  status <- mpi.remote.exec(localConstructRnormMatrix, objName, objPos, n1, n2, h1, h2, ret = TRUE)
  if("try-error" %in% sapply(status, class)) 
    stop("remoteConstructRnormMatrix: error on slaves:\n", status)
  invisible(NULL)
}

localConstructRnormMatrix <- function(objName, objPos, n1, n2, h1, h2){
  if(is.null(n2)){  # vectors are stored differently than one-column matrices
    status <- try( {
      bsr <- (n1 + .bigGP$D*h1 - 1) %/% (.bigGP$D*h1)
      if( .bigGP$I == .bigGP$J ) {
        assign(objName, rnorm(bsr*h1), pos = eval(as.name(objPos)))
      } else {
        assign(objName, numeric(0), pos = eval(as.name(objPos)))
      }
    })
  } else{
    status <- try( {
      bsr <- (n2 + .bigGP$D*h2 - 1) %/% (.bigGP$D*h2)
      bsc <- (n1 + .bigGP$D*h1 - 1) %/% (.bigGP$D*h1)
      if( .bigGP$I == .bigGP$J ) {
        assign(objName, rnorm(bsr*bsc*h1*h2), pos = eval(as.name(objPos)))
      } else {
        assign(objName, rnorm(bsr*bsc*h1*h2*2), pos = eval(as.name(objPos)))
      }
    })
  }
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

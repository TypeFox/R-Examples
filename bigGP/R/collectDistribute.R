push <- function(.tmp, objName = deparse(substitute(.tmp)), objPos = ".GlobalEnv"){
  ' insert doc-string '
  mpi.bcast.Robj2slave(.tmp)
  status <- mpi.remote.exec(localAssign, objName, '.tmp', objPos, ret = TRUE)
  remoteRm(.tmp)
  if("try-error" %in% sapply(status, class))
    stop("push: Push to pos ", objPos, " failed:\n", status)
  invisible(NULL)
}

localAssign <- function(nameToAssign, currentName, objPos = ".GlobalEnv"){
  status <- try(assign(nameToAssign, eval(as.name(currentName)), pos = eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}


pull <- function(objName, objPos = ".GlobalEnv", tag = 1){
  # couldn't figure out mpi.gather.Robj, so not used here
  ' insert doc-string '
  if(!is.character(objName)) stop("pull: 'objName' must be a character string.")
  status <- mpi.remote.exec(localPullTest, objName, objPos, ret = TRUE)
  if("try-error" %in% sapply(status, class)) {
    stop("pull: Pull from pos '", objPos, "' failed:\n", status)
  } else{
    mpi.remote.exec(localPull, objName, objPos, tag = tag, ret = FALSE)
    results <- list()
    for(src in 1:.bigGP$P)
      results[[src]] <- mpi.recv.Robj(src, tag)
    return(results)
  }
}

localPullTest <- function(objName, objPos){
  status <- try(get(objName, pos = eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

localPull <- function(objName, objPos, tag = 1){
  mpi.send.Robj(get(objName, pos = eval(as.name(objPos))), 0, tag)
}

collectDiagonal <- function(objName, objPos = '.GlobalEnv', n, h = 1){
  result <- rep(0, n)
  status <- mpi.remote.exec(localCollectDiagonalTest, objName, objPos, n, h, ret = TRUE)
  if("try-error" %in% sapply(status, class)) {
    stop("collectDiagonal: Collection from pos '", objPos, "' failed:\n", status)
  } else{
    mpi.remote.exec(localCollectDiagonal, objName, objPos, n, h, ret = FALSE)
    .Call("collect_diagonal_wrapper", as.double(result), as.double(0), as.integer(n), as.integer(h), as.integer(0), as.integer(0), .bigGP$D, PACKAGE="bigGP")
    return(result)
  }
}

localCollectDiagonalTest <- function(objName, objPos, n, h){
  status <- try(get(objName, eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}


localCollectDiagonal <- function(objName, objPos, n, h){
  input <- get(objName, eval(as.name(objPos)))
  .Call("collect_diagonal_wrapper", as.double(0), as.double(input), as.integer(n), as.integer(h), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
  invisible(NULL)
}


collectVector <- function(objName, objPos = '.GlobalEnv', n, h = 1){
  result <- rep(0, n)  # allocate space for full vector
  status <- mpi.remote.exec(localCollectVectorTest, objName, objPos, n, h, ret = TRUE) 
  if("try-error" %in% sapply(status, class)) {
    stop("collectVector: Collection from pos '", objPos, "' failed:\n", status)
  } else {
    mpi.remote.exec(localCollectVector, objName, objPos, n, h, ret = FALSE) 
    .Call("collect_vector_wrapper", as.double(result), as.double(0), as.integer(n), as.integer(h), as.integer(0), as.integer(0), .bigGP$D, PACKAGE="bigGP")  # do it on the master; this should put the full vector in 'result' in the correct order 
    return(result)
  }
}

localCollectVectorTest <- function(objName, objPos, n, h){
  status <- try(get(objName, eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

localCollectVector <- function(objName, objPos, n, h){
  input <- get(objName, eval(as.name(objPos)))
  .Call("collect_vector_wrapper", as.double(0), as.double(input), as.integer(n), as.integer(h), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
    invisible(NULL)
}

collectTriangularMatrix <- function(objName, objPos = '.GlobalEnv', n, h = 1){
  result <- rep(0, n*n)
  status <- mpi.remote.exec(localCollectTriangularMatrixTest, objName, objPos, n, h, ret = TRUE) 
  if("try-error" %in% sapply(status, class)) {
    stop("collectTriangularMatrix: Collection from pos '", objPos, "' failed:\n", status)
  } else {
    mpi.remote.exec(localCollectTriangularMatrix, objName, objPos, n, h, ret = FALSE)
    .Call("collect_triangular_matrix_wrapper", as.double(result), as.double(0), as.integer(n), as.integer(h), as.integer(0), as.integer(0), .bigGP$D, PACKAGE="bigGP")
    result <- matrix(result, n, n)
    result <- result + t(result)
    diag(result) <- diag(result)/2
    return(result)
  }
}

localCollectTriangularMatrixTest <- function(objName, objPos, n, h){
  status <- try(get(objName, eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)  
}

localCollectTriangularMatrix <- function(objName, objPos, n, h){
  input <- get(objName, eval(as.name(objPos)))
  .Call("collect_triangular_matrix_wrapper", as.double(0), as.double(input), as.integer(n), as.integer(h), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
  invisible(NULL)
}

collectRectangularMatrix <- function(objName, objPos = '.GlobalEnv', n1, n2, h1 = 1, h2 = 1){
  result <- rep(0, n1*n2)
  status <- mpi.remote.exec(localCollectRectangularMatrixTest, objName, objPos, n1, n2, h1, h2, ret = TRUE) 
  if("try-error" %in% sapply(status, class)) {
    stop("collectRectangularMatrix: Collection from pos '", objPos, "' failed:\n", status)
  } else {
    mpi.remote.exec(localCollectRectangularMatrix, objName, objPos, n1, n2, h1, h2, ret = FALSE) 
    .Call("collect_rectangular_matrix_wrapper", as.double(result), as.double(0), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), as.integer(0), as.integer(0), .bigGP$D, PACKAGE="bigGP")
    return(t(matrix(result, nrow = n2, ncol = n1)))
  }
}

localCollectRectangularMatrixTest <- function(objName, objPos, n1, n2, h1, h2){
  status <- try(get(objName, eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)  
}

localCollectRectangularMatrix <- function(objName, objPos, n1, n2, h1, h2){
  input <- get(objName, eval(as.name(objPos)))
  .Call("collect_rectangular_matrix_wrapper", as.double(0), as.double(input), as.integer(n1), as.integer(n2), as.integer(h1), as.integer(h2), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
  invisible(NULL)
}

distributeVector <- function(obj, objName = deparse(substitute(obj)), objPos = '.GlobalEnv', n, h = 1){
  input <- as.double(obj) # having this here catches cases when 'obj' argument doesn't exist, before going into Rmpi calls that can hang
  status <- mpi.remote.exec(localDistributeVectorTest, objName, objPos, n, h, ret = TRUE)
  if("try-error" %in% sapply(status, class)) {
    stop("distributeVector: Distribution to pos ", objPos, " failed:\n", status)
  } else {
    mpi.remote.exec(localDistributeVector, objName, objPos, n, h, ret = FALSE)
    .Call("distribute_vector_wrapper", as.double(0), as.double(obj), as.integer(n), as.integer(h), as.integer(0), as.integer(0), .bigGP$D, PACKAGE="bigGP")
    invisible(NULL)
  }
}

localDistributeVectorTest <- function(objName, objPos, n, h){
  status <- try(assign(objName, 0, eval(as.name(objPos))))
  if(class(status) == "try-error") invisible(status) else invisible(NULL)  
}

localDistributeVector <- function(objName, objPos, n, h){
  if(.bigGP$I == .bigGP$J) len <- h * ceiling(n / (.bigGP$D * h)) else len <- 0
  output <- rep(0, len)
  .Call("distribute_vector_wrapper", as.double(output), as.double(0), as.integer(n), as.integer(h), .bigGP$I, .bigGP$J, .bigGP$D, PACKAGE="bigGP")
  assign(objName, output, eval(as.name(objPos)))
  invisible(NULL)
}


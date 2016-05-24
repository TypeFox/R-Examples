.print.glmnri <- function(x, ...)
{
  fun_name = if (is.null(attr(x,"function"))) "glm" else attr(x,"function")
  
  mc <- as.character(attr(x, "call"))
  
  cat(paste("Call: ", fun_name, "(", sep = ""))
  cat(paste(mc[2], mc[1]))
  cat(paste("",mc[3]))
  cat(")\n\n")
  
  dim1 <- vector(mode = "numeric", length = 0)
  dim2 <- vector(mode = "numeric", length = 0)
  dim3 <- vector(mode = "numeric", length = 0)
  cat("Models contain following parameters:\n")
  for (i in 1:length(names(x)))
  {
    cat(paste("[[", i, "]] ", names(x)[i], "\n", sep=""))
    dim1 <- c(dim1, dim(x[[i]])[1])
    dim2 <- c(dim2, dim(x[[i]])[2])
    dim3 <- c(dim3, dim(x[[i]])[3])
  }
  
  if (any(c(any(dim1 != dim1[1]),any(dim2 != dim2[1]),any(dim2 != dim2[1]))))
  {
    cat("Dimensions of parameter:\n")
    for (i in 1:length(names(x)))
      cat(paste("  ", names(x)[i], ": ", dim1[i], ", ", dim2[i], ", ", dim3[i], "\n", sep = ""))
      
  } else {
    cat(paste("Dimension of each parameter: ", dim(x[[1]])[1], ", ",
              dim(x[[1]])[2], ", ", dim(x[[1]])[3], "\n", sep=""))
  }
  if (!is.null(dimnames(x[[1]])[[1]]))
  {
    cat("  Term labels\n")
    nam <- dimnames(x[[1]])[[1]]
    for (i in 1:length(nam))
    {
      cat(paste("    (", i, ") ", nam[i], "\n", sep=""))
    }
  }
  if (!is.null(dimnames(x[[1]])[[2]]))
  {
    cat("  Variable names\n")
    nam <- dimnames(x[[1]])[[2]]
    for (i in 1:length(nam))
    {
      cat(paste("    (", i, ") ", nam[i], "\n", sep=""))
    }
  }
}



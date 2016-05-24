getFiniteNri <- function(x)
{
  if (class(x) != "Nri")
    stop("'x' must be of class Nri")
  
  
  result <- apply(x$nri, 1, FUN = function(i) return(all(is.finite(i))))
  
  
  result <- distMat3D(result*1, ncol = ncol(x@nri), nlyr = 1)
  
  result <- as.matrix(result)

  wl1 <- matrix(data=rep.int(c(1:dim(result)[1]),dim(result)[1]),
                nrow=dim(result)[1], ncol=dim(result)[1], byrow = FALSE)
  wl2 <- matrix(data=rep.int(c(1:dim(result)[1]),dim(result)[1]),
                nrow=dim(result)[1], ncol=dim(result)[1], byrow = TRUE)



  relevant <-  data.frame(dim1=as.vector(wl1[result==1 & is.finite(result) & lower.tri(result)]),
                          dim2=as.vector(wl2[result==1 & is.finite(result) & lower.tri(result)]))
  relevant <- data.frame(Band_1=x$wavelength[relevant$dim1], Band_2=x$wavelength[relevant$dim2])
  return(list(Indices=relevant, Models=NULL))
}
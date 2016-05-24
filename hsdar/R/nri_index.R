setMethod("[", "Nri",
          function(x, i, ...)
{  
  x@nri <- distMat3D(x@nri[,,i], lower_tri = TRUE)
  attribute(x) <- attribute(x)[i,]
  return(x)  
})

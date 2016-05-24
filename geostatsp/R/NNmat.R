nbMat = rbind(
    c(1,0,0),# self
    c(2,0,1), #up
    c(2,0,-1), #down
    c(2,1,0), #right,
    c(2,-1,0), #left
    c(3,-1,-1),
    c(3,-1,1),
    c(3,1,-1),
    c(3,1,1),
    c(4,0,2), #up
    c(4,0,-2), #down
    c(4,2,0), #right,
    c(4,-2,0), #left
    c(5,-1,-2),
    c(5,-1,2),
    c(5,1,-2),
    c(5,1,2),
    c(5,-2,-1),
    c(5,-2,1),
    c(5,2,-1),
    c(5,2,1),
    c(6,0,3), #up
    c(6,0,-3), #down
    c(6,3,0), #right,
    c(6,-3,0) #left
)


NNmat = function(N,Ny=N,nearest=3) {
  UseMethod("NNmat")	
}

NNmat.Raster = function(N,Ny=N,nearest=3) {
  res = NNmat(ncol(N),nrow(N),nearest)
  
  attributes(res)$raster= raster(N)
  
  res
}	

NNmat.default = function(N, Ny=N,nearest=3) {
  
  if(nearest<3){
    nbMat = nbMat[nbMat[,1]<=nearest,] 
  }
  
  Nx = N
 

  Ncell = Nx * Ny
    
    xMat = values(raster(matrix(1:Nx, byrow=TRUE, ncol=Nx, nrow=Ny)))
    yMat = values(raster(matrix(1:Ny, byrow=FALSE, ncol=Nx, nrow=Ny)))
    
    
    xNb = outer(nbMat[,2],xMat, FUN='+')   
    yNb = outer(nbMat[,3],yMat, FUN='+')
    nbCode = matrix(nbMat[,1], nrow(xNb), ncol(xNb))
    
    xNb[xNb<1] = NA
    xNb[xNb>Nx] = NA
    yNb[yNb<1] = NA
    yNb[yNb>Ny] = NA
    
    nbIndex = xNb + Nx*(yNb-1)
    nbCol = matrix(1:ncol(xNb), nrow(xNb), ncol(xNb), byrow=TRUE)
    notNa = !(is.na(xNb) | is.na(yNb))
    
    result = sparseMatrix(
        i=nbIndex[notNa],
        j=nbCol[notNa], 
        x=nbCode[notNa]
    )
    result=forceSymmetric(result)
  
  attributes(result)$Nx = Nx
  attributes(result)$Ny = Ny
  
  
  return(result)
}
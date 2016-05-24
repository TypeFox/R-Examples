figures8landm <- function(figure,data){
  
 if(figure == "cube"){
  cube <- data
  cube_arr <- array(as.matrix(cube), dim = c(dim(cube)[1], 3, 1))
   
  x <- cube_arr[,,1] ; type = "p" ; color = 2   
  xt <- array(0, c(dim(x), 1))
  xt[, , 1] <- x
  x <- xt
  
  rgl.open()
  rgl.bg(color = "white")
  k <- dim(x)[1]
  sz <- centroid.size(x[, , 1])/sqrt(k)/30
  joinline = c(4,6,5,2,4,3,8,6,8,7,5,7,1,3,1,2) 
  plotshapes3d(x, type = type, color = color, size = sz, joinline = joinline)
  axes3d(color = "black")
  title3d(xlab = "x", ylab = "y", zlab = "z", color = "black")
  
  rgl.texts(x = x[,,1][,1][1] + 1, y = x[,,1][,2][1] + 0.5, z = x[,,1][,3][1] + 0.5, "1", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][2] + 0, y = x[,,1][,2][2] + 0.5, z = x[,,1][,3][2] + 0.5, "2", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][3] + 1.2, y = x[,,1][,2][3] + 0.5, z = x[,,1][,3][3] + 0.5, "3", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][4] + 1, y = x[,,1][,2][4] + 0.5, z = x[,,1][,3][4] + 0.5, "4", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][5] + 0.8, y = x[,,1][,2][5] + 0.5, z = x[,,1][,3][5] + 0.5, "5", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][6] + 1, y = x[,,1][,2][6] + 0.5, z = x[,,1][,3][6] + 0.5, "6", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][7] + 1, y = x[,,1][,2][7] + 0.5, z = x[,,1][,3][7] + 0.5, "7", adj = c(1,1), cex = 1, col = "black")
  rgl.texts(x = x[,,1][,1][8] + 0, y = x[,,1][,2][8] + 0.5, z = x[,,1][,3][8] + 0.5, "8", adj = c(1,1), cex = 1, col = "black")
 }else if(figure == "paral"){  
   paral <- data
   paral_arr <- array(as.matrix(paral), dim = c(dim(paral)[1], 3, 1))
  
   x <- paral_arr[,,1] ; type = "p" ; color = 2  
   xt <- array(0, c(dim(x), 1))
   xt[, , 1] <- x
   x <- xt
  
   rgl.open()
   rgl.bg(color = "white")
   k <- dim(x)[1]
   sz <- centroid.size(x[, , 1])/sqrt(k)/30
   joinline = c(8,6,5,7,8,3,4,2,5,2,1,7,1,3,4,6)
   plotshapes3d(x, type = type, color = color, size = sz, joinline = joinline)
   axes3d(color = "black")
   title3d(xlab = "x", ylab = "y", zlab = "z", color = "black")
    
   rgl.texts(x = x[,,1][,1][1] + 2, y = x[,,1][,2][1] + 1, z = x[,,1][,3][1] + 1, "1", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][2] + 1, y = x[,,1][,2][2] + 1, z = x[,,1][,3][2] + 1, "2", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][3] + 2, y = x[,,1][,2][3] + 1, z = x[,,1][,3][3] + 1, "3", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][4] + 2, y = x[,,1][,2][4] + 1, z = x[,,1][,3][4] + 1, "4", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][5] + 1, y = x[,,1][,2][5] + 1, z = x[,,1][,3][5] + 1, "5", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][6] + 2, y = x[,,1][,2][6] + 1, z = x[,,1][,3][6] + 1, "6", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][7] + 2, y = x[,,1][,2][7] + 1, z = x[,,1][,3][7] + 1, "7", adj = c(1,1), cex = 1, col = "black")
   rgl.texts(x = x[,,1][,1][8] + 2, y = x[,,1][,2][8] + 1, z = x[,,1][,3][8] + 1, "8", adj = c(1,1), cex = 1, col = "black")
  }else{
    stop("Sorry,that figure does not belong to the package")  
   }   
}
# Author: Carlos A. Silva. Uses code by Remko Duursma in Maeswrap package, see Plotstand.
interpol<- function(input,col) {
  surf.3d <- t(convhulln(input,options = "QJ")) 
  rgl.triangles(input[surf.3d,1],input[surf.3d,2],input[surf.3d,3],col=col,alpha = c(1.0),
                lit = TRUE,ambient = "black",specular = "white",emission = "black",shininess = 50.0,
                smooth = TRUE, texture = NULL,front = "fill",back ="fill",fog = F) 
}

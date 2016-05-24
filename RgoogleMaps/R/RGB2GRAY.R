`RGB2GRAY` <-structure(function#translates an RGB image matrix to gray scale
 ###  This function translates the rgb values of the array myTile into a scalar matrix with just one gray value per pixel.
 (
   myTile ##<< rgb image matrix, usually array with 3 dimensions
 ){
   #Gray scale intensity = 0.30R + 0.59G + 0.11B
   stopifnot(attr(myTile, "type") != "gray");
   f = NULL
   ##details<< Gray scale intensity defined as  0.30R + 0.59G + 0.11B
   if (class(myTile)[1] == 'nativeRaster'){
     # comparing nativeRaster values is not easy, so let's do write/read again
     f = file.path(tempdir(), "tmpMap.png")
     writePNG(myTile, f)
     myTile = readPNG(f)
   }
   
   if (class(myTile)[1] == 'array'){
     #       if (0) {
     #         tmp = matrix(0.7, ncol=dim(myTile)[2], nrow=dim(myTile)[1] )
     #         for (i in 1:nrow(tmp))
     # 		for (j in 1:ncol(tmp))
     # 			tmp[i,j] = .3*myTile[i,j,1] + .59*myTile[i,j,2] + .11*myTile[i,j,3]
     #          myTile = tmp
     #        }
     #or in vectorized form:
     ncol=dim(myTile)[2]; nrow=dim(myTile)[1]
     dim(myTile) <- c(prod(dim(myTile)[1:2]),dim(myTile)[3])
     tmp = .3*myTile[,1] + .59*myTile[,2] + .11*myTile[,3]
     myTile = matrix(tmp, ncol=ncol, nrow=nrow )
     
   } 
   if (!is.null(f)){
     writePNG(myTile, f)
     myTile = readPNG(f, native = TRUE)
   }
   return(myTile);	
   ### image tile
}, ex=function(){
  if (interactive()){
    BrooklynLatLon = getGeoCode("Brooklyn")
    mapBrooklyn <- GetMap(center=BrooklynLatLon, destfile = file.path(tempdir(), "Brooklyn.png"), 
                   zoom=11, size = c(240,240))
    mapBrooklynBW$myTile = RGB2GRAY(mapBrooklyn$myTile)
    PlotOnStaticMap(mapBrooklynBW)
  }
})


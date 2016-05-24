`pltomo` <-
function(x,y,MOD,i, colmap=rainbow(100), zlim=NULL, bkgr="DarkSlateGray4", ...)
  {
    
    if(missing(colmap)) { colmap=tomo.colors(100) }
    if(missing(bkgr)) { bkgr="DarkSlateGray4" }
    if(missing(zlim)) { zlim=NULL }

    Z = MOD[[i]]

    
    plot(range(x), range(y), asp=1,   axes=FALSE, ann=FALSE,  type='n', ...)

    
    rect(min(x), min(y), max(x), max(y), col=bkgr)

    if(is.null(Z)) {
      print(paste(sep=" ", "This LAYER is NULL:", i) )
      return(-1)
    }
    
    if(all(is.na(Z))) {
      print(paste(sep=" ", "This LAYER is all NA: ", i) )
      return(-1)
    }

  ###  print("plotting")
    
    if(is.null(zlim))
      {
        image(x,y,Z, asp=1, col=colmap,  xlab="km", ylab="km", add=TRUE)
      }
    else
      {

        if(length(zlim)==2)
          {
        image(x,y,Z, asp=1, col=colmap,  xlab="km", ylab="km", zlim=zlim, add=TRUE)
      }
        else
          {

            image(x,y,Z, asp=1, col=colmap,  xlab="km", ylab="km", breaks=zlim, add=TRUE)
          }
      }

    return(1)
  }


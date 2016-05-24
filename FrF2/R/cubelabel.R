cubelabel <- function(cdobj, lablev=list(c(-1,1),c(-1,1),c(-1,1)),labs=NULL,
          pos="outab", offs=1, adj=c(0,0),col="blue",cex=1.5){
    ### function to write labels onto the corners of the cube
    ### lablev is a list of factor levels that is automatically expanded to the right order
    ###       for the eight cube points; lablev is ignored, if labs is given
    ### labs is a vector of eight labels in the order of the points
    ###       for labelling only part of the points, give empty labels for unlabelled points

    ### pos can be NULL, 1(below),2(left),3(above),4(right),
    ###       and special words outab (outer above/below), in ab (inner above/below)
    ###                         outlr (outer left/right),  in lr (inner left/right)
    ### for labelling the corners inside plot symbols, use pos=NULL

    if(!is.null(pos)){
      if (pos=="outab") pos <- c(1,1,1,1,3,3,3,3)
      else if (pos=="inab") pos <- c(3,3,3,3,1,1,1,1)
      else if (pos=="outlr") pos <- c(2,4,2,4,2,4,2,4)
      else if (pos=="inlr") pos <- c(4,2,4,2,4,2,4,2)
    }
    x <- cdobj$res3d$xyz.convert(cdobj$cub)$x
    y <- cdobj$res3d$xyz.convert(cdobj$cub)$y
    if (is.null(labs)){
        lc <- expand.grid(lablev)  ## labels for corners
        if (all(sapply(lc,is.numeric)))
           labs <- paste(lc[,1],lc[,2], lc[,3],sep="/")
        else
           labs <- paste(as.character(lc[,1]),as.character(lc[,2]),
               as.character(lc[,3]),sep="")
    }
    text(cdobj$res3d$xyz.convert(cdobj$cub), labels=labs, col=col,
         cex=cex, pos=pos, offset=offs, adj=adj)
}

spodtSpatialLines <- function(object,data){
  x<- coordinates(data)[,1]
  y<- coordinates(data)[,2]
  if (class(object@racine)=="f.spodt"){#root is a leaf
    d2<-matrix(ncol=4,nrow=4)
    if (!is.null(x)  |  !is.null(y)){
        minX <- min(x) - (max(x) - min(x))/50
        maxX <- max(x) + (max(x) - min(x))/50
        minY <- min(y) - (max(y) - min(y))/50
        maxY <- max(y) + (max(y) - min(y))/50
        d2 <- rbind(c(minX,minY,minX,maxY), c(minX,maxY,maxX,maxY), c(maxX,maxY,maxX,minY), c(maxX,minY,minX,minY))
    }
    colnames(d2)[c(3,4)] <- colnames(d2)[c(1,2)]
    rownames(d2)<-1:dim(d2)[1]

  }
  else{
    if((is.na(object@cl.grf)[1])){  #no graft
        data.coord <- recursif.carte.sp(object@racine)
        ind.na <- which(is.na(data.coord[,1]))
        d2<- data.coord[-ind.na,]
        if (!is.null(x)  |  !is.null(y)){
            minX <- min(x) - (max(x) - min(x))/50
            maxX <- max(x) + (max(x) - min(x))/50
            minY <- min(y) - (max(y) - min(y))/50
            maxY <- max(y) + (max(y) - min(y))/50
            d2 <- rbind(d2, c(minX,minY,minX,maxY), c(minX,maxY,maxX,maxY), c(maxX,maxY,maxX,minY), c(maxX,minY,minX,minY))
        }
        colnames(d2)[c(3,4)] <- colnames(d2)[c(1,2)]

        points.d2 <- rbind(d2[,1:2],d2[,3:4])
        points.d2 <- unique(points.d2)
        nbLig <- dim(points.d2)[1]
        nbLig.d2 <- dim(d2)[1]
        d2<-recurs.d2(points.d2,d2)
        rownames(d2)<-1:dim(d2)[1]

      }
      else{ #graft
        d2 <- object@brd[,1:4]
        points.d2 <- rbind(d2[,1:2],d2[,3:4])
        points.d2 <- unique(points.d2)
        nbLig <- dim(points.d2)[1]
        nbLig.d2 <- dim(d2)[1]
        nd2<-recurs.d2(points.d2,d2)
        rownames(nd2)<-1:dim(nd2)[1]
        d1 <- object@sgmts.grf
        d2 <- rm.dupl(d1,nd2)
        d1_inv <- cbind(d1[,3:4], d1[,1:2])
        d2 <- rm.dupl(d1_inv,d2)

      }
   }


    spodtL <- list()

    spodtL <- sapply(1:(dim(d2)[1]), function(i){
                       x<-rbind(d2[i,1:2],d2[i,3:4])
                       spodtL[[i]] <- Lines(list(Line(x)),ID=rownames(d2)[i])})
    spodtSL<-SpatialLines(spodtL)
     proj4string(spodtSL)<-proj4string(data)
    return(spodtSL)

}
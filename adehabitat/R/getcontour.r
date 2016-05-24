"getcontour" <- function(x)
  {
      ## Verifications
      if (!inherits(x, "asc"))
          stop("should be an object of class asc")

      xyc<-getXYcoords(x)

      ## The function rajfond adds an empty line of 0 at the top
      ## and the bottom of the map, and an empty column of 0 at
      ## the left and the right of an map

      rajfond<-function(x)
      {
          nr<-nrow(x)
          nc<-ncol(x)

          f<-rep(0,nr)
          x<-cbind(f,x,f)
          f<-rep(0,nc+2)
          x<-rbind(f,x,f)
      }

      ## transform the map into a matrix of 0 (missing values) and 1 (mapped)
      x[!is.na(x)]<-1
      x[is.na(x)]<-0

      ## adds empty lines and columns at the border of the map with rajfond
      x<-rajfond(x)

      ## sequential labelling of the connex components thanks to a call
      ## to the C function "seqeticorr"
      toto<-.C("seqeticorr", as.double(t(x)), as.integer(nrow(x)),
               as.integer(ncol(x)), PACKAGE="adehabitat")

      ## gets the resulting map
      etiquete<-matrix(toto[[1]], nrow=nrow(x), byrow=TRUE)

      ## removes the additionnal lines and columns of 0s
      etiquete<-etiquete[-c(1,nrow(etiquete)),-c(1,ncol(etiquete))]


      ###########################
      ## Now, use this map to get the contour of the connex components


      ## Bases
      entree<-list()
      sorties<-c(0, 0, 0)
      lev<-levels(factor(toto[[1]]))
      lev<-lev[lev!="0"]

      ## for each connex component
      for (i in lev) {

          ## prepares the data
          j<-as.numeric(i)
          tmp<-etiquete
          tmp[tmp!=j]<-0
          tmp[tmp==j]<-1

          ## adds empty lines and columns
          tmp<-rajfond(tmp)

          ## Components of 3 pixels min are needed
          if (sum(as.vector(tmp))>=3)
          {
              ## computes the number of vertices of the connex components
              toto<-.C("lcontour", as.double(t(tmp)), as.integer(nrow(tmp)),
                       as.integer(ncol(tmp)),  as.integer(0),
                       PACKAGE="adehabitat")[[4]]

              ## computes the connex components
              pol<-.C("getcontour", as.double(t(tmp)), as.integer(nrow(tmp)),
                      as.integer(ncol(tmp)), integer(toto), integer(toto),
                      as.integer(toto), PACKAGE="adehabitat")

              ## output
              xt<-c(0,xyc$x,0)
              yt<-c(0,xyc$y,0)
              x<-xt[pol[[4]]]
              y<-yt[pol[[5]]]
              sorties<-rbind(sorties, cbind(rep(j,length(x)), x, y))
          }
      }
      ## output as area
      sorties<-sorties[-1,]
      row.names(sorties)<-1:nrow(sorties)
      sorties<-as.data.frame(sorties)
      sorties[,1]<-factor(sorties[,1])
      names(sorties)<-c("id","x", "y")
      sorties<-as.area(sorties)
      return(sorties)
  }


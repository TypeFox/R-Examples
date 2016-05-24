"mahasuhab" <- function(x, pts, type=c("distance", "probability"))
  {
      ## Verifications
      type<-match.arg(type)
      if (!inherits(x, "SpatialPixelsDataFrame"))
          stop("should be an object of class SpatialPixelsDataFrame")
      gridded(x) <- TRUE
      gr <- gridparameters(x)
      if (nrow(gr) > 2)
          stop("x should be defined in two dimensions")
      if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
          stop("the cellsize should be the same in x and y directions")
      if (!inherits(pts, "SpatialPoints"))
          stop("should inherit from class \"SpatialPoints\"")

      ## Computation of the variance-covariance matrix of the used points:
      hihi<-join(pts, x)
      hihi <- hihi[!is.na(hihi[,1]),]
      used<-list()

      ## factors are transformed into dummy variables
      for (i in 1:ncol(hihi)) {
          if (is.factor(hihi[,i]))
              used[[i]]<-acm.disjonctif(data.frame(hihi[,i]))[,-1]
          else
              used[[i]]<-hihi[,i]
      }
      used[[i+1]]<-rep(1, nrow(hihi))
      hihi<-as.data.frame(used)
      hihi<-hihi[!is.na(hihi[,1]),]
      mu<-apply(hihi,2, function(x) mean(x, na.rm=TRUE))
      varcov<-t(as.matrix(hihi))%*%as.matrix(hihi)/nrow(hihi)


      ## habitat Availability
      kasc <- slot(x,"data")
      ava<-list()

      ## factors are transformed into dummy variables
      for (i in 1:ncol(kasc)) {
          if (is.factor(kasc[,i]))
              ava[[i]]<-acm.disjonctif(data.frame(kasc[,i]))[,-1]
          else
              ava[[i]]<-kasc[,i]
      }

      ava[[i+1]]<-rep(1, nrow(kasc))
      df<-as.data.frame(ava)

      ## computation of the Mahalanobis distances
      map<-mahalanobis(as.matrix(df), mu, varcov)
      if (type=="probability")
          map<-1-pchisq(map, ncol(hihi)-1)
      map <- data.frame(MD=map)
      coordinates(map) <- coordinates(x)
      gridded(map) <- TRUE
      return(map)
  }


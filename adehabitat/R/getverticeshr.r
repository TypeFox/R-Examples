"getverticeshr" <- function(x, lev=95)
  {
      ## Verifications
      if ((!inherits(x,"khr")))
          stop("non convenient data-type")
      if (inherits(x,"khrud"))
          x<-getvolumeUD(x)
      if (inherits(x,"kbbhrud"))
          x<-getvolumeUD(x)
      if (length(lev)>1)
          stop("lev should be of length 1")
      ## output list
      contour<-list()

      ## for each animal
      for (i in 1:length(x)) {

          ## gets the UD and keep areas upper than lev
          ud<-x[[i]]$UD

          ## gets the contour of the connected features
          xyl <- getXYcoords(ud)
          re <- contourLines(x = xyl$x,
                             y = xyl$y,
                             ud, nlevels = 1,
                             levels = lev)
          so <- do.call("rbind", lapply(1:length(re), function(i) {
              so <- data.frame(fac=rep(i, length(re[[i]]$x)),
                               x=re[[i]]$x,
                               y=re[[i]]$y)
              return(so)
          }))
          so[,1] <- factor(so[,1])

          contour[[i]]<-as.area(so)
      }
      ## output of class "kver"
      names(contour)<-names(x)
      class(contour) <- "kver"
      return(contour)
  }


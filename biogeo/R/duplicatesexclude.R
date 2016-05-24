duplicatesexclude <-
function (dat,res=10) 
  {
    fieldsmissing(dat, fields = c("ID", "x", "y", "Species", "Exclude","Reason"))
    rst <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = res/60)
    rst[] <- 1:ncell(rst)
    spp <- as.character(unique(dat$Species))
    xy <- data.frame(coord2numeric(dat$x), coord2numeric(dat$y))
    cid <- cellFromXY(rst, xy)
    #dat<-data.frame(dat,cid)
    #e <- extract(rst, xy)
    for (i in 1:length(spp)) {
      f1 <- which(as.character(dat$Species) == spp[i])
      #e1 <- e[f1]
      xf <- (duplicated(cid[f1])) * 1
      fx<-which(xf==1)
      f2<-f1[fx]
      dat$Exclude[f2] <- 1
      dat$Reason <- as.character(dat$Reason)
      dat$Reason[f2] <- 'Duplicated'
    }
    return(dat)
  }

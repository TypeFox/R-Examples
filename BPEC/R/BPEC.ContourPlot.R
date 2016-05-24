BPEC.ContourPlot <- function(MCMCout,CoordsLocs,GoogleEarth=0,colorcode=c(7,5,6,3,2,8,4,9,10))
{
  writeLines("Creating geographical contour plot...")
  
  Rootlocs = numeric(3)
  Rootlocs[1] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-1+1)[2]
  Rootlocs[2] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-2+1)[2]
  Rootlocs[3] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-3+1)[2]
  
  CoordsLocs=CoordsLocs[,1:2]
  NCP = 15
  NCPR = 2 * NCP + 1
  
  MeanSamples = MCMCout$SampleMeansR
  CovSamples = MCMCout$SampleCovsR
  NoClusters=dim(MeanSamples)[2]
  FullMeanSamples = MeanSamples
  FullCovSamples = CovSamples
  
  SubSeq=seq(1,length(CovSamples[1,1,1,]),length.out=10)
  MeanSamples = MeanSamples[,,SubSeq]
  CovSamples =CovSamples[,,,SubSeq]
  
  fullclust=numeric(NoClusters)
  for(i in 1:NoClusters)
  {
    if(length(which(!is.na(MeanSamples[1,i,])))>0)
      #                  if(length(which(!is.na(MCMCout$SampleMeansR[1,i,])))>length(MCMCout$SampleMeansR[1,i,])/2)
    {
      fullclust[i]=1
    }
    else
    {
      fullclust[i]=0
    }
  }
  
 
  xtot=seq(min(min(MeanSamples[2,fullclust==1,]-sqrt(CovSamples[2,2,fullclust==1,]),na.rm=TRUE)),max(max(MeanSamples[2,fullclust==1,]+sqrt(CovSamples[2,2,fullclust==1,]),na.rm=TRUE)),length.out=NCPR)    
  ytot=seq(min(min(MeanSamples[1,fullclust==1,]-sqrt(CovSamples[1,1,fullclust==1,]),na.rm=TRUE)),max(max(MeanSamples[1,fullclust==1,]+sqrt(CovSamples[1,1,fullclust==1,]),na.rm=TRUE)),length.out=NCPR)    
  
  plot(xtot,ytot, col='white',xlab='',ylab='',xlim=range(xtot),ylim=range(ytot),axes=FALSE,asp=1)
  
  nv = length(CoordsLocs[, 1])
  newCoordsLocs = array(NA, c(nv, 2))
  for(i in 1:nv) 
  {
    newCoordsLocs[i, 1] = CoordsLocs[i, 2] 
    newCoordsLocs[i, 2] = CoordsLocs[i, 1]
  }
  
  AT=TRUE
  for(i in 1:NoClusters)  
  {     
    for(it in 1:dim(MeanSamples)[3])
    {
      newmeans = array(NA, c(2, NoClusters))
      newcovs = array(NA, c(2, 2, NoClusters))
      for(j in 1:NoClusters) 
      {
        newmeans[1, j] = MeanSamples[2, j,it]
        newmeans[2, j] =  MeanSamples[1, j,it]
        newcovs[1, 1, j] = CovSamples[2, 2, j,it]
        newcovs[1, 2, j] =   CovSamples[2, 1, j,it]
        newcovs[2, 1, j] =   CovSamples[1, 2, j,it]
        newcovs[2, 2, j] = CovSamples[1, 1, j,it]
      }
      z = array(NA, c(NCPR, NCPR, NoClusters))
      if(is.na(newcovs[1,1,i])==0)
      {
        x=seq(min(min(MeanSamples[2,i,],na.rm=TRUE))-2*sqrt(max(max(CovSamples[2,2,i,],na.rm=TRUE))),max(max(MeanSamples[2,i,],na.rm=TRUE))+2*sqrt(max(max(CovSamples[2,2,i,],na.rm=TRUE))),length.out=NCPR)        
        y=seq(min(min(MeanSamples[1,i,],na.rm=TRUE))-2*sqrt(max(max(CovSamples[1,1,i,],na.rm=TRUE))),max(max(MeanSamples[1,i,],na.rm=TRUE))+2*sqrt(max(max(CovSamples[1,1,i,],na.rm=TRUE))),length.out=NCPR)
        for(k in 1:NCPR) 
        {
          for(j in 1:NCPR) 
          {                     
            z[k, j, i] = dmvnorm(c(x[k],y[j]),newmeans[,i], newcovs[,,i])            
          }
        }
        ColCont=col2rgb(colorcode[i], alpha = FALSE)
        
        Conti=contourLines(x,y,z[,,i], levels=0.5*max(z[,,i]))
        polygon(Conti[[1]]$x,Conti[[1]]$y,col=c(rgb(ColCont[1]/255,ColCont[2]/255,ColCont[3]/255,alpha=0.1)), border=NA,new =AT,lwd=2, xlab = "", ylab = "")
        AT=FALSE
      }
    }
    #cat ("Press [enter] to continue")
    #line <- readline()
  }
  for(i in 1:NoClusters) 
  {    
    newmeans[1,i] = mean(MeanSamples[2, i,],na.rm=TRUE)
    newmeans[2,i]= mean(MeanSamples[1, i,],na.rm=TRUE)
    
    newcovs[1, 1, i] = mean(CovSamples[2, 2, i,],na.rm=TRUE)
    newcovs[1, 2, i] = mean(CovSamples[2, 1, i,],na.rm=TRUE)
    newcovs[2, 1, i] = mean(CovSamples[1, 2, i,],na.rm=TRUE)
    newcovs[2, 2, i] = mean(CovSamples[1, 1, i,],na.rm=TRUE)        
  }
  x = array(NA, c(NCPR, NoClusters))
  y = array(NA, c(NCPR, NoClusters))
  z = array(NA, c(NCPR, NCPR, NoClusters))
  
  for(i in 1:NoClusters)
  {
    if(fullclust[i]>0)
    {
      x[, i] = newmeans[1, i] + (( - NCP:NCP) * 6 * sqrt(newcovs[1, 1,i]))/NCPR
      y[, i] = newmeans[2, i] + (( - NCP:NCP) * 6 * sqrt(newcovs[2, 2,i]))/NCPR
      
      for(k in 1:NCPR) 
      {
        for(j in 1:NCPR) 
        {
          z[k, j, i] = dmvnorm(c(x[k, i], y[j, i]),newmeans[, i], newcovs[,  , i])
        }
      }
    }
  }
  
  for(i in 1:NoClusters) 
  {
    if(fullclust[i]>0)
    {
      contour(x[,i],y[,i],z[,,i],col=colorcode[i],levels=0.5*max(z[,,i]), xlim = c(min(xtot,na.rm=TRUE),max(xtot,na.rm=TRUE)),ylim =c(min(ytot,na.rm=TRUE),max(ytot,na.rm=TRUE)),add = TRUE, lwd=1,drawlabels=FALSE,axes = FALSE, xlab = "", ylab = "")                          
    }
  }           
  
  
  
  
  #  points(newCoordsLocs, pch = 16, cex = 0.3)
  points(newCoordsLocs, pch = 17, cex = 0.5, col=1)
  
  points(newCoordsLocs[Rootlocs[1], 1], newCoordsLocs[Rootlocs[1], 2], pch = 17, cex = 1.1, col = 1)
  points(newCoordsLocs[Rootlocs[2], 1], newCoordsLocs[Rootlocs[2], 2], pch = 17, cex = 0.9, col = 1)
  points(newCoordsLocs[Rootlocs[3], 1], newCoordsLocs[Rootlocs[3], 2], pch = 17, cex = 0.7, col = 1)
  
  world(add=TRUE)
  
  if(GoogleEarth==1)
  {
    # second section copied from contrgen, slightly modified
    # initialize 
    reslist = list()
    for(i in 1:NoClusters) 
    {
      if(i > 1) 
      {
        AT = TRUE
      }
      reslist[[i]] = list(x[, i], y[, i], z[, , i]); 
    }
    # end of sections copied from contrgen
    
    #### this section produces Google Earth output
    ##### FUNCTION for kmlPoints output ######
    #[R-sig-Geo] Preparing KML files for GoogleMaps/Earth:
    #Why no kmlPoint() in maptools package?
    #Robert J. Hijmans r.hijmans at gmail.com
    kmlPoints = function(obj, kmlFile, kmlname = "") {
      if (kmlname=="") {kmlname = basename(kmlFile)}
      kmlHeader = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<kml xmlns=\"http://earth.google.com/kml/2.2\">", "<Document>",paste("<name>", kmlname, "</name>", sep = ""))
      kml = ""
      obj[,3] = gsub('&', 'and', obj[,3])
      for (p in 1:length(obj[,1])) {
        #		kml = append(kml, "<Folder>")
        kml = append(kml, "<Placemark>")
        kml = append(kml, paste("<name>", obj[p,3], "</name>", sep = ""))
        kml = append(kml, paste("<Point><coordinates>", obj[p,1], ',',obj[p,2], ",0</coordinates></Point>", sep = ""))
        kml = append(kml, "</Placemark>")
        #		kml = append(kml, "</Folder>")
      }
      kmlFooter = c("</Document>", "</kml>")
      cat(paste(c(kmlHeader, kml, kmlFooter), sep = "", collapse ="\n"), "\n", file = kmlFile, sep = "")
    }
    
    #############################################################
    ##### MODIFIED FUNCTION for kmlRoots output ######
    root.color = "ff00ff00"
    kmlRoots = function(obj2, kmlFile, kmlname = "")
    {
      if (kmlname=="") {kmlname = basename(kmlFile)}
      kmlHeader = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<kml xmlns=\"http://earth.google.com/kml/2.2\">", "<Document>",paste("<name>", kmlname, "</name>", sep = ""))
      kml = ""
      obj[,3] = gsub('&', 'and', obj[,3])
      for (p in 1:length(obj2[,1])) {
        #		kml = append(kml, "<Folder>")
        kml = append(kml, "<Placemark>")
        kml = append(kml, paste("<name>", obj2[p,3], "</name>", sep = ""))
        kml = append(kml, "<Style>")
        kml = append(kml, "<IconStyle>")
        kml = append(kml, paste("<color>",root.color,"</color>", sep = ""))
        kml = append(kml, paste("<scale>", 1.2, "</scale>", sep = ""))
        kml = append(kml, "<Icon>")
        kml = append(kml, paste("<href>","http://maps.google.com/mapfiles/kml/shapes/arrow.png","</href>", sep = ""))
        kml = append(kml, "</Icon>")
        kml = append(kml, "</IconStyle>")
        kml = append(kml, "</Style>")
        kml = append(kml, paste("<Point><coordinates>", obj2[p,1], ',',obj2[p,2], ",0</coordinates></Point>", sep = ""))
        kml = append(kml, "</Placemark>")
        #		kml = append(kml, "</Folder>")
      }
      kmlFooter = c("</Document>", "</kml>")
      cat(paste(c(kmlHeader, kml, kmlFooter), sep = "", 
                collapse ="\n"), "\n", file = kmlFile, sep = "")
    }
    
    
    ##### export to GoogleEarth GE, load library with drivers:
    # - for Windows: give temporary files the extension .kml, 
    #then double-click the *.kml files in the temp directory to start GE
    
    #### GENERATION OF PLOTS
    ## loop over lists     
    for(i in 1:NoClusters) 
    {
      if(fullclust[i]>0)
      {
        x = unlist(reslist[[i]][[1]]);
        #x
        y = unlist(reslist[[i]][[2]]);
        #y
        z = matrix(unlist(reslist[[i]][[3]]), nrow=length(x), ncol=length(y));
        # z
        #generate contourLines
        l = contourLines(x,y,z); 
        ## set geographic projection
        ## here only geographic coordinates are allowed 
        llCRS = CRS("+proj=longlat +datum=WGS84 +towgs84=0,0,0")
        ################################################################################
        # points
        newCoordsLocs.sp = SpatialPoints(newCoordsLocs, llCRS)
        newCoordsLocs.spdf  =  SpatialPointsDataFrame(newCoordsLocs, proj4string=llCRS,data.frame(matrix(seq(1:length(newCoordsLocs[,1])),ncol=1))) ##### export sampling points to GoogleEarth
        #  str(newCoordsLocs.spdf)
        # newCoordsLocs.spdf@data
        obj = cbind(coordinates(newCoordsLocs.sp),newCoordsLocs.spdf@data)
        # obj
        #########
        # contours
        ## create empty list to hold lines
        ll = vector("list", length(l))
        ## loop over list of contours to create Lines object
        for (k in 1:length(l))
          ll[[k]] = Lines(list(Line(cbind(l[[k]]$x, l[[k]]$y))), as.character(k))
        #### make those lines Spatial, and generate a SpatialLinesDataFrame
        ll = SpatialLines(ll, proj4string=llCRS)
        #### return the length of the lines slot, how many lines are there?
        out  =  sapply(slot(ll, "lines"),  function(x) {kmlLine(x,name=slot(x, "ID"), col=colorcode[i], lwd=1.5,description=paste("contours_[i]",  slot(x, "ID"))) })
        ##### generates R temporary files
        #tf  =  tempfile()
        tf  =  paste("GoogleEarthContour",as.character(i),".kml",sep = "")
        kmlFile  =  file(tf,  "w")
        cat(kmlLine(kmlname=paste("Contours",i,sep="_"),  kmldescription="<i>BPC</i>")$header,file=kmlFile,  sep="\n")
        cat(unlist(out["style",]),  file=kmlFile,  sep="\n")
        cat(unlist(out["content",]),  file=kmlFile,  sep="\n")
        cat(kmlLine()$footer,  file=kmlFile,  sep="\n")
        close(kmlFile)
      }
    }
    #############
    ### plot points
    #tf  =  tempfile()
    tf  =  "GoogleEarthPoints.kml"
    kmlFile  =  file(tf,  "w")
    kmlPoints(obj, tf)
    
    llCRS = CRS("+proj=longlat +datum=WGS84 +towgs84=0,0,0")
    newCoordsLocs.sp = SpatialPoints(newCoordsLocs[Rootlocs,], llCRS)
    newCoordsLocs.spdf  =  SpatialPointsDataFrame(newCoordsLocs[Rootlocs,], proj4string=llCRS,data.frame(matrix(seq(1:length(newCoordsLocs[Rootlocs,1])),ncol=1))) ##### export sampling points to GoogleEarth
    #  str(newCoordsLocs.spdf)
    # newCoordsLocs.spdf@data
    obj2 = cbind(coordinates(newCoordsLocs.sp),newCoordsLocs.spdf@data)
    close(kmlFile)
    tf  =  "GoogleEarthRoots.kml"
    kmlFile  =  file(tf,  "w")
    kmlRoots(obj2, tf)
    close(kmlFile)
    #############
  }
}

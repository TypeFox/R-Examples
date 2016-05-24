#'Individual tree LiDAR metrics
#'
#'@description This funtion computes LiDAR-based statistic metrics for individual trees detected from the 3D LiDAR point cloud. 
#'
#'@usage rMetrics(xyziId)
#'
#'@param xyziId A 5-column matrix with the x, y, z coordinates, intensity and the tree id classification for the LiDAR point cloud.
#'@return Returns A matrix of the LiDAR-based metrics for the individual tree detected.
#'@author Carlos Alberto Silva
#'@details
#' 
#'# List of the individual tree LiDAR metrics:
#'
#'# Total number of returns (TotalReturns);   
#'# UTM Easting coordinate of the tree top (Etop);
#'# UTM Northing coordinate of the tree top (Ntop);
#'# Minimum UTM Easting coordinate (Emin);
#'# Minimum UTM Northing coordinate (Nmin);
#'# Maximum UTM Easting coordinate (Emax);
#'# Maxmium UTM Northing coordinate (Nmax);
#'# Tree crown width 01 (Ewidth);
#'# Tree crown width 02 (Nwidth);
#'# Maximum Height (hmax);
#'# Mean height (hmean);
#'# Standard deviation of height (hsd);
#'# Coefficient of variation of height (hcv);
#'# Mode of height (hmod);
#'# 5th percentile of height (h5);	
#'# 10th percentile of height (h10);	
#'# 20th percentile of height (h20);
#'# 25th percentile of height (h25);	
#'# 30th percentile of height (h30);	
#'# 40th percentile of height (h40);	
#'# 50th percentile of height (h50);	
#'# 60th percentile of height (h60);
#'# 70th percentile of height (h70);
#'# 75th percentile of height (h75);
#'# 80th percentile of height (h80);
#'# 90th percentile of height (h90);
#'# 95th percentile of height (h95);
#'# 99th percentile of height (h99);
#'# Maximum intensity (imax);
#'# Mean intensity (imean);
#'# Standard deviation of intensity (isd);
#'# Coefficient of variation of intensity (icv);
#'# Mode of intensity (imod);
#'# 5th percentile of intensity (i5);
#'# 10th percentile of intensity (i10);
#'# 20th percentile of intensity (i20);
#'# 25th percentile of intensity (i25);
#'# 30th percentile of intensity (i30);
#'# 40th percentile of intensity (i40);
#'# 50th percentile of intensity (i50);
#'# 60th percentile of intensity (i60);
#'# 70th percentile of intensity (i70);
#'# 75th percentile of intensity (i75);
#'# 80th percentile of intensity (i80);
#'# 90th percentile of intensity (i90);
#'# 95th percentile of intensity (i95);
#'# 99th percentile of intensity (i99).
#'
#'@examples
#'
#'#=======================================================================#
#'# Individual tree detection using K-means cluster
#'#=======================================================================#
#'# Importing LAS file:
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Reading LAS file
#'LAS<-readLAS(LASfile,short=TRUE)
#'
#'# Setting the xyz coordinates and subsetting the data
#'xyzi<-subset(LAS[,1:4],LAS[,3] >= 1.37)
#'
#'# Finding clusters (trees)
#'clLAS<-kmeans(xyzi[,1:2], 32)
#'
#'# Set the tree id vector
#'Id<-as.factor(clLAS$cluster)
#'
#'# Combining xyzi and tree id 
#'xyziId<-cbind(xyzi,Id)
#'
#'#=======================================================================#
#'#  Computing individual tree LiDAR metrics 
#'#=======================================================================#
#'
#'TreesMetrics<-rMetrics(xyziId)
#'head(TreesMetrics)
#'
#'@export
rMetrics<-function(xyziId) {  
  
  # ----from moments package: Lukasz Komsta et al.(2015) ---#
  "skewness" <-
    function (x, na.rm = FALSE) 
    {
      if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
      else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)] 
        n <- length(x)
        (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
      }
      else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
      else skewness(as.vector(x), na.rm = na.rm)
    }
  
  "kurtosis" <-
    function (x, na.rm = FALSE) 
    {
      if (is.matrix(x)) 
        apply(x, 2, kurtosis, na.rm = na.rm)
      else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)] 
        n <- length(x)
        n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2)
      }
      else if (is.data.frame(x)) 
        sapply(x, kurtosis, na.rm = na.rm)
      else kurtosis(as.vector(x), na.rm = na.rm)
    }
  #-----------------------------------------------------------#
  MetricsList<-matrix(,ncol=68)[-1,]
  nlevels<-as.numeric(levels(factor(xyziId[,5])))
  
  for ( i in nlevels){
    #print(i)
    cat (".");flush.console()
    
    xyz.c<-subset(xyziId[,1:4],xyziId[,5]==i)
    
    if (nrow(xyz.c) <= 1) { 
      xRange<-round(range(xyz.c[,1]), digits=2)
      yRange<-round(range(xyz.c[,2]), digits=2)
      MaxZ<-max(xyz.c[,3])  # fild the max point
      XY<-as.data.frame(subset(xyz.c,xyz.c[,3]==MaxZ)) # get the x and y from the max point
      maxPoint<-round(XY[1,1:2],digits=2)
      
      Metrics<-c(
        npoits<-round(nrow(xyz.c), digits=2),
        maxPoint,
        xRangeMin<-xRange[1],
        xRangeMax<-xRange[2],
        yRangeMin<-yRange[1],
        yRangeMax<-yRange[2],
        xWidth<-round(xRangeMax-xRangeMin,digits=2),
        yWidth<-round(yRangeMax-yRangeMin,digits=2),
        rep(0,61))
      MetricsList<-rbind(MetricsList,c(i,Metrics))
    } else {
      
      MaxZ<-max(xyz.c[,3])  # fild the max point
      XY<-as.data.frame(subset(xyz.c,xyz.c[,3]==MaxZ)) # get the x and y from the max point
      maxPoint<-round(XY[1,1:2],digits=2)
      xRange<-round(range(xyz.c[,1]), digits=2)
      yRange<-round(range(xyz.c[,2]), digits=2)
      
      Metrics<-c( 
        
        # Number of points
        npoits<-round(nrow(xyz.c), digits=2),
        maxPoint,
        # Range UTM E,N
        xRangeMin<-xRange[1],
        xRangeMax<-xRange[2],
        yRangeMin<-yRange[1],
        yRangeMax<-yRange[2],
        xWidth<-round(xRangeMax-xRangeMin,digits=2),
        yWidth<-round(yRangeMax-yRangeMin,digits=2),
        
        # hieght metrics
        hmax=round(max(xyz.c[,3]), digits=2),
        hmin=round(min(xyz.c[,3]), digits=2),
        hmean=round(mean(xyz.c[,3]),digits=2),
        hmedian=round(median(xyz.c[,3]),digits=2),
        hmode = round(as.numeric(names(table(xyz.c[,3]))[which.max(table(xyz.c[,3]))]), digits=2),
        hvar=round(var(xyz.c[,3]),digits=2),
        hsd=round(sd(xyz.c[,3]),digits=2),
        hcv=round((sd(xyz.c[,3])/mean(xyz.c[,3]))*100,digits=2),
        hkurtosis=round(kurtosis(xyz.c[,3]),digits=2),
        hskewness=round(skewness(xyz.c[,3]),digits=2),
        h5=round(quantile(xyz.c[,3],0.05),digits=2),
        h10=round(quantile(xyz.c[,3],0.1),digits=2),
        h15=round(quantile(xyz.c[,3],0.15),digits=2),
        h20=round(quantile(xyz.c[,3],0.20),digits=2),
        h25=round(quantile(xyz.c[,3],0.25),digits=2),
        h30=round(quantile(xyz.c[,3],0.30),digits=2),
        h35=round(quantile(xyz.c[,3],0.35),digits=2),
        h40=round(quantile(xyz.c[,3],0.40),digits=2),
        h45=round(quantile(xyz.c[,3],0.45),digits=2),
        h50=round(quantile(xyz.c[,3],0.50),digits=2),
        h55=round(quantile(xyz.c[,3],0.55),digits=2),
        h60=round(quantile(xyz.c[,3],0.60),digits=2),
        h65=round(quantile(xyz.c[,3],0.65),digits=2),
        h70=round(quantile(xyz.c[,3],0.70),digits=2),
        h75=round(quantile(xyz.c[,3],0.75),digits=2),
        h80=round(quantile(xyz.c[,3],0.85),digits=2),
        h90=round(quantile(xyz.c[,3],0.90),digits=2),
        h95=round(quantile(xyz.c[,3],0.95),digits=2),
        h99=round(quantile(xyz.c[,3],0.99),digits=2),
        imax=round(max(xyz.c[,4]), digits=2),
        imin=round(min(xyz.c[,4]), digits=2),
        imean=round(mean(xyz.c[,4]),digits=2),
        imedian=round(median(xyz.c[,4]),digits=2),
        imode = round(as.numeric(names(table(xyz.c[,4]))[which.max(table(xyz.c[,4]))]), digits=2),
        ivar=round(var(xyz.c[,4]),digits=2),
        isd=round(sd(xyz.c[,4]),digits=2),
        icv=round((sd(xyz.c[,4])/mean(xyz.c[,4]))*100,digits=2),
        ikurtosis=round(kurtosis(xyz.c[,4]),digits=2),
        iskewness=round(skewness(xyz.c[,4]),digits=2),
        i5=round(quantile(xyz.c[,4],0.05),digits=2),
        i10=round(quantile(xyz.c[,4],0.1),digits=2),
        i15=round(quantile(xyz.c[,4],0.15),digits=2),
        i20=round(quantile(xyz.c[,4],0.20),digits=2),
        i25=round(quantile(xyz.c[,4],0.25),digits=2),
        i30=round(quantile(xyz.c[,4],0.30),digits=2),
        i35=round(quantile(xyz.c[,4],0.35),digits=2),
        i40=round(quantile(xyz.c[,4],0.40),digits=2),
        i45=round(quantile(xyz.c[,4],0.45),digits=2),
        i50=round(quantile(xyz.c[,4],0.50),digits=2),
        i55=round(quantile(xyz.c[,4],0.55),digits=2),
        i60=round(quantile(xyz.c[,4],0.60),digits=2),
        i65=round(quantile(xyz.c[,4],0.65),digits=2),
        i70=round(quantile(xyz.c[,4],0.70),digits=2),
        i75=round(quantile(xyz.c[,4],0.75),digits=2),
        i80=round(quantile(xyz.c[,4],0.85),digits=2),
        i90=round(quantile(xyz.c[,4],0.90),digits=2),
        i95=round(quantile(xyz.c[,4],0.95),digits=2),
        i99=round(quantile(xyz.c[,4],0.99),digits=2))
      
      MetricsList<-rbind(MetricsList,c(i,Metrics))
    }
  }
  
  colnames(MetricsList)<-c("Tree","TotalReturns","Etop","Ntop","Emin","Nmin","Emax","Nmax","Ewidth","Nwidth","hmax","hmin","hmean","hmedian","hmode",
                           "hvar","hsd","hcv","hkurtosis","hskewness","h5","h10","h15","h20","h25","h30","h35","h40",
                           "h45","h50","h55","h60","h65","h70","h75","h80","h90","h95","h99","imax","imin","imean","imedian","imode",
                           "ivar","isd","icv","ikurtosis","iskewness","i5","i10","i15","i20","i25","i30","i35","i40",
                           "i45","i50","i55","i60","i65","i70","i75","i80","i90","i95","i99")
  return(data.frame(MetricsList))
}

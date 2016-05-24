boxplot_geog<-function(area,data,fence,groups){
  if (missing(fence)) fence=1
  if (missing(groups)) groups=1
  nbox<-nlevels(groups)
  plot(x=c(1:(2.5*nbox-1)),type='n',ylab="data",ylim=c(min(data),max(data)),main="Boxplot",xlab="Geographic Boxplot")
  #xaxt='n'
  for (i in 1:nbox){
    Summary<-summary_geog(area[as.numeric(groups)==i],data[as.numeric(groups)==i],fence)$Summary
    Outliers<-summary_geog(area[as.numeric(groups)==i],data[as.numeric(groups)==i],fence)$Outliers
    Spatial_weighted_mean<-summary_geog(area[as.numeric(groups)==i],data[as.numeric(groups)==i],fence)$Spatial_weighted_mean
    y<-Summary;outliers<-Outliers;ywmean<-Spatial_weighted_mean;x<-c(1:(2*nbox));ix<-2*i
    if(length(outliers)==0){
      outliers=1
      xout<-kronecker(matrix(1,length(outliers)),2*i)
      lines(xout,outliers) 
    }else
    {
      xout<-kronecker(matrix(1,length(outliers)),2*i)
      points(xout,outliers,cex=1.5)
    }
    xinc<-(x[2]-x[1])/20
    xbox<-(x[2]-x[1])/5
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[1],y[1]),col='blue')
    lines(c(x[ix],x[ix],x[ix]-xbox,x[ix]-xbox,x[ix]+xbox,x[ix]+xbox,x[ix]),c(y[1],y[2],y[2],y[4],y[4],y[2],y[2]))
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[5],y[5]),col='blue')
    lines(c(x[ix],x[ix]),c(y[4],y[5]))
    lines(c(x[ix]-xbox,x[ix]+xbox),c(y[3],y[3]),col='blue')
    lines(x[ix],ywmean,pch=16,type='p')
  }
}

summary_geog<-function(area,data,fence){
  if (missing(fence)) fence=1
  isort<-order(data)
  datasort<-data[order(data)]
  csarea=cumsum(area[isort])/sum(area)
  csareapct<-csarea*100
  
  pcts<-c(25,50,75)
  qnts<-pcts/100
  ywmean<-sum(area*data)/sum(area)
  
  y<-c(NA,0,0,0,NA)
  i25<-which(csareapct>25)
  y[2]<-datasort[i25[1]]
  i50<-which(csareapct>50)
  y[3]<-datasort[i50[1]]
  i75<-which(csareapct>75)
  y[4]<-datasort[i75[1]]
  
  iqr<-y[4]-y[2]
  innerfence<-y[2]-1.5*iqr
  outerfence<-y[4]+1.5*iqr
  outliers<- data[data<innerfence | data>outerfence]
  
  if (fence!=1){
    y[1]<-min(datasort)
    y[5]<-max(datasort)}else{
      # outliers
      if (innerfence<min(datasort)){
        y[1]<-min(datasort)} else {
          y[1]<-min(datasort[datasort>=innerfence])}
      if (is.na(y[1])){
        y[1]<-y[2]}
      
      if (outerfence>max(datasort)){
        y[5]<-max(datasort)} else {
          y[5]<-max(datasort[datasort<=outerfence])}
      if (is.na(y[5])){
        y[5] = y[4]}
    }
  names<-c("Min.","1st Quan.","Med.","3rd Quan.","Max.")
  names(y)<-names
  return(list(Summary=c(y),Outliers=c(outliers),Spatial_weighted_mean=c(ywmean)))
}


boxplot_geog_example<-function(Area,Popdents,State){
    Area<-Area;Popdents<-Popdents;State<-State
    Summary<-summary_geog(Area,Popdents,1)$Summary
    Outliers<-summary_geog(Area,Popdents,1)$Outliers
    Spatial_weighted_mean<-summary_geog(Area,Popdents,1)$Spatial_weighted_mean
    y<-Summary;outliers<-Outliers;ywmean<-Spatial_weighted_mean
    xout<-kronecker(matrix(1,length(outliers)),4)
    iter<-length(Outliers)
    plot(xout,outliers,type='n',ylim=c(0,max(outliers)+50),xaxt='n',xlab="Geographic Boxplot  and   Traditional",main='Boxplot of Population Density for 2000 by State for Lower 48',ylab='Population Density (people/km^2)')
    xout<-kronecker(matrix(1,length(outliers)),3)
    points(xout,outliers,cex=1.5,col=c(1:iter))
    x<-c(1:3);ix<-3
    xinc<-(x[2]-x[1])/20
    xbox<-(x[2]-x[1])/5
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[1],y[1]),col='blue')
    lines(c(x[ix],x[ix],x[ix]-xbox,x[ix]-xbox,x[ix]+xbox,x[ix]+xbox,x[ix]),c(y[1],y[2],y[2],y[4],y[4],y[2],y[2]))
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[5],y[5]),col='blue')
    lines(c(x[ix],x[ix]),c(y[4],y[5]))
    lines(c(x[ix]-xbox,x[ix]+xbox),c(y[3],y[3]),col='blue')
    lines(x[ix],ywmean,pch=16,type='p')
    
    for(i in 1:iter){
        arrows(3.2, Outliers[i], 3.1, Outliers[i], length = 0.1, angle = 23.5,
        code = 2, col = i)
        text(3.4,Outliers[i],State[which(Popdents==outliers[i])],cex=0.7,col=i)
    }
    
    summary<-summary_geog(array(1,length(Popdents)),Popdents,1)$Summary
    outliers<-summary_geog(array(1,length(Popdents)),Popdents,1)$Outliers
    mean<-summary_geog(array(1,length(Popdents)),Popdents,1)$Spatial_weighted_mean
    y<-summary;outliers<-outliers;ywmean<-mean
    iters<-length(outliers)
    xout<-kronecker(matrix(1,length(outliers)),5)
    points(xout,outliers,cex=1.5,col=c(1:iters))
    x<-c(1:5);ix<-5
    xinc<-(x[2]-x[1])/20
    xbox<-(x[2]-x[1])/5
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[1],y[1]),col='blue')
    lines(c(x[ix],x[ix],x[ix]-xbox,x[ix]-xbox,x[ix]+xbox,x[ix]+xbox,x[ix]),c(y[1],y[2],y[2],y[4],y[4],y[2],y[2]))
    lines(c(x[ix]-xinc,x[ix]+xinc),c(y[5],y[5]),col='blue')
    lines(c(x[ix],x[ix]),c(y[4],y[5]))
    lines(c(x[ix]-xbox,x[ix]+xbox),c(y[3],y[3]),col='blue')
    lines(x[ix],ywmean,pch=16,type='p')
    
    for(i in 1:iters){
        arrows(4.8, outliers[i], 4.9, outliers[i], length = 0.1, angle = 23.5,
        code = 2, col = i)
        text(4.6,outliers[i],State[which(Popdents==outliers[i])],cex=0.7,col=i)
    }
}

genhistogram<-function(x,breaks=7,input="openair",pollutant=NA,delta=0,plot=TRUE,distr="norm",gap=0.05,columns=2,col="#A52375",emboss=0,xlab="Concentration",ylab="Number of samples",main=NA) {
  
  ## Selection of relevant compound(s).
  # If x is defined as genasis-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="genasis") {
    compounds<-unique(x[,2])
  }
  
  # If x is defined as openair-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="openair") {
    compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind","note")))]
  }
  
  # If x is defined as numeric vector, set compounds=NA.
  if (class(x)!="data.frame") {
    compounds<-NA
    if(!is.na(pollutant)&pollutant!="") {
      compounds<-pollutant[1]
    }
  }
  
  if (length(pollutant)==1) {
    if(is.na(pollutant)|pollutant=="") {
      pollutant<-compounds
    }
  }
  comp<-pollutant[which(is.element(pollutant,compounds))] # Compounds, which will be used.
  if (length(pollutant)>length(comp)) {
    warning(paste0("One or more pollutants (",pollutant[which(!is.element(pollutant,comp))],") was not recognized."))
  }
  
  # Settings of plotting output.
  if (plot==TRUE&length(comp)>1) {
    par(mfrow=c(ceiling(length(comp)/columns),columns))
  }
  if (plot==TRUE&length(comp)<=1) {
    par(mfrow=c(1,1))
  }
  if (length(comp)>6) {
    warning(paste0("There is a lot of compounds to plot in the dataset x (",length(comp),"). The plot could be crowded on some devices."))
  }
  
  # Output variables.
  gborders<-as.data.frame(matrix(NA,0,breaks+1))
  gdistrib<-as.data.frame(matrix(NA,0,breaks))
  
  ## Cycle over compound(s).
  for (compound in comp) {
    
    ## Preparation of the data for the plot.
    # If x is defined as genasis-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="genasis") {     
      valu<-x[which(x[,2]==compound),1]
    }
    
    # If x is defined as openair-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="openair") {
      valu<-x[,compound]
    }
    
    # If x is defined as numeric vector, assignes the variables.
    if (class(x)!="data.frame") {
      valu<-x
    } 
    
    
    # Main headers.    
    if (is.na(main)) {
      plotmain<-as.character(compound)
    } else {
      if (class(x)=="data.frame") {
        plotmain<-paste0(main," (",as.character(compound),")")
      } else
        plotmain<-main     
    }
    
    ## Plotting.
    cells<-breaks+1
    borders=c(0:cells)*((ceiling(((ceiling((1+delta)*(max(valu,na.rm=TRUE)-delta*min(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))-(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1))))/(cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2))))*cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2)))/cells)+(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))
    if (length(valu)<2) {
      borders<-c(0:(cells)*2*valu/cells)
    }
    
    distrib<-c()
    for (i in 1:cells) {
      distrib<-c(distrib,length(valu[which(valu[]>=(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))+((ceiling(((ceiling((1+delta)*(max(valu,na.rm=TRUE)-delta*min(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))-(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1))))/(cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2))))*cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2)))/cells)*(i-1) & valu[]<((floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))+((ceiling(((ceiling((1+delta)*(max(valu,na.rm=TRUE)-delta*min(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))-(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1))))/(cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2))))*cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2)))/cells)*i))]))
    }
    distrib[cells]<-distrib[cells]+length(valu[which(valu[]==((floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))+((ceiling(((ceiling((1+delta)*(max(valu,na.rm=TRUE)-delta*min(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))-(floor((1+delta)*(min(valu,na.rm=TRUE)-delta*max(valu,na.rm=TRUE))/(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1)))*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-1))))/(cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2))))*cells*(10^(floor(log10(max(valu,na.rm=TRUE)-min(valu,na.rm=TRUE)))-2)))/cells)*i))])
    if (length(valu)<2) {
      distrib[length(which(borders<valu))+1]<-1
    }
    
    if (plot==TRUE) {
      column<-(borders[length(borders)]-borders[1])/cells
      column.y<-max(distrib)/cells*1.5
      
      if (distr=="") {
        maximum<-max(distrib)
      } else {
        ddistr<-paste0("d",distr)
        nodes<-seq(from=min(borders)-delta*(borders[length(borders)]-borders[1]),to=max(borders)+delta*(borders[length(borders)]-borders[1]),by=(max(borders)-min(borders))/100)
        
        if (distr=="lnorm") {
          meanx<-mean(log(valu),na.rm=TRUE)
          sdx<-sd(log(valu),na.rm=TRUE)
        } else {
          meanx<-mean(valu,na.rm=TRUE)
          sdx<-sd(valu,na.rm=TRUE)
        }
        maximum<-max(c(distrib,unlist(lapply(nodes,FUN=ddistr,meanx,sdx))*length(valu)*(borders[length(borders)]-borders[1])/cells),na.rm=TRUE)
      }
      
      plot(borders[-length(borders)],distrib,xlim=c(borders[1]-delta*(borders[length(borders)]-borders[1]),borders[length(borders)]+delta*(borders[length(borders)]-borders[1])),ylim=c(0,maximum),cex=0,xlab=xlab,ylab=ylab,xaxt="n",main=plotmain)
      axis(1,labels=borders,at=borders)
      for (i in 1:cells) {
        ld<-c(borders[i]+column*gap,0)
        rd<-c(borders[i+1]-column*gap,0)
        ru<-c(borders[i+1]-column*gap,distrib[i])
        lu<-c(borders[i]+column*gap,distrib[i])
        rect(ld[1],ld[2],ru[1],ru[2],col=col)
        rx<-column*(1-gap)/(3*emboss)
        ry<-column.y*(1-gap)/(3*emboss)
        if (emboss==1) {
          if (distrib[i]>0.23*column.y) {
            polygon(c(borders[i]+column*gap+column*0.08,
                      borders[i]+column*gap+column*0.08-rx*(cos(1*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(2*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(3*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(4*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(4*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(3*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(2*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(1*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15,
                      borders[i+1]-column*gap-column*0.15,
                      borders[i+1]-column*gap-column*0.15+rx*(cos(1*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(2*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(3*pi/8)-1),
                      borders[i+1]-column*gap-column*0.15+rx*(cos(4*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(4*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(3*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(2*pi/8)-1),
                      borders[i]+column*gap+column*0.08-rx*(cos(1*pi/8)-1),
                      borders[i]+column*gap+column*0.08),
                    c(0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(1*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(2*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(3*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(4*pi/8)-1),
                      0+column.y*0.08,
                      0+column.y*0.08,
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(4*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(3*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(2*pi/8)-1),
                      0+column.y*0.08-min(ry,(distrib[i]-0.23*column.y)/2)*(sin(1*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(1*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(2*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(3*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(4*pi/8)-1),
                      distrib[i]-column.y*0.15,
                      distrib[i]-column.y*0.15,
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(4*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(3*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(2*pi/8)-1),
                      distrib[i]-column.y*0.15+min(ry,(distrib[i]-0.23*column.y)/2)*(sin(1*pi/8)-1)),
                    col=rgb(1,1,1,alpha=0.4),border=NA)
          }
        }
        if (emboss==2) {
          ldr<-c(borders[i]+column*gap+rx,min(ry,distrib[i]/2))
          rdr<-c(borders[i+1]-column*gap-rx,min(ry,distrib[i]/2))
          rur<-c(borders[i+1]-column*gap-rx,max(distrib[i]-ry,distrib[i]/2))
          lur<-c(borders[i]+column*gap+rx,max(distrib[i]-ry,distrib[i]/2))
          
          polygon(c(lur[1],rur[1],ru[1],lu[1]),c(lur[2],rur[2],ru[2],lu[2]),col=rgb(1,1,1,alpha=0.4),border=NA)
          polygon(c(rdr[1],rd[1],ru[1],rur[1]),c(rdr[2],rd[2],ru[2],rur[2]),col=rgb(0,0,0,alpha=0.4),border=NA)
          polygon(c(ld[1],rd[1],rdr[1],ldr[1]),c(ld[2],rd[2],rdr[2],ldr[2]),col=rgb(0,0,0,alpha=0.4),border=NA)
          polygon(c(ld[1],ldr[1],lur[1],lu[1]),c(ld[2],ldr[2],lur[2],lu[2]),col=rgb(1,1,1,alpha=0.4),border=NA)
        }
        if (emboss==3) {
          ldr<-c(borders[i]+column*gap+rx,min(ry,distrib[i]/2))
          rdr<-c(borders[i+1]-column*gap-rx,min(ry,distrib[i]/2))
          rur<-c(borders[i+1]-column*gap-rx,max(distrib[i]-ry,distrib[i]/2))
          lur<-c(borders[i]+column*gap+rx,max(distrib[i]-ry,distrib[i]/2))
          
          polygon(c(lur[1],rur[1],ru[1],lu[1]),c(lur[2],rur[2],ru[2],lu[2]),col=rgb(1,1,1,alpha=0.4),border=NA)
          polygon(c(rdr[1],rd[1],ru[1],rur[1]),c(rdr[2],rd[2],ru[2],rur[2]),col=rgb(0,0,0,alpha=0.4),border=NA)
          polygon(c(ld[1],rd[1],rdr[1],ldr[1]),c(ld[2],rd[2],rdr[2],ldr[2]),col=rgb(0,0,0,alpha=0.4),border=NA)
          polygon(c(ld[1],ldr[1],lur[1],lu[1]),c(ld[2],ldr[2],lur[2],lu[2]),col=rgb(1,1,1,alpha=0.4),border=NA)
          
          if (distrib[i]>0.60*column.y) {
            polygon(c(borders[i]+column*gap+column*0.24,
                      borders[i]+column*gap+column*0.24-rx*(cos(1*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(2*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(3*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(4*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(4*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(3*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(2*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(1*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36,
                      borders[i+1]-column*gap-column*0.36,
                      borders[i+1]-column*gap-column*0.36+rx*(cos(1*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(2*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(3*pi/8)-1),
                      borders[i+1]-column*gap-column*0.36+rx*(cos(4*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(4*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(3*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(2*pi/8)-1),
                      borders[i]+column*gap+column*0.24-rx*(cos(1*pi/8)-1),
                      borders[i]+column*gap+column*0.24),
                    c(0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(1*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(2*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(3*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(4*pi/8)-1),
                      0+column.y*0.24,
                      0+column.y*0.24,
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(4*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(3*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(2*pi/8)-1),
                      0+column.y*0.24-min(ry,(distrib[i]-0.60*column.y)/2)*(sin(1*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(1*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(2*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(3*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(4*pi/8)-1),
                      distrib[i]-column.y*0.36,
                      distrib[i]-column.y*0.36,
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(4*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(3*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(2*pi/8)-1),
                      distrib[i]-column.y*0.36+min(ry,(distrib[i]-0.60*column.y)/2)*(sin(1*pi/8)-1)),
                    col=rgb(1,1,1,alpha=0.4),border=NA)
          }
        }
        rect(borders[i]+column*gap,0,borders[i+1]-column*gap,distrib[i],col=NA)
      }
      
      if (distr!="") {
        ddistr<-paste0("d",distr)
        nodes<-seq(from=min(borders)-delta*(borders[length(borders)]-borders[1]),to=max(borders)+delta*(borders[length(borders)]-borders[1]),by=(max(borders)-min(borders))/100)
        lines(nodes,unlist(lapply(nodes,FUN=ddistr,meanx,sdx))*length(valu)*(borders[length(borders)]-borders[1])/cells,lwd=3)      
      }
    }
    gborders<-rbind(gborders,borders)
    gdistrib<-rbind(gdistrib,distrib)  
  }
  invisible(list(borders=gborders,distrib=gdistrib))
}
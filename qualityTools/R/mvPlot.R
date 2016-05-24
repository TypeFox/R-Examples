mvPlot <-
function(response,fac1,fac2,fac3=NA,fac4=NA,sort=TRUE,col,pch,cex.txt=1,las=1,labels=FALSE,quantile=TRUE,FUN=NA)
{
 old.par <- par(no.readonly = TRUE)                                             #save old par settings
   n=4                                                                          #n number of factors
 if(is.na(fac3) || is.na(fac4))
   n=3
 if(is.na(fac3) && is.na(fac4))
   n=2      
  
 if(is.na(fac3[1]))
  fac3=0
 if(is.na(fac4[1]))
  fac4=0  
 entire_data=data.frame(fac1,fac2,fac3,fac4,response)                           #entire data.frame
 names(entire_data)=c(deparse(substitute(fac1)),deparse(substitute(fac2)),      #set original names
                      deparse(substitute(fac3)),deparse(substitute(fac4)),
                      deparse(substitute(response)))
 original_entire_data=entire_data                                               #save original data
 original_names=names(entire_data)                                              #save original names
 if(sort==TRUE)
 {
  fac_order=c(1,2,3,4)
  names(fac_order)=c(length(unique(fac1)),length(unique(fac2)),
                     length(unique(fac3)),length(unique(fac4)))
  fac_order=fac_order[order(names(fac_order),decreasing=TRUE)]                  #create order-vector to sort factors
  entire_data[,1:4]=entire_data[,fac_order]                                     #sort data by order-vector
  names(entire_data)[1:4]=names(entire_data)[fac_order]                         #sort row names too
  original_names=names(entire_data)                                             #save original names
 } 
 xvalues=as.factor(entire_data[,1])                                             #get values for x-axis
 xlabels=levels(xvalues)                                 
 labels_legend=as.vector(unique(entire_data[,2]))                               #labels for the future legend
 number_of_plots=length(unique(entire_data[,3]))*length(unique(entire_data[,4]))#number of plots needed
 xaxis_lab=c(xlabels,rep(c(NA,NA,xlabels),                                      #create labels for x-axis out of entire data
                         times=length(unique(unique(entire_data[,2])))-1))
 tckx=1:length(xaxis_lab)                                                       #points for tickmarks
 tckx[is.na(xaxis_lab)]=NA  
 levels(xvalues)=1:nlevels(xvalues)                                             #set values for x-axis as numeric
  
 entire_data=cbind(entire_data,xvalues)
 subdata=list(); subsubdata=list(); temp=list() ; plotdata=list(); m=0          #create list of subdata for single plots
 names(entire_data)=c("x-lab","legend","a","b","y","x")
 for(i in 1:length(unique(entire_data[,4])))
 {
  subdata[[i]]=subset(entire_data,entire_data$b==unique(entire_data$b)[i])
  for(k in 1:length(unique(entire_data[,3])))
  {
   subsubdata[[k]]=subset(subdata[[i]],subdata[[i]]$a==unique(subdata[[i]]$a)[k])
   for(j in 1:length(unique(entire_data[,2])))
   { 
    m=m+1
    temp[[m]]=subsubdata[[k]]                                                   #temp is list to make title in plot
    plotdata[[m]]=subset(subsubdata[[k]],                                       #create lists of data for every single plot
                         subsubdata[[k]]$legend==unique(entire_data$legend)[j]) 
    plotdata[[m]]=plotdata[[m]][order(plotdata[[m]][,1]),]                      #order plotdata
   }
  }                                                           
 }
 if(n==3)                                                                       #set device settings
 {
  if(length(as.character(unique(entire_data[,3])))==1)
   par(mfrow=c(1,1)) 
  if(length(as.character(unique(entire_data[,3])))==2)
   par(mfrow=c(1,2))
  if(length(as.character(unique(entire_data[,3])))==3)
   par(mfrow=c(1,3))
  if(length(as.character(unique(entire_data[,3])))==4)
   par(mfrow=c(2,2))
  if(length(as.character(unique(entire_data[,3])))==5)
   par(mfrow=c(2,3))
  if(length(as.character(unique(entire_data[,3])))>=6)
   par(mfrow=c(length(as.character(unique(entire_data[,3]))),1)) 
 }
 if(n==4)
  par(mfrow=c(length(unique(entire_data[,3])),length(unique(entire_data[,4]))))  
 if(missing(col))                                                               #set color (if missing)
  col=c(seq(3,length(unique(fac2))+2))
 if(missing(pch))                                                               #set pch (if missing)
  pch=c(seq(1,length(unique(fac2))),8,16)
  
 m=0
 for(i in 1:(length(plotdata)/length(unique(entire_data[,2]))))
 {
  plot(x=1:length(xaxis_lab),                                                   #empty plot
       y=seq(min(entire_data[,5]),max(entire_data[,5]),length=length(xaxis_lab)),
       axes=FALSE,ylab="",xlab="",col="transparent",
       xlim=c(min(tckx,na.rm=TRUE),max(tckx,na.rm=TRUE)),
       ylim=c(min(entire_data[,5],na.rm=TRUE),max(entire_data[,5],na.rm=TRUE)
              +0.1*(max(entire_data[,5],na.rm=TRUE)-min(entire_data[,5],na.rm=TRUE))))
  axis(1,at=tckx,labels=xaxis_lab,las=las)                                      #add axes
  axis(2)
  if(n==4)
   title(main=paste(names(original_entire_data[3]),"=",as.character(unique(temp[[i*length(unique(entire_data[,2]))]][,3]))," & ",
                    names(original_entire_data[4]),"=",as.character(unique(temp[[i*length(unique(entire_data[,2]))]][,4]))),
         cex.main=1,font.main=1)
  if(n==3)
   title(main=paste(names(original_entire_data[3]),"=",as.character(unique(temp[[i*length(unique(entire_data[,2]))]][,3]))),
         cex.main=1,font.main=1)                                    
  for(j in 1:length(unique(entire_data[,2])))
  {
   m=m+1
   points(x=as.numeric(plotdata[[m]][,6])+(j-1)*(2+nlevels(xvalues)),           #plot single points
          y=plotdata[[m]][,5],col=col[j],pch=pch[j]) 
   if(labels==TRUE && length(row.names(plotdata[[m]]))!=0)                                             #labels
    {
     text(x=as.numeric(plotdata[[m]][,6])+(j-1)*(2+nlevels(xvalues)),
          y=plotdata[[m]][,5], labels=row.names(plotdata[[m]]),cex=0.75,pos=4,xpd=TRUE)
    }
   if(quantile==TRUE)                                                                                  #quantiles
    {
     points(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
            y=rep(quantile(plotdata[[m]][,5],probs=0.5),2),col="gray",pch=3)
     lines(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
           y=rep(quantile(plotdata[[m]][,5],probs=0.5),2),col="gray",lty=2)
     points(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
            y=rep(quantile(plotdata[[m]][,5],probs=0.00135),2),col="gray",pch=3)
     lines(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
           y=rep(quantile(plotdata[[m]][,5],probs=0.00135),2),col="gray",lty=2)
     points(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
            y=rep(quantile(plotdata[[m]][,5],probs=0.99865),2),col="gray",pch=3)
     lines(x=c(0.5,length(unique(entire_data[,1]))+0.5)+(j-1)*(length(unique(entire_data[,1]))+2),
           y=rep(quantile(plotdata[[m]][,5],probs=0.99865),2),col="gray",lty=2)
    }
   if(is.function(FUN))
    {                                                                                                   #insert group information
     points(x=((length(unique(entire_data[,1])))/2+(j-1)*(length(unique(entire_data[,1]))+2))+0.5,
            y=FUN(plotdata[[m]][,5]),col="darkred",pch=15)
     if((m+1)<=length(plotdata) && identical(plotdata[[m]][1:min(length(plotdata[[m]][,3]),length(plotdata[[m+1]][,3])),3],plotdata[[m+1]][1:min(length(plotdata[[m]][,3]),length(plotdata[[m+1]][,3])),3]) 
                                && identical(plotdata[[m]][1:min(length(plotdata[[m]][,4]),length(plotdata[[m+1]][,4])),4],plotdata[[m+1]][1:min(length(plotdata[[m]][,4]),length(plotdata[[m+1]][,4])),4]))
     {
     lines(x=c(((length(unique(entire_data[,1])))/2+(j)*(length(unique(entire_data[,1]))+2))+0.5,
              ((length(unique(entire_data[,1])))/2+(j-1)*(length(unique(entire_data[,1]))+2))+0.5),
           y=c(FUN(plotdata[[m+1]][,5]),FUN(plotdata[[m]][,5])),col="darkred",lty=1)
     }
     fun_val=numeric()
     for(k in 1:length(unique(plotdata[[m]][,1])))
     {
      fun_val[k]=FUN(subset(plotdata[[m]],plotdata[[m]][,1]==unique(plotdata[[m]][,1])[k])[,5])
      points(x=unique(as.numeric(plotdata[[m]][,6])+(j-1)*(2+nlevels(xvalues)))[k],
             y=fun_val[k],col="red",pch=16)
      lines(x=c(unique(as.numeric(plotdata[[m]][,6])+(j-1)*(2+nlevels(xvalues)))[k],unique(as.numeric(plotdata[[m]][,6])+(j-1)*(2+nlevels(xvalues)))[k-1]),
            y=c(fun_val[k],fun_val[k-1]),col="red",lty=2)
     }
    }
   text(x=((length(unique(entire_data[,1])))/2+(j-1)*(length(unique(entire_data[,1]))+2))+0.5,          #legend
        y=max(entire_data[,5],na.rm=TRUE)+0.1*(max(entire_data[,5],na.rm=TRUE)-
        min(entire_data[,5],na.rm=TRUE)),
        labels=as.character(unique(entire_data[,2])[j]),cex=cex.txt)
   } 
 }
 par(old.par)                                                                   #reset par
 invisible(original_entire_data)
}

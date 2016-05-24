Ridler <-
function(Original,pixel,p=TRUE)
{  
	pix<-function(pixel)
{
	missing(pixel)
}
          if(p==TRUE)
        {    
        	display(Original,all= FALSE)
        }
        Original<-(0.1+Original)^2    #accentuate contrast (important especially when pictures are too bright
        Original<-Image(Original,colormode=Grayscale)   #switch to gray scale
        Gray<-Original
        if(p==TRUE)
        {
        	display(Gray,all= FALSE)
        }
        o<-as.numeric(Gray)
        m<-hist(o,50,plot=FALSE)   #Histogram of color
        mu=cumsum(m$counts);
        T1=(sum(m$counts*m$density))/mu[length(mu)] #First threshold
        T1=round(T1)
        mu2=cumsum(m$counts[1:T1])
        MBT=sum(m$density[1:T1]*m$counts[1:T1])/mu2[length(mu2)]   #Color below threshold (background)
        mu3=cumsum(m$counts[T1:length(mu)])
        MFT=sum(m$density[T1:length(mu)]*m$counts[T1:length(mu)])/mu3[length(mu3)] #Color above threshold (foreground)
        T2<-round((MFT+MBT)/2)  #New threshold between the background and foreground
        T<-rep(0,length(mu))
        T[1]<-T1
        T[2]<-T2
        a<-c()
        for(i in 3:length(mu))    #Repeat the calculation of new threshold from background and foreground until the threshold barely changes
        {
            mu2=cumsum(m$counts[1:T[i-1]])
            MBT=sum(m$density[1:T[i-1]]*m$counts[1:T[i-1]])/mu2[length(mu2)]
            mu3=cumsum(m$counts[T[i-1]:length(mu)])
            MFT=sum(m$density[T[i-1]:length(mu)]*m$counts[T[i-1]:length(mu)])/mu3[length(mu3)]
            T[i]=round((MFT+MBT)/2)
            a<-c(a,abs(T[i]-T[i-1]))
        }
        thres<-T[3:length(mu)]
        Thres<-thres[a<=1]
        Final<-Gray
        Final@.Data[Final@.Data<= (Thres[1]-1)/(max(T)-1)]<-0   #Change image into a black and white according to the threshold
        Final@.Data[Final@.Data> (Thres[1]-1)/(max(T)-1)]<-1
        if(p==TRUE)
        {
        display(Final,all=FALSE)
        }
        if(pix(pixel)==TRUE)
        {
            pixel<-length(o)/2 #assumes the half of the entire pictures is the studied area
        }
        Result<-length(as.numeric(Final)[as.numeric(Final)==1])/pixel  #standardized threshold   
        return(Result)
  }

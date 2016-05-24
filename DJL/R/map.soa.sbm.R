map.soa.sbm <-
function(xdata,ydata,date,rts,orientation="n",sg="ssm",mk="dmu"){
  
  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(orientation,c("n","i","o")))){stop('orientation must be "n", "i" or "o".')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  if(is.na(match(mk,c("dmu","eff")))){stop('mk must be either "dmu" or "eff".')}
  
  # Subset index
  till<-function(x,y){
    t<-0
    while(x[t+1]<=y&&t<nrow(x)){t<-t+1}
    return(t)
  }
  
  # Parameters
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);date<-as.matrix(date) # format input data as matrix
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  o<-matrix(c(1:n),ncol=1) # original data order
  
  # Sort data ascending order
  x<-matrix(c(xdata[order(date),]),ncol=m)
  y<-matrix(c(ydata[order(date),]),ncol=s)
  d<-matrix(c(date[order(date),]),ncol=1)
  o<-matrix(c(o[order(date),]),ncol=1)
  
  # max map size
  c<-nrow(unique(d)) 
  ud<-unique(d)
  
  # map frame
  fanta<-matrix(c(NA),nrow=n,ncol=c);colnames(fanta)<-ud
  
  # generate the map
  for(i in 1:c){
    # subset data
    e<-till(d,ud[i])
    x_s<-matrix(x[1:e,],nrow=e)
    y_s<-matrix(y[1:e,],nrow=e)
    
    # run distance measure
    dj<-dm.sbm(x_s,y_s,rts,orientation,se=0,sg)
    
    # soa set
    soa<-which(round(dj$eff,8)==1)
    #soa<-intersect(which(round(dj$eff,8)==1), which(cbind(dj$xslack,dj$yslack)==0))
    
    # fill the map
    if(mk=="dmu"){
      j<-sum(soa>0)
      q<-1
      for(k in 1:j){
        if(ud[i]==ud[1]){fanta[k,1]<-o[soa[k],]}
        else{
          l<-which(fanta[,i-1]==o[soa[k],])
          if(length(l)>0){fanta[l,i]<-o[soa[k],]}
          else{
            p<-n
            while(is.na(fanta[p,i-1])){p<-p-1}
            fanta[p+q,i]<-o[soa[k],]
            q<-q+1
          }
        }
      }
    }
    if(mk=="eff"){
      if(i==1){gsoa<-NULL}
      gsoa<-union(gsoa,soa);l<-length(gsoa)
      fanta[1:l,i]<-dj$eff[gsoa,]
    }
  }
  p<-n;while(is.na(fanta[p,i])){p<-p-1}
  fanta<-fanta[1:p,]
  if(mk=="dmu"){rownames(fanta)<-na.omit(unique(c(fanta)))}
  if(mk=="eff"){rownames(fanta)<-c(o[gsoa,])}
  print(fanta)
}

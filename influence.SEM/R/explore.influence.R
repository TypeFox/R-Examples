explore.influence=function(x,cut.offs="default",plot=TRUE,cook=FALSE)
{  
  # definizione cut-offs
  if ( (length(cut.offs) )==1 && (cut.offs=="default")){
    q25=quantile(x,.25,na.rm=TRUE) 
    q75=quantile(x,.75,na.rm=TRUE)
    if(q75<q25){
      cut.low=q75-(q25-q75)*1.5
      cut.upp=q25+(q25-q75)*1.5
    }
    else{cut.low=q25-(q75-q25)*1.5
         cut.upp=q75+(q75-q25)*1.5}
  }
  else if ( (is.numeric(cut.offs)==TRUE) &&  (length(cut.offs)==2) && (sum(is.na(cut.offs))==0) && (cut.offs[1]<cut.offs[2]) ){
    cut.low=cut.offs[1]
    cut.upp=cut.offs[2] 
  }
  else stop ("\"cut.offs\" must be a vector of 2 numeric elements, with the first element less than the second element")
  if (cook==TRUE) cut.low=max(0,cut.low)
  # plot
  if (plot==TRUE) {
    plot(x,xlab="observations",ylab="influence",ylim=c(min(cut.low,min(x,na.rm=TRUE)),max(cut.upp,max(x,na.rm=TRUE))))
    if (cook==FALSE) abline(h=cut.low,lty=2)
    if ((cook==TRUE) && (cut.low>0)) abline(h=cut.low,lty=2)
    abline(h=cut.upp,lty=2)  
  }
  # output
  ris=NULL
  n=length(x)
  id.row=c(1:n)
  not.allowed=id.row[is.na(x)==TRUE]
  less.cut.low=id.row[(is.na(x)==FALSE) & (x<=cut.low)]
  greater.cut.upp=id.row[(is.na(x)==FALSE) & (x>=cut.upp)]
  ris=list(n=n,cook=cook,cut.low=as.numeric(cut.low),cut.upp=as.numeric(cut.upp),not.allowed=not.allowed,less.cut.low=less.cut.low,greater.cut.upp=greater.cut.upp)  
  return(ris)
}
BaseflowSeparation <-
function(streamflow, filter_parameter=0.925, passes=3){
suppressWarnings(Ends<-c(1,length(streamflow))*rep(1,(passes+1))) # Start and end values for the filter function
suppressWarnings(AddToStart<-c(1,-1)*rep(1,passes))
btP<-streamflow##Previous pass's baseflow approximation
qft<-vector(length=length(streamflow))
bt<-vector(length=length(streamflow))
bt[1]<-if(streamflow[1]<quantile(streamflow,0.25)) streamflow[1] else mean(streamflow)/1.5
##Guess baseflow value in first time step.  
for(j in 1:passes){
for (i in (Ends[j]+AddToStart[j]):Ends[j+1]){
if ((filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]]))>btP[i]){
bt[i]<-btP[i]
} else bt[i]<-filter_parameter*bt[i-AddToStart[j]]+((1-filter_parameter)/2)*(btP[i]+btP[i-AddToStart[j]])
qft[i]<-streamflow[i]-bt[i]
}
if (j<passes){
btP<-bt
bt[Ends[j+1]]<-if(streamflow[Ends[j+1]]<mean(btP))streamflow[Ends[j+1]]/1.2 else mean(btP)
##Refines the approximation of end values after the first pass
}
}
f <- data.frame(bt,qft)
return(f)
}

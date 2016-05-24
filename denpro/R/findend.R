findend<-function(inde,left,right,low,upp,N){
#
lenn<-length(inde)
current<-1
dep<-1
if ((left[current]==0) && (right[current]==0)){
   exists<-FALSE
   }
else{
   exists<-TRUE
}
while ((exists) && ((left[current]>0) || (right[current]>0))){
     mid<-(low[current]+upp[current])/2 
     direc<-depth2com(dep,N)$direc 
     if (inde[direc]<=mid){
           if (left[current]>0){ 
                current<-left[current]
                dep<-dep+1
           }
           else{ 
                exists<-FALSE
           }
     }
     else{   
           if (right[current]>0){
                  current<-right[current]
                  dep<-dep+1
           }
           else{
               exists<-FALSE
           }
     }
}
return(list(exists=exists,location=current,dep=dep))
}

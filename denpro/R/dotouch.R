dotouch<-function(inde1,inde2,direction){
#
d<-length(inde1)
touch<-TRUE
for (i in direction:d){
     if ((inde1[i]>inde2[i]+1) || (inde1[i]<inde2[i]-1)){
           touch<-FALSE
     }
}
return(touch)
}

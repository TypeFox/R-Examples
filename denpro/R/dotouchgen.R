dotouchgen<-function(indelow1,indeupp1,indelow2,indeupp2,direction){
#
epsi<-0
d<-length(indelow1)
touch<-TRUE
i<-1
while (i<=d){
     if ((i != direction) &&
         ((indelow1[i]>indeupp2[i]+epsi) || (indeupp1[i]<indelow2[i]-epsi))){
              touch<-FALSE
     }
     i<-i+1
}
return(touch)
}

mtest<-function(profile,n){
#
parents<-profile$parent
volumes<-profile$volume
levels<-profile$level
#
nodelkm<-length(parents)
lowmasses<-matrix(1,nodelkm,1)
for (i in 1:nodelkm){
   par<-parents[i]
   if (par==0){
      lowmasses[i]<-levels[i]*volumes[i]
   }
   else{
      lowmasses[i]<-levels[par]*volumes[i]
   }
}
testcrit<-sqrt(lowmasses/n)
#
return(t(testcrit))
#return(t(lowmasses))
}







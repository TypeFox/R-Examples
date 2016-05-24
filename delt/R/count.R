count<-function(beg,end){
#Laskee hilaan kuuluvien pisteiden lkm:n
#
#hila on 2-vector
#
#returns integer >1
#
if (is.na(beg)) ans<-0              #or end==NA
else if ((beg==0) && (end==0)) ans<-0
else ans<-end-beg+1
return(ans)
}

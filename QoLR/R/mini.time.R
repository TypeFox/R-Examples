mini.time <-
function(vector){
m1=matrix(vector,nrow=length(vector),ncol=length(vector),byrow=TRUE)
m1[upper.tri(m1)] <- NA
vector2=apply(m1,1,min,na.rm=TRUE)
vector2
}

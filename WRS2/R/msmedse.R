msmedse <-
function(x, sewarn = TRUE){
#
# Compute  standard error of the median using method
# recommended by McKean and Shrader (1984).
#
x=elimna(x)
chk=sum(duplicated(x))
if (sewarn) {
if(chk>0){
warning("Tied values detected. Estimate of standard error might be inaccurate.")
}}
y<-sort(x)
n<-length(x)
av<-round((n+1)/2-qnorm(.995)*sqrt(n/4))
if(av==0)av<-1
top<-n-av+1
sqse<-((y[top]-y[av])/(2*qnorm(.995)))^2
sqse<-sqrt(sqse)
sqse
}

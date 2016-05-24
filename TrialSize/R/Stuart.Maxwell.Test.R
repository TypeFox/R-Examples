Stuart.Maxwell.Test <-
function(noncen,p.ij,p.ji,r){
for (i in 1:1500){
b<-pchisq(qchisq(0.05,df=3*(3-1)/2,ncp=0,lower.tail=FALSE),df=3*(3-1)/2,ncp=i/100)
print(c(i,b))
}
# noncen=10.89

p.ij<-rbind(c(3,4,4),c(2,3,3),c(1,2,3))/25
p.ji<-rbind(c(3,4,4),c(2,3,3),c(1,2,3))/25

sum=0
n<-noncen*(
for (j in 1:r){
for (i in 1:j){
sum=sum+(p.ij[i,j]-p.ji[j,i])^2/(p.ij[i,j]+p.ji[j,i])
}}
)^-1
}

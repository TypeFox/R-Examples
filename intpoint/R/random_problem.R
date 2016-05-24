random_problem <-
function(n){
for(i in 1:n){
if(i==1) A<-runif(n) else
A<-rbind(A,runif(n))
}
c<-runif(n)
v<-runif(n)
b<-array(0,c(n,1))
for(i in 1:n){
for(j in 1:n){
b[i]<-b[i]+A[i,j]*v[j]
}
}
sol<-interior_point(1,c,bm=b,m=A)
if(n<=5)
return(list("THE RANDOM PROBLEM GENERATED IS:","Constraints matrix A"=A,"
Vector of coefficients of the objective function c"=c,"Vector of right hand side constants b"=b,"SOLUTION:",sol))
else return(sol)
}

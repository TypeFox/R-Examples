Nucleolus <-function(x,type="Gains"){

####
n<-x$Ag
as.matrix(x$Binary)->z
as.vector(z)->V
#####

vec <- c(0, 1)
lst <- lapply(numeric(n), function(x) vec)
Amat<-as.matrix(expand.grid(lst))[-1,]
N<-dim(Amat)[1]

if (type=="Gains"){
Amat<-cbind(Amat,c(rep(-1,N-1),0))
lprec<-make.lp(0, n+1)
set.objfn(lprec, c(rep(0,n),-1))
add.constraint(lprec, c(rep(0,n),1), ">=", 0) 
for (i in 2: N-1){
 add.constraint(lprec, Amat[i,], ">=", V[i])
}

add.constraint(lprec, Amat[N,], "=",V[N])
Lim<-NULL
for (i in seq(1,n)){
	Lim<-c(Lim,V[2^(i-1)])
}

Lim<-c(Lim,0)

set.bounds(lprec, lower = c(rep(0,n),0), upper = c(rep(Inf,n),Inf))
name.lp(lprec, "Nucleolus of a gains game ")
lp.control(lprec,sense='min',verbose="normal",simplextype=c("dual","primal"))
print(lprec)

FLP<-lprec


nucleolus<-NULL

for (i in seq(1,n)) {
	solve(lprec); 
    select.solution(lprec)
	X<-   get.variables(lprec)[i]
	E<- - get.variables(lprec)[n+1]
	S<-c(V[2^(i-1)],X,E)
	nucleolus<-rbind(nucleolus,S)
	set.mat(lprec,2^(i-1)+1,n+1,0)
	set.mat(lprec,2^n,i,0)
	set.mat(lprec,2^n-2^(i-1),n+1,0)
	set.constr.type(lprec,"=", 2^(i-1)+1)
	RHS<-as.matrix(get.rhs(lprec))
	RHS[2^(i-1)+1,]<- X
	RHS[2^n-2^(i-1)] <- RHS[2^n-2^(i-1)] + E
	RHS[2^n] <- RHS[2^n] - X
	set.rhs(lprec,RHS)
}
}

if (type=="Cost"){
	Amat<-cbind(Amat,rep(1,N),rep(-1,N))
	lprec<-make.lp(0, n+2)
set.objfn(lprec, c(rep(0,n),1,-1))

print(lprec)

add.constraint(lprec, c(rep(0,n),0,1), ">=", 0) 



for (i in 2: N-1){
 add.constraint(lprec, Amat[i,], "<=", V[i])
}

add.constraint(lprec, c(rep(1,n),0,0), "=",V[N])


Lim<-NULL
for (i in seq(1,n)){
	Lim<-c(Lim,V[2^(i-1)])
	
}
Lim<-c(Lim,0)

set.bounds(lprec, lower = c(rep(0,n),0,0), upper = c(rep(Inf,n),Inf,Inf))
name.lp(lprec, "Nucleolus of a cost game ")
lp.control(lprec,sense='max',verbose="normal",simplextype=c("primal","dual"))
print(lprec)
FLP<-lprec


nucleolus<-NULL

for (i in seq(1,n)) {
	
	solve(lprec); 
    print(get.variables(lprec))
	X<-  get.variables(lprec)[i]
	E<-  get.objective(lprec)
	S<-c(V[2^(i-1)],X,E)
	nucleolus<-rbind(nucleolus,S)

	
	set.mat(lprec,2^(i-1)+1,n+1,0)
	set.mat(lprec,2^(i-1)+1,n+2,0)
	
	set.mat(lprec,2^n,i,0)
	set.mat(lprec,2^n-2^(i-1),n+1,0)
	set.mat(lprec,2^n-2^(i-1),n+2,0)
	
	set.constr.type(lprec,"=", 2^(i-1)+1)
	
	RHS<-as.matrix(get.rhs(lprec))
	RHS[2^(i-1)+1,]<- X
	RHS[2^n-2^(i-1)] <- RHS[2^n-2^(i-1)] - E
	RHS[2^n] <- RHS[2^n] - X
	
	set.rhs(lprec,RHS)
	
		
}
	
	}

delete.constraint(FLP, 1)
nucleolus<-as.matrix(nucleolus)
colnames(nucleolus)<-c("v(S)","x(S)","Ei")
output<-list(nucleolus=nucleolus,mode=type)
class(output)<-"Nucleolus"
return(output)
}
BAGIDIS_2D = function(sig.out.1,sig.out.2,wk,lambdaV,lambdaH){
# lambdaV and lambda H must be in 0;0.5 if they are equal 

sig.out.1=as.matrix(sig.out.1)  
sig.out.2=as.matrix(sig.out.2)  
N= nrow(sig.out.1)

#wk = c(1,1,rep(0,N-2))

dist=rep(0,N)
for (i in 1:N){
dist[i]= wk[i]*sqrt( 
               lambdaV*sum( (sig.out.1[i,c(1,3)]-sig.out.2[i,c(1,3)])^2  )
		  +lambdaH*sum( (sig.out.1[i,c(2,4)]-sig.out.2[i,c(2,4)])^2  )  
              +(1-lambdaV-lambdaH)*(sig.out.1[i,5]-sig.out.2[i,5])^2     )   }
dist=sum(dist)
return(dist)}

#============================================================================================


semimetric.BAGIDIS_2D_BD=function(list.sig.1,list.sig.2=list.sig.1,wk,lambdaV,lambdaH){
  
DistMatrix =matrix(NA,nrow=length(list.sig.1),ncol=length(list.sig.2))

for (i in 1:length(list.sig.1)){
  for (j in 1:length(list.sig.2)){
  DistMatrix[i,j]=BAGIDIS_2D(list.sig.1[[i]],list.sig.2[[j]],wk,lambdaV,lambdaH)
  }
}
return(DistMatrix)
}

#============================================================================================

semimetric.BAGIDIS_2D_Array=function(Array1,Array2=Array1,wk,lambdaV,lambdaH){

out.list.1 =apply(Array1,3,BUUHWE_2D)
list.sig.1=lapply(out.list.1,Signature_2D)

if (!identical(Array1,Array2)){
out.list.2 =apply(Array2,3,BUUHWE_2D)
list.sig.2=lapply(out.list.2,Signature_2D)} else{list.sig.2=list.sig.1}
  
DistMatrix = semimetric.BAGIDIS_2D_BD(list.sig.1,list.sig.2,wk,lambdaV,lambdaH)
return(DistMatrix)
}

#===================================================================

#==================================================================
# interface on semimetric.BAGIDIS_2D with images encoded as vectors, 
# stored row by row.
#==================================================================

semimetric.BAGIDIS_2D_Vec=function(Matrix1,Matrix2=Matrix1,NbROW, wk,lambdaV,lambdaH){
  
Array1 = NULL
for (i in 1:nrow(Matrix1)){
  Matrix= matrix(Matrix1[i,], nrow=NbROW)
  Array1 = abind(Array1, Matrix, along =3)
}
if (!identical(Matrix1,Matrix2)){
Array2 = NULL
for (i in 1:nrow(Matrix2)){
  Matrix= matrix(Matrix2[i,], nrow=NbROW)
  Array2 = abind(Array2, Matrix, along =3)
}} 
else
{Array2 =Array1}
DistMatrix = semimetric.BAGIDIS_2D_Array(Array1,Array2,wk,lambdaV,lambdaH)
return(DistMatrix)
}

#=================================================================================
# Final Interface for semimetric.BAGIDIS
#================================================================================

semimetric.BAGIDIS_2D = function(Data1,Data2=Data1,
                                 NbROW = NULL, wk=NULL,lambdaV=NULL,lambdaH=NULL, Type =c('Array', 'Vector', 'BD','Flex')){
                                 
                                 
if(is.null(wk)){
if(Type=='Array'){
N = nrow(Data1[,,1])* ncol(Data1[,,1]) -1
wk = log(N+1-(1:N))/log(N+1)}
if(Type=='Vector'){
N = ncol(Data1)-1
wk = log(N+1-(1:N))/log(N+1)}
if(Type=='BD'){
N = nrow(Data1[[1]])
wk = log(N+1-(1:N))/log(N+1)}
else if(Type=='Flex'){
wk= matrix(1, nrow= nrow(Data1[[1]]), ncol =ncol(Data1[[1]]))}
}

if (Type!='Flex'){

if (!is.null(lambdaH)){
if (lambdaH <0){stop("lambdaH cannot be negative")}
if (lambdaH >1){stop("lambdaH cannot be >1")} }
if (!is.null(lambdaV)){
if (lambdaV <0){stop("lambdaV cannot be negative")}
if (lambdaV >1){stop("lambdaV cannot be >1")} }

if (is.null(lambdaV)){
  if (is.null(lambdaH)){
    lambdaV=1/3
    lambdaH=1/3}
  else{
  lambdaV= max(1-lambdaH)/2
}}

if (is.null(lambdaH)){
  if (!is.null(lambdaV)){
  lambdaH= (1-lambdaV)/2
}}

if (lambdaV+lambdaH >1){stop("lambdaV+lambdaH cannot be >1")} 

}

                                 
                                 
if (Type=='Vector') {  
DistMatrix = semimetric.BAGIDIS_2D_Vec(Data1,Data2,NbROW, wk,lambdaV,lambdaH)
return(DistMatrix)  
} 
else if (Type=='Array'){
DistMatrix = semimetric.BAGIDIS_2D_Array(Data1,Data2,wk,lambdaV,lambdaH)
return(DistMatrix)    
}
else if (Type=='BD'){
DistMatrix = semimetric.BAGIDIS_2D_BD(Data1,Data2,wk,lambdaV,lambdaH)
return(DistMatrix)    
}
else if (Type=='Flex'){
DistMatrix = semimetric.BAGIDIS_2D_Flex(Data1,Data2,wk)
return(DistMatrix)    
}

}

#=================================================================================
# Flexible semimetric
#================================================================================


BAGIDIS_2D_Flex= function(sig.out.1,sig.out.2,wk){

sig.out.1=as.matrix(sig.out.1)  
sig.out.2=as.matrix(sig.out.2)  

dist= sum(wk*abs(sig.out.1 -sig.out.2))

return(dist)}




semimetric.BAGIDIS_2D_Flex=function(list.sig.1,list.sig.2=list.sig.1,wk){


# with Data1 = list(x1,y1,x2,y2,Details)
# wk = matrix / row index = rank; col index = x1,y1,x2,y2,Details
  
DistMatrix =matrix(NA,nrow=length(list.sig.1),ncol=length(list.sig.2))

for (i in 1:length(list.sig.1)){
  for (j in 1:length(list.sig.2)){
  DistMatrix[i,j]=BAGIDIS_2D_Flex(list.sig.1[[i]],list.sig.2[[j]],wk)
  }
}
return(DistMatrix)
}

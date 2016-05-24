buildSb <-
function(bA = NULL, A = NULL, bm = NULL,
	m = NULL, bM = NULL, M = NULL)
{
#We check if A is not null
if(!is.null(A))
{
	if(is.vector(A)){
	aux<-array(0,c(1,length(A)))
	for(i in 1:length(A))
	aux[i]<-A[i]
	aux2<-array(0,c(1,length(A)))
	for(i in 1:length(A))
	aux2[i]<-(-A[i])
	S<-array(0,2,2)
	S<-rbind(aux,aux2)
	b1<-c(bA,-bA)
	J1<-joint(S,b1,m,bm)
	S<-J1[[1]]
	b<-J1[[2]]
	J2<-joint2(S,b,M,bM)
	S<-J2[[1]]
	b<-J2[[2]]
	}
	else  if(is.array(A)){
	nr<-nrow(A)
	S<-array(0,(2*nr),2)
	for (i in 1:nr) S[i,]<-A[i,]
	for (i in (nr+1):(2*nr)) S[i,]<-(-A[i,])
	b1<-c(bA,-bA)
	J1<-joint(S,b1,m,bm)
	S<-J1[[1]]
	b<-J1[[2]]
	J2<-joint2(S,b1,M,bM)
	S<-J2[[1]]
	b<-J2[[2]]
	}
}
# If A is null then we check if m is not null
if ( is.null(A) & !is.null(m) )
	{
	S<-m
	b<-bm
	J2<-joint2(S,b,M,bM)
	S<-J2[[1]]
	b<-J2[[2]]
	}
if( is.null(A) & is.null(m) & !is.null(M) )
# If A and m are both null, then M  can not be null
{
S<-(-M)
b<-(-bM)
}
return( list(S1=S, b1=b))
}

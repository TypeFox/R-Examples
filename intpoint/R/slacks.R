slacks <-
function(m,M){
s<-0
n<-0
if(!is.null(m)){
if(is.vector(m)){
n<-1
}
else if(is.array(m)){
for(i in 1:length(m))
s<-s+abs(m[i])
if (s!=0) n<-nrow(m) else n<-0
}
}
s<-0
if(!is.null(M)){
if(is.vector(M)){
n<-n+1
}
else if(is.array(M)){
for(i in 1:length(M))
s<-s+abs(M[i])
if (s!=0) n<-n+nrow(M)
}
}
return(n)
}

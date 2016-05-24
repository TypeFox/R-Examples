fillcomb <-
function(r){

out<-matrix(0,ncol(r),max(r))

for(i in 1:ncol(r)){
	out[i,r[,i]]<-1
}

out
}


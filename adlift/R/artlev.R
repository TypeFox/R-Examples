`artlev` <-
function(y,rem){

l<-length(rem)

remaining<-1:l
end<-0
n<-1
q<-.5
p<-list()

while(!end){
	r<-which(y<=quantile(y,q))
	s<-intersect(remaining,r)
	if(length(s)<10){
		p[[n]]<-rem[remaining]
		end<-1
	}
	else{
		p[[n]]<-rem[sort(s)]
		remaining<-setdiff(remaining,s)
		q<-(1+q)/2
		n<-n+1
	}
	if(length(remaining)==0){
		end<-1
	}
}

p
}


corsim<-function(x,lambda,mu,missing,told=0,tyoung=0){
x<-sort(x,decreasing=TRUE)

#make vector ranks and times
n<-length(x)+1
if (told == 0){
	rho<- 1-missing/(missing+n)
	origin<-0
	#sample origin conditioned on n and x_1
	while (x[1]>origin){
		r <- runif(1,0,1)
		if (lambda>mu) {
			origin <- log((-lambda * rho - lambda * r^(1/n) + mu * r^(1/n) + lambda * rho * r^(1/n))/(lambda * rho * (-1 + r^(1/n)))) / (lambda - mu)
		} else {
			origin<- -(r^(1/n)/(lambda *(-1 + r^(1/n)* rho)))
		}
	}
	if (tyoung==0){
		ranks<-0:(length(x)+1)
		times<-c(origin,x,0)
	} else {
		missyoung<-length(which(x<tyoung))
		times<-c(0,x[1:(length(x)-missyoung)],tyoung)
		ranks<-1:length(times)
		ranks<-ranks-1
	}
} else {
	missold<-length(which(x>told))
	missyoung<-length(which(x<tyoung))
	if (missold<length(x)){
		times<-c(told,x[(missold+1):(length(x)-missyoung)],tyoung)
	} else {
		times<-c(told,tyoung)	
	}
	ranks<-1:length(times)
	ranks<-ranks+missold-1
	}
#after times[i] we have ranks[i]+1 lineages

while(missing>0){
#distrranks[i]: prob insert between ranks[i] and ranks[i+1]
if (length(ranks)>2){
distrranks<-vector()
for (i in 2:length(ranks)){
	temp <- ranks[i] * (intp1(times[i-1],lambda,mu) - intp1(times[i],lambda,mu))
	distrranks<-c(distrranks,temp)
	}

distrranks<-distrranks/sum(distrranks)
for (i in 2:length(distrranks)){distrranks[i]<-distrranks[i]+distrranks[i-1]}
r <- runif(1,0,1)
addrank<-min(which(distrranks>r))
# addrank=k means adding between ranks[k] and ranks[k+1] in time
# means adding to ranks[k+1] lineages
} else {addrank<-1}
r <- runif(1,0,1)
const<-intp1(times[addrank],lambda,mu) - intp1(times[(addrank+1)],lambda,mu)
temp<-  intp1(times[(addrank+1)],lambda,mu)/const
xnew<- 1/(mu-lambda)*log((1-(r+temp)*const*lambda)/(1-(r+temp)*const*mu))
x<-c(x,xnew)
x<-sort(x,decreasing=TRUE)
missing<-missing-1
}
x
}
accept.reject <-
function(prob,a,b,n){
seq.unif<-seq(a:b)
c<-max(prob)/(1/length(seq.unif))
samp<-NULL
while (length(samp)<n){
y<-sample(a:b,1,replace = TRUE)
if (a==0){ 
y2<-y+1
}else{
y2<-y
}
u<-runif(1)
if (u<=(prob[y2]/(c*(1/length(seq.unif))))) samp<-c(samp,y)
}
samp
}


`audps` <-
function(evaluation, dates, type = "absolute") {
if(!(is.matrix(evaluation) | is.data.frame(evaluation))){
evaluation<-rbind(evaluation)
}
n<-length(dates)
k<-ncol(evaluation)
if (n!=k) {
cat("Error:\nThe number of dates of evaluation \nmust agree with the number of evaluations\n")
return()
}
d1<-(dates[2]-dates[1])/2
d2<-(dates[n]-dates[n-1])/2
d<-d1+d2+dates[n]-dates[1]
audps<-0
for (i in 1:(n-1)) {
audps<- audps + evaluation[,i]*(dates[i+1]-dates[i])
}
audps<- audps + evaluation[,n]*(dates[n]-dates[n-1])
if (type =="relative" ) audps <-audps/(d*100)
if (type =="absolute" | type =="relative" ) {
return(audps)
}
else cat("Error: type is 'absolute' or 'relative'\n\n")
}


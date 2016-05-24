`audpc` <-
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
audpc<-0
area.total<- 100*(dates[n]-dates[1])
for (i in 1:(n-1)) {
audpc<- audpc + (evaluation[,i]+evaluation[,i+1])*(dates[i+1]-dates[i])/2
}
if (type =="relative" ) audpc <-audpc/area.total
if (type =="absolute" | type =="relative" ) {
return(audpc)
}
else cat("Error: type is 'absolute' or 'relative'\n\n")
}


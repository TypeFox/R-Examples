buildMulti <-
function(dat,copy,multi){
 if (missing(dat)) stop("Need the data frame")
 if (missing(copy)) stop("Need column numbers to copy")
 if (missing(multi)) stop("Need multi list(vector of numbers,vector of names)")
col_mult<- matrix(0,ncol=1,nrow=length(multi[[1]]))  #columns numbers
col_sum<-  matrix(0,ncol=1,nrow=length(multi[[1]]))  #counts in columns
for (i in 1:length(multi[[1]])) {
  col_mult[i,]<- which(colnames(dat)==multi[[2]][i]) } #colnames in list
for (i in 1:length(multi[[1]])) {
  col_sum[i,]<- sum(dat[col_mult[i,]]) }   #get length in counts for each column
org<- do.call("rbind", replicate(length(multi[[1]]),dat[,copy],simplify=F))
tim<-  c()
 for (i in 1:length(multi[[1]])) { tim<- c(tim,dat[,col_mult[i]]) }
status<- rep(multi[[1]],each=length(dat[,1]))
org<- cbind(tim,status,org)  #column is named status
expa<- org[rep(1:length(org[,1]),org[,1][1:length(org[,1])]),]   #tim=,1
expa<- expa[,-c(1)] #remove tim
  invisible(expa) #does not print
}

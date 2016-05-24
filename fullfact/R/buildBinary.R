buildBinary <-
function(dat,copy,one,zero){  #label object, 'one' column, 'zero' column
 if (missing(dat)) stop("Need the data frame")
 if (missing(copy)) stop("Need column numbers to copy")
 if (missing(one)) stop("Need column name of 'one' counts")
 if (missing(zero)) stop("Need column name of 'zero' counts")
  labels_one<- cbind(1,dat[,copy]);colnames(labels_one)[1]<- "status"
  labels_zero<- cbind(0,dat[,copy]);colnames(labels_zero)[1]<- "status"
one2<- which(colnames(dat)==one)
zero2<- which(colnames(dat)==zero)
  dat_one<- labels_one[rep(1:length(dat[,one2]),dat[,one2][1:length(dat[,one2])]),]
  dat_zero<- labels_zero[rep(1:length(dat[,zero2]),dat[,zero2][1:length(dat[,zero2])]),]
  expa<-rbind(dat_one,dat_zero)
  invisible(expa) #does not print
}

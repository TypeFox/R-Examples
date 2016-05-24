selby<-function(m,grpc,coln){
#
#
#  A commmon situation is to have data stored in an n by p matrix where
#  one or more of the columns are  group identification numbers.
#  This function groups  all values in column coln according to the
#  group numbers in column grpc and stores the  results in list mode.
#
#  More than one column of data can sorted
#
# grpc indicates the column of the matrix containing group id number
#
if(is.null(dim(m)))stop("Data must be stored in a matrix or data frame")
if(is.na(grpc[1]))stop("The argument grpc is not specified")
if(is.na(coln[1]))stop("The argument coln is not specified")
if(length(grpc)!=1)stop("The argument grpc must have length 1")
x<-vector("list")
grpn<-sort(unique(m[,grpc]))
it<-0
for (ig in 1:length(grpn)){
for (ic in 1:length(coln)){
it<-it+1
flag<-(m[,grpc]==grpn[ig])
x[[it]]<-m[flag,coln[ic]]
}}
list(x=x,grpn=grpn)
}
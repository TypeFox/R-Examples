selby2<-function(m,grpc,coln=NA){
# Create categories according to the grpc[1] and grpc[2] columns
# of the matrix m. The function puts the values in column coln into
# a vector having list mode.
#
if(is.na(coln))stop("The argument coln is not specified")
if(length(grpc)>4)stop("The argument grpc must have length less than or equal to 4")
x<-vector("list")
ic<-0
if(length(grpc)==2){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2])]==c(cat1[i1],cat2[i2]))==2){
it<-it+1
temp[it]<-m[i,coln]
}
}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2]),1,2)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2]))
}
}}
}
if(length(grpc)==3){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
cat3<-selby(m,grpc[3],coln)$grpn
x<-vector("list")
ic<-0
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
for (i3 in 1:length(cat3)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2],grpc[3])]==c(cat1[i1],cat2[i2],cat3[i3]))==3){
it<-it+1
temp[it]<-m[i,coln]
}}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3]),1,3)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3]))
}}}}
}
if(length(grpc)==4){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
cat3<-selby(m,grpc[3],coln)$grpn
cat4<-selby(m,grpc[4],coln)$grpn
x<-vector("list")
ic<-0
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
for (i3 in 1:length(cat3)){
for (i4 in 1:length(cat4)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2],grpc[3],grpc[4])]==c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))==4){
it<-it+1
temp[it]<-m[i,coln]
}}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]),1,4)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))
}}}}}
}
list(x=x,grpn=grpn)
}

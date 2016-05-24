format_p.ade <-
function(x, pgits=4, digits=2){
x<-unlist(x)
if(!is.vector(x) & !is.matrix(x)) stop('x must be vector or matrix')


#########################################################
# split in single numbers
if(length(x)>1){

if(is.vector(x)){
y<-x
for(i in 1:length(x)){
y[i]<- format_p.ade(x[i], pgits, digits)
}
return(y)
}

if(is.matrix(x)){
M<-x
for(i in 1:dim(x)[1]){
for(j in 1:dim(x)[2]){
M[i, j]<- format_p.ade(x[i, j], pgits, digits)
}}
return(M)
}

}
#########################################################




####################################
# if single number
if(length(x)==1){

if(is.na(x))   return(x)
x<-as.numeric(x)

if(digits%%1==0){
# Meine Rundung
if(pgits%%1==0){
y<-base::format(x, digits=digits, scientific = FALSE)
if(y!=1){
for(i in 1:(digits-1)){
yy<-sub('0.', '', y)
sdigs<- nchar(as.numeric(yy))
if(sdigs<digits)
y<-paste(y, '0', sep='')
}
y[as.numeric(y)<(1/(10^pgits))]<-paste('<', base::format((1/(10^pgits)), scientific = FALSE), sep='')
}
if(y==1){
y<-'1.'
for(i in 1:digits){  y<-paste(y, '0', sep='')
}
}
}


# Halbe Rundung
if(pgits%%1==0.5){
y<-base::format(x, digits=digits, scientific = FALSE)
if(y!=1){
for(i in 1:(digits-1)){
yy<-sub('0.', '', y)
sdigs<- nchar(as.numeric(yy))
if(sdigs<digits) y<-paste(y, '0', sep='')
}
y[as.numeric(y)<(1/(10^trunc(pgits)))/2]<-paste('<', base::format((1/(10^trunc(pgits)))/2, scientific = FALSE), sep='')
}
if(y==1){
y<-'1.'
for(i in 1:digits){  y<-paste(y, '0', sep='')
}
}
}
}


# Schnabel Rundung
if(digits==1.5){
if(pgits%%1==0){
if(as.numeric(x)>0.01)  y<-base::format(x, digits=2, scientific = FALSE)
if(as.numeric(x)<=0.01) y<-base::format(x, digits=1, scientific = FALSE)
if(y!=1){
if(y>0.01){
yy<-sub('0.', '', y)
sdigs<- nchar(as.numeric(yy))
if(sdigs<2) y<-paste(y, '0', sep='')
}
y[as.numeric(y)<(1/(10^pgits))]<-paste('<', base::format((1/(10^pgits)), scientific = FALSE), sep='')
}
if(y==1) y<-'1.0'
}

}

return(y)
}
####################################

}

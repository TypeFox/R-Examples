format_n.ade <-
function(x, digits=2, scientific = FALSE){
x<-unlist(x)
if(!is.vector(x) & !is.matrix(x)) stop('x must be vector or matrix')

#########################################################
# split in single numbers
if(length(x)>1){

if(is.vector(x)){
y<-x
for(i in 1:length(x)){
y[i]<- format_n.ade(x[i], digits, scientific)
}
return(y)
}

if(is.matrix(x)){
M<-x
for(i in 1:dim(x)[1]){
for(j in 1:dim(x)[2]){
M[i, j]<- format_n.ade(x[i, j], digits, scientific)
}}
return(M)
}


}
#########################################################



if(length(x)==1){

#########################################################
# pass flat zero
if(is.na(x))  return(x)
if(x==0)      return(as.character(x))
#########################################################


#########################################################
# if scientific representation
if(scientific){
fi3<-base::format.info(x)[3]
if(!is.na(fi3)){
if(fi3==0)      scientific<-FALSE
}
if(is.na(fi3))  scientific<-FALSE
}

if(scientific){
return(base::format(x, digits=digits, scientific=TRUE))
}
#########################################################



#########################################################
# without scientific representation
if(!scientific){
x<- as.numeric(x)


###########################
# if no scientific representation used
if(!is.na(base::format.info(x)[3])) {


if(base::format.info(x)[3]==0) {

# calc r to round
digs  <-base::format.info(x)[2]
if(is.na(digs))  digs<-0
ints<- base::format.info(trunc(abs(x)))[1]
if(trunc(abs(x))==0) ints<-0
if(ints>=digits)  r<-0
if(ints< digits)  r<- digits-ints
if(trunc(abs(x))==0) {
y<-base::format(abs(x),  scientific = FALSE)
digs<-nchar(base::format(as.numeric(gsub('^0.', '', y)), scientific = FALSE))
undigs<-nchar(gsub('^0.', '', y))-digs
r<-undigs+digits
}


y<- base::format(base::round(x,digits=r),nsmall=r, scientific = FALSE)
return(y)
}
###########################



###########################
# if scientific representation is used
if(base::format.info(x)[3]>0){
if(abs(x)<1){
y<-base::format(x,  scientific = FALSE)
ycut<-gsub('-', '', y)
ycut<-gsub('^0.', '', ycut)
digs<-nchar(base::format(as.numeric(ycut), scientific = FALSE))
undigs<-nchar(ycut)-digs

y<-base::format(base::round(abs(x), digits=(undigs+digits)) , scientific = FALSE)
digs<-nchar(base::format(as.numeric(gsub('^0.', '', y)), scientific = FALSE))
undigs<-nchar(gsub('^0.', '', y))-digs
if(digs<digits){
for(i in 1:(digits-digs))
y<-paste(y, '0', sep='')
}
if(sign(x)==-1)   y<-paste('-', y,  sep='')
return(y)
}

#  if a large number
if(abs(x)>=1) {
y<- base::format(base::round(x), scientific = FALSE)
return(y)

}
}   }
###########################




}
#########################################################
}

}

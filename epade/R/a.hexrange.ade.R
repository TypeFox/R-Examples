a.hexrange.ade <-
function(x,hmin,hmid, hmax, rangemy=c(0,1), abs=FALSE){

################################################################################
hexr <- function(x,hmin,hmid, hmax, rangemy=c(0,1)){
if (!is.na(x)){

if(x>= rangemy[2]) return(hmax)
if(x<= rangemy[1]) return(hmin)
if(x== mean(rangemy, na.rm=TRUE)) return(hmid)
if(x> rangemy[2]) x<-rangemy[2]
if(x< rangemy[1]) x<-rangemy[1]

x<-(x-(rangemy[1]+(rangemy[2]-rangemy[1])/2))/(rangemy[2]-rangemy[1])
x<-x*2

hmin<- gsub('#','',hmin)
hmid<- gsub('#','',hmid)
hmax<- gsub('#','',hmax)

######################
##################
#mod function
mod.ade<-function(x,z){
mod <-x-floor(x/z)*z
return(mod)
}
##################
######################

######################
##################
hex.to.dec <- function(x){
n<- nchar(x)
d<-0
for(i in 1:n){
z<-substr(x, start=i, stop=i)
if(z=='A') z<-10
if(z=='B') z<-11
if(z=='C') z<-12
if(z=='D') z<-13
if(z=='E') z<-14
if(z=='F') z<-15
z<-as.numeric(z)
d=d+z*16^(n-i)
}
return(d)
}
##################
######################

######################
##################
dec.to.hex <- function(x){

hex<-NULL
i=1
while(x!=0){
hex[i]<- x-(trunc(x/16)*16)
x<-trunc(x/16)
i=i+1
}
n<-length(hex)
hex[hex==10] <-'A'
hex[hex==11] <-'B'
hex[hex==12] <-'C'
hex[hex==13] <-'D'
hex[hex==14] <-'E'
hex[hex==15] <-'F'


string<-NULL
for(i in 1:n){
string <- paste(string, hex[(n+1-i)], sep='')
}
return(string)
}
##################
######################

rmin<- substr(hmin, start=1, stop=2)
rmid<- substr(hmid, start=1, stop=2)
rmax<- substr(hmax, start=1, stop=2)
gmin<- substr(hmin, start=3, stop=4)
gmid<- substr(hmid, start=3, stop=4)
gmax<- substr(hmax, start=3, stop=4)
bmin<- substr(hmin, start=5, stop=6)
bmid<- substr(hmid, start=5, stop=6)
bmax<- substr(hmax, start=5, stop=6)

rr1<-hex.to.dec(rmid)-hex.to.dec(rmin)
gr1<-hex.to.dec(gmid)-hex.to.dec(gmin)
br1<-hex.to.dec(bmid)-hex.to.dec(bmin)

rr2<-hex.to.dec(rmax)-hex.to.dec(rmid)
gr2<-hex.to.dec(gmax)-hex.to.dec(gmid)
br2<-hex.to.dec(bmax)-hex.to.dec(bmid)

if(x < 0 ){
xrr<-round((x)*rr1) + hex.to.dec(rmid)
xgr<-round((x)*gr1) + hex.to.dec(gmid) 
xbr<-round((x)*br1) + hex.to.dec(bmid)
}

if(x > 0 ){
xrr<-round((x)*rr2)+ hex.to.dec(rmid)
xgr<-round((x)*gr2)+ hex.to.dec(gmid)
xbr<-round((x)*br2)+ hex.to.dec(bmid)
}


if(x == 0 ){
xrr<-hex.to.dec(rmid)
xgr<-hex.to.dec(gmid)
xbr<-hex.to.dec(bmid)
}

rx<-dec.to.hex(xrr)
gx<-dec.to.hex(xgr)
bx<-dec.to.hex(xbr)


if(nchar(rx)==1) rx<-paste('0', rx ,sep='')
if(nchar(gx)==1) gx<-paste('0', gx ,sep='')
if(nchar(bx)==1) bx<-paste('0', bx ,sep='')


out<-paste('#', rx, gx, bx, sep='')
###########
}

if(is.na(x))  out<-''


return(out)
}
################################################################################


###################
if(length(x)==1){
if(!is.numeric(x))  x<-as.numeric(x)
if(abs) x<-abs(x)
  out<-hexr(x,hmin=hmin, hmid=hmid ,hmax=hmax,rangemy=rangemy)
  }
##################

if(length(x)>1)
{

if(is.matrix(x)){
if(!is.numeric(x))  x<-matrix(as.numeric(x),  dim(x)[1], dim(x)[2])
if(abs) x<-abs(x)
out<-  matrix(sapply(x,hexr, hmin=hmin, hmid=hmid ,hmax=hmax, rangemy=rangemy, simplify =TRUE), dim(x)[1], dim(x)[2])
}
  
if(is.vector(x)){
if(!is.numeric(x))  x<-as.numeric(x)
if(abs) x<-abs(x)
out<-  sapply(x,hexr, hmin=hmin, hmid=hmid  ,hmax=hmax, rangemy=rangemy, simplify =TRUE)
}
    
}

return(out)
}

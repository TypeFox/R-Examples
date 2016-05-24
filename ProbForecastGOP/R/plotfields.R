plotfields <-
function(field,x.lim,y.lim,country.outline="US",title){

if(missing(country.outline)){
country.outline <- "US"
}

# Input check
dim.f <- dim(field)

if(length(dim.f)==0){
stop("Invalid input for field")
}

if(length(dim.f)!=2){
stop("Invalid input for field")
}

if(dim.f[1]!=65 | dim.f[2]!=65){
stop("Invalid input for field")
}

dim.f2 <- dim.f^2

if(sum(is.numeric(field))==0){
stop("Invalid input for field")
}

## input check for the x vector
l.limx <- length(x.lim)
if(l.limx < 2){
stop("Invalid input for x.lim")
}

if(l.limx > 2){
stop("Invalid input for x.lim")
}

if(l.limx==2){
if(sum(is.numeric(x.lim)==rep(TRUE,2))<2){
stop("Invalid input for x.lim")
}
if(x.lim[1] >= x.lim[2]){ 
stop("Invalid input for x.lim")
}
}


## input check on the lims for the y vector 
l.limy <- length(y.lim)
if(l.limy < 2){
stop("Invalid input for y.lim")
}

if(l.limy > 2){
stop("Invalid input for y.lim")
}

if(l.limy==2){
if(sum(is.numeric(y.lim)==rep(TRUE,2))<2){
stop("Invalid input for y.lim")
}
if(y.lim[1] >= y.lim[2]){ 
stop("Invalid input for y.lim")
}
}

# input check on the country.outline field
if(is.character(country.outline)==FALSE){
stop("country.outline can only be equal to US, world or both")
}

if(country.outline!="US" & country.outline!="world" & country.outline!="both"){
stop("country.outline can only be equal to US, world or both")
}

# here the function starts
lims <- c(min(field,na.rm=TRUE),max(field,na.rm=TRUE))
size=65  
hor.crd <- seq(x.lim[1],x.lim[2],length = size)  
ver.crd <- seq(y.lim[1],y.lim[2],length = size)  
n.level <- 5
US.map <- 0
par(ask=FALSE)
image.plot(hor.crd, ver.crd, field, xlim=c(min(hor.crd), max(hor.crd)), ylim = c(min(ver.crd),max(ver.crd)),zlim=lims,legend.width=.015,legend.shrink=.8,
   main=title,xlab="",ylab="",offset = 0.02, col=rainbow(100,start=0,end=0.85)[100:1])

if(min(hor.crd) <= -124.7 & max(hor.crd) >= -124.7){
US.map <- 1
}

if(min(hor.crd) >= -124.7 & max(hor.crd) <= -67.1){
US.map <- 1
}

if(min(hor.crd) <= -67.1 & max(hor.crd) >= -67.1){
US.map <- 1
}

if(min(ver.crd) <= 25.2 & max(ver.crd) >= 25.2){
US.map <- 1 + US.map
}

if(min(ver.crd) >= 25.2 & max(ver.crd) <= 49.4){
US.map <- 1 + US.map
}

if(min(ver.crd) <= 49.4 & max(ver.crd) >= 49.4){
US.map <- 1 + US.map
}

if(US.map==2 & country.outline=="US"){
US(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
} 

if(US.map==2 & country.outline!="US"){
US(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)     
}  

if(US.map < 2 & country.outline=="world"){
world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
}

if(US.map < 2 & (country.outline=="US" || country.outline=="both")){
world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
print("The area delimited by the latitude and longitude specified is not 
  contained in the US")
}
}


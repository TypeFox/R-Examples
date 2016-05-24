a_cb_brackets <-
function(phi=1, ticks=0.5, type=1){
if(phi<0.01) phi <-0.01
n      <- 1000
d_start <- 1
d_end   <- 1
add_s   <- 0
add_e   <- 0
if(!is.null(ticks)){
ticks  <- unique(ticks)
tsigns <- sign(ticks)
o      <- order(abs(ticks))
ticks  <- abs(ticks)[o]
tsigns <- tsigns[o]
if(ticks[1]==0){
add_s <- 2
d_start<- -1
ticks<- ticks[-1]
tsigns <- tsigns[-1]
}
if(length(ticks)>0){
if(ticks[length(ticks)]==1){
add_e <- 2
d_end<- -1
ticks  <- ticks[-length(ticks)]
tsigns <- tsigns[-length(tsigns)]
}
}
}
nt     <- length(ticks)
np     <- 2+2*nt
rx     <- (1/np)*phi
md     <- min(diff(c(0, ticks, 1)))
if(md<(rx*2)) rx<- md/2
if(type==1) p <- -rev(exp(seq(0,5,length.out=round(n*rx))))
if(type==2) p <-  sqrt(seq(0,5,length.out=round(n*rx)))
if(type==3) p <- seq(0,1,length.out=round(n*rx))
if(type==4) p <- c(0,rep(1, round(n*rx)-1))
if(type==5) p <- -rev((seq(0,5,length.out=round(n*rx)))^2)
p      <- a_st(p)
pb     <- length(p)
sy <- c(p*d_start+add_s)
location <- pb
if(nt>0){
for(i in 1:nt){
on <- round(n*ticks[i])
add<- 2
if(tsigns[i]==-1) add<- 0
sy <- c(sy, rep(1, on-location-pb), tsigns[i]*(-rev(p)+add), tsigns[i]*(-p+add))
location <- on+pb
}
}
sy <- c(sy, rep(1, n-location-pb), rev(p)*d_end+add_e)
sy <- a_st(sy)
sx <- seq(0, 1, length.out=length(sy))
rbind(sx, sy)
}

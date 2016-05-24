grid.on <-
function(h=TRUE,v=TRUE,col="gray",nv=1,nh=1,lty=3){
is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  
abs(x - round(x)) < tol

if(!is.wholenumber(nh))stop("Error in 'nh'!")
if(!is.wholenumber(nv)) stop("Error in 'nv'!")
  
ve<-par()
n.col<-ve$xaxp
n.row<-ve$yaxp

#ve luoi ngang
if(h==TRUE){
for(i in -nh:(nh*n.row[3]+nh))
abline(h=n.row[1]+i*((n.row[2]-n.row[1])/(n.row[3]*nh)),col=col,lty=lty)
}

#ve luoi dung
if(v==TRUE){
for(i in -nv:(nv*n.col[3]+nv))
abline(v=n.col[1]+i*((n.col[2]-n.col[1])/(n.col[3]*nv)),col=col,lty=lty)
}

}

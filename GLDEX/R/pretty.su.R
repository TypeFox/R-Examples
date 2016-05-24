"pretty.su" <-
function(x,nint=5){
r<-range(x)
x<-seq(r[1],r[2],length=nint+1)
return(x)
}


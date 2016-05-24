residue <-
function(f,z0){
res<-paste("F<-function(",names(formals(f)),"){(",names(formals(f)),"-(",z0,"))*(",body2string(f),")}",collapse="",sep="")
eval(parse(text=res))
X<-delta<-NULL
dir<-seq(0,7,1)
for(i in 1:8){
X<-c(X,limit(F,z0,z0+exp(1i*pi*dir[i]/4)))
}
for(i in 1:20){if(length(unique(signif(abs(Re(X)),i)))>1) break};Re.digit<-i-1
for(j in 1:20){if(length(unique(signif(abs(Im(X)),j)))>1) break};Im.digit<-j-1
Re.sign<-sign(sum(sign(Re(X))))
Im.sign<-sign(sum(sign(Im(X))))
#result<-signif(abs(Re(X[1])),Re.digit)*Re.sign+signif(abs(Im(X[1])),Im.digit)*Im.sign*1i
result<-mean(abs(Im(c(X[1],X[5]))))*1i*Im.sign+mean(abs(Re(c(X[1],X[5]))))*Re.sign
return(result)
}

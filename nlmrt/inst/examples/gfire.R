# 3000 PRINT "GFIRE - FIRE PROTECTION MODEL - 870512"
model<-"yy ~ 10*b1*xx^b2 + 10*b3"
xx<-c(2.632, 1.216, 1.736, 0.388, 1.769, 0.476, 1.172, 1.546, 1.769, 0.978, 
       1.991, 0.954, 3.9, 2.606)
yy<-c(185, 127, 158, 53, 135, 61, 110, 139, 166, 118, 152, 102, 236, 198)

require(nlmrt)
ans<-nlsmnq(model, start=c(b1=1,b2=1, b3=1), data=data.frame(yy=yy, xx=xx), trace=TRUE, control=list(watch=TRUE))
print(ans)



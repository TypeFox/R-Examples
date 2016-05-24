signifCC<-
function(tableau){
X2<-chisq.test(tableau)$statistic
df<-chisq.test(tableau)$parameter
p.value<-chisq.test(tableau)$p.value
N<-sum(tableau)
k<-min(dim(tableau))-1
w<-sqrt(X2/N)
V<-w/sqrt(k)

p1<-paste("w=",round(w,3)," , V=",round(V,3),sep="")
p2<-paste("p-value=",signif(p.value),sep="")
p3<-"[Cohen] : petit (w=0.1), moyen (w=0.3), grand (w=0.5)\n"
p4<-paste("Note : X2=",round(X2,3)," , df=",df, " ,N=",N,sep="")
cat(p1,p3,p2,"\n",p4,sep="\n")
}

signifNC<-
function(Maov){
eta2<-summary.lm(Maov)$r.squared
s<-summary.lm(Maov)$fstatistic
p.value<-pf(q=s[1], df1=s[2], df2=s[3],lower.tail = FALSE)
p1<-paste("eta2=",round(eta2,4),sep="")
p2<-paste("p-value=",signif(p.value),sep="")
p3<-"[Cohen] : petit (eta2=0.01), moyen (eta2=0.05), grand (eta2=0.15)\n"
cat(p1,p3,p2,sep="\n")
}

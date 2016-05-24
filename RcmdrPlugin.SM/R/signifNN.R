signifNN<-
function(LM){
r2<-summary.lm(LM)$r.squared
s<-summary.lm(LM)$fstatistic
p.value<-pf(q=s[1], df1=s[2], df2=s[3],lower.tail = FALSE)
p1<-paste("r2=",round(r2,4),sep="")
p2<-paste("p-value=",signif(p.value),sep="")
p3<-"[Cohen] : petit (r2=0.01), moyen (r2=0.10), grand (r2=0.25)\n"
cat(p1,p3,p2,sep="\n")
}

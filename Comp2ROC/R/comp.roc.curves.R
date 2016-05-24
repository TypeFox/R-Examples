comp.roc.curves <-
function(result,ci.flag=FALSE,graph.flag=FALSE,nome) {

x<-result$diffareas
mask=c()
for(i in 1:length(x))
{
if (x[i]<0)  mask=c(mask,-1) else mask=c(mask,1)
}
cc=0
# funcao para permutation
soma <- function(x,d) {
      cc <<- cc+1
      if (cc == 1)
        return(sum(x[d]))
      else
      {
        perm <- mask*abs(x[d])
        return(sum(perm))
      }
}
# resampling
b=boot(x, soma, R=10000)
#print(b)
if (graph.flag==TRUE){
  # graficos
  plot(b)
  mtext(paste("PERMUTATION - MEDIAN :",nome),line=3)
}
# p-value unilateral
pp1=sum(b$t <= b$t0)/b$R
pp2=sum(b$t >= b$t0)/b$R
p1=sum(b$t <= 0)/b$R
p2=sum(b$t >= 0)/b$R
if  (p1>p2)  p.value=p2 else p.value=p1
p.value2=2*p.value
#cat(nome, "\n p-value (one-sided)= ",p.value,"\n")
#cat(nome, "\n p-value (two-sided)= ",p.value2,"\n")
if  (pp1>pp2)  pp.value=pp2 else pp.value=pp1
pp.value2=2*pp.value
#cat(nome, "\n pp-value (one-sided)= ",pp.value,"\n")
#cat(nome, "\n pp-value (two-sided)= ",pp.value2,"\n")
# intervalos de confianca
if (ci.flag==TRUE)
  ci=boot.ci(b,type=c("basic","perc"))
else
  ci=NA

answer=list(boot=b,pvalue=p.value2,pvalue2=pp.value2,ci=ci)
# grafico jackknife
#windows()
#jack.after.boot(b)
#mtext(paste("PERMUTATION - MEDIAN :",nome),line=3)
return(answer)
}

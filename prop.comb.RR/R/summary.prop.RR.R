summary.prop.RR <-
function(object, ...) {

 
    cat("", "\n")
    cat("      ","relative risk", "\n")
   
cat("","\n")

P=cbind(1:2,object$x,object$n,object$estimate)

colnames(P)=c("sample","x","n","prop")
rownames(P)=rep("",dim(P)[1])

cat("x: number of successes", "\n")
cat("n: number of trials", "\n")
cat("prop: proportion sample estimates", "\n")
cat("","\n")
print(P)

cat("","\n")


cat("realtive risk: RR=p2/p1","\n")
cat("estimated RR:",object$RR, "\n")

cat("","\n")

borrar = switch(object$alternative, 
          two.sided = "is not equal to", 
          greater = "is greater than", 
          less = "is less than" )

cat("alternative hypothesis: true relative risk R", borrar, object$rho, "\n")
cat("","\n")
cat("               ",round(100*object$conf.level),"percent confidence interval", "\n")
print(object$inference)
cat("","\n")
cat("Recommendation:","\n")
cat( object$recomen,"\n")
}
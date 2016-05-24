`reg.homog` <-
function(trt,x,y) {
sumx<-function(x)(nrow(x)-1)*var(x)
datos<-data.frame(trt,x,y)
sumxy<-by(datos[,-1],trt,function(x) sumx(x))
t.vec<-as.character(unique(datos$trt))
r.total <-nrow(datos)
n.trt<-length(t.vec)
sx<-0; sy<-0; sxy<-0;residual<-0
for ( i in 1:n.trt) {
a1<-data.frame(sumxy[t.vec[i]])
a2<-as.matrix(a1)
residual <- residual+a2[2,2]-a2[1,2]^2/ a2[1,1]
sx<- sx + a2[1,1]
sy<- sy + a2[2,2]
sxy<- sxy + a2[1,2]
}
# suma de las regresiones
B<-sy-sxy^2/sx 
# suma de los residuales
A<- residual
# diferencia de homogeneidad
diff<-B-A
# Prueba de la homogenidad de regresiones
gl.trt <- n.trt-1
gl.r <- r.total - 2* n.trt
f.cal <- ( (B-A)/gl.trt ) /  ( A/gl.r)
p.value <- 1- pf( f.cal, gl.trt, gl.r)
resp<-"homogeneity of regressions exists"
if(p.value <= 0.05 ) resp<-"homogeneity of regressions does not exists"
# Imprime resultados
cat("\nTest of Homogeneity of regressions\n\n")
cat("Total of simple regressions: ", n.trt ,"\n")
cat("Total of residual          : ", A, "\n")
cat("Difference for homogeneity : ", B-A, "\n\n")
cat("D.f. for the homogeneity   : ", gl.trt,"\n")
cat("Residual degrees of freedom: ", gl.r,"\n")
cat("F calculated value         : ", f.cal,"\n")
cat("P.value                    : ", p.value,"\n\n")
cat("Criterion                  : ", resp,"\n\n")
output<-list(regressions=n.trt,residual=A,Difference=B,DF.homgeneity=gl.trt,DF.Residual=gl.r,F.value=f.cal,P.value=p.value,Criterion=resp)
return(output)
}


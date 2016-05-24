iconfianza.percentil.boot <-
function(datos,q=0.50,nivel=0.95,ic=T,tipo.boot="norm",iteraciones.boot=10000,colas=2){
	# Esto lo hago para que luego el !all me funcione, dos valores aparecen como 100, pero realmente son distintos, y el comparador logico lo detecta, pero el boot NO!!!
	x<-round(datos[!is.na(datos)],6)
  n<-length(x)
	if (n!=0){
  	quantile.fun<-function(d,i){
  		a<-quantile(d[i],probs=q)
  		n<-length(i)
  		b<-(n-1)*var(d[i])/n^2
  		return(c(a,b))
  	}
  	if (!ic){
  		return(rep(quantile(datos,probs=q),3))
  	} else {
  		datos.boot<-boot(x,quantile.fun,iteraciones.boot)
  		if (!all(min.fix.na(x)==x,na.rm=T)){
  			inter.boot<-boot.ci(datos.boot,type=c("norm","basic","perc","stud"),conf=((3-colas)*nivel-2+colas))
  			if (tipo.boot=="norm" || tipo.boot=="normal") return(c(inter.boot$normal[2],inter.boot$t0[1],inter.boot$normal[3]))
  			else if (tipo.boot=="basic") return(c(inter.boot$basic[4],inter.boot$t0[1],inter.boot$basic[5]))
  			else if (tipo.boot=="stud" || tipo.boot=="student") return(c(inter.boot$student[4],inter.boot$t0[1],inter.boot$student[5]))
  			else if (tipo.boot=="perc" || tipo.boot=="percent") return(c(inter.boot$percent[4],inter.boot$t0[1],inter.boot$percent[5]))
  			#else if (tipo.boot=="bca") return(c(inter.boot$bca[4],inter.boot$t0[1],inter.boot$bca[5]))
  			#else return(rep(inter.boot$t0[1],3))
  			else return(c(NA,inter.boot$t0[1],NA))
  		}else{
  		  #return(rep(datos.boot$t0[1],3))
  		  return(c(NA,inter.boot$t0[1],NA))
  		}
  	}
	} else return(rep(NA,3))
}

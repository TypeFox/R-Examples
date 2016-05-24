
#                       F U N C I O N E S
#
#####################################################################
## FUNCION PARA IDENTIFICAR PUNTOS ATIPICOS EN VARIABLES EN DONDE  
## HAY CAMPOS CON ALGUN VALOR NA "NO APLICA" 
## LA FUNCIÓN REALIZA COMBINACIONES DE VARIABLES  
## PARÁMETROS DE ENTRADA
## datos  matriz con los datos para identificar puntos extremos 
##		incluyendo campos llave, estos deben estar al inicio de 
##		los registros.
## ncllav	número de campos llave
## SALIDA	
## regID	Campos de identificación del registro
## valores 	Valor de las variables identificadas en grupo como atípicas
## norm.estan  valor de la normal estándar del punto respecto a la variable
## vars.ext  variables atípicas univariadas.
#####################################################################

PEconNA<-function(datos,ncllav,nsig){
  nas<-sum(as.numeric(is.na(datos)))
  stopifnot(dim(datos)[1]>20,dim(datos)[2]>ncllav+2,nas>0)
  TAM<-dim(datos)
  TAM[2]<-TAM[2]-ncllav
  datos1<-matrix(rep(0,TAM[1]*TAM[2]),ncol=TAM[2])
  pvalores<-rep(0,TAM[2])
  Tipo_transforma<-rep(NA,TAM[2])
  Normalidad<-rep(NA,TAM[2])
  for(i in 1:TAM[2]){
    paso <- MejorTr(datos[,i+ncllav])
    datos1[,i] <- paso$vartra
    pvalores[i] <- paso$pvalor
    Tipo_transforma[i] <- paso$trans
    ifelse(pvalores[i]>=.05,Normalidad[i]<-"Excelente",
           ifelse(pvalores[i]<.05&pvalores[i]>=.01,Normalidad[i]<-"Muy buena",
                  ifelse(pvalores[i]<.01&pvalores[i]>=.001,Normalidad[i]<-"Buena",
                         ifelse(pvalores[i]<.001&pvalores[i]>1e-5,Normalidad[i]<-"Regular",
                                Normalidad[i]<-"No se normaliz?"))))
  }
  Trans<-cbind(names(datos[ncllav+1:TAM[2]]),Normalidad,Tipo_transforma,pvalores)
  ##############################################################
  ##										 #
  ## Combinaci?n de variables en donde cada registro tenga     #
  ## informaci?n diferente de cero.   n = n?m de variables     #
  ##										 #
  ##############################################################
  comple<-matrix(rep(0,3),ncol=3)
  NC<-0
  Ncombina<-0
  for(i in seq(dim(datos1)[2],3, by = -1)){   
    combina<-combn(dim(datos1)[2],i)          # combinaciones n de i 
    NC<-NC+Ncombina					
    Ncombina<-choose(dim(datos1)[2],i)		# n?mero de combinaciones
    compleA<-matrix(rep(0,Ncombina*3),ncol=3) 
    comple<-rbind(comple,compleA)
    for(j in 1:dim(combina)[2]){
      datos1b<-datos1[,combina[,j]]
      # calcula NA por registro (rengl?n)
      SNA<-rep(0,dim(datos1b)[1])
      for(k in 1:dim(datos1)[1])
        SNA[k]<-sum(as.numeric(is.na(datos1b[k,])))
      comple[NC+j,]<-c(i,j,length(SNA[SNA==0]))
    }
  }
  # Selecci?n de combinaciones que contienen mas de 20 registros con informaci?n sin NA
  ComValid<-comple[comple[,3]>20,]
  PE<-rep(0,TAM[1])
  # identifica PE para cada combinaci?n v?lida
  for(i in 1:dim(ComValid)[1]){
    combOK<-combn(dim(datos1)[2],ComValid[i,1])[,ComValid[i,2]]
    datos1a<-datos1[,combOK]
    # separa los datos sin NA para las variables seleccionadas
    MV<-var(datos1a,na.rm=TRUE)
    invMV<-solve(MV)
    distM <- rep(0,TAM[1])
    medias<-apply(datos1a,2,mean,na.rm=TRUE) 
    # Genera distancias de Mahalanobis
    for (j in 1:TAM[1])
      distM[j]=as.numeric(datos1a[j,]-medias)%*%(invMV)%*%as.numeric(t(datos1a[j,]-medias))
    #Compara con el valor de la Ji-cuadrada al 95% dos colas
    mal<-length(na.omit(distM[distM>qchisq(nsig,dim(datos1a)[2]-1)]))
    if (mal>0)
      PE[distM>qchisq(nsig,dim(datos1a)[2]-1)]<-1
  }
  NPE<-length(PE[PE==1])	
  if (NPE == 0)
    stop("No hay puntos extremos")
  medias<-apply(datos1,2,mean,na.rm=TRUE) 
  desvestan<-apply(datos1,2,sd,na.rm=TRUE)
  pextrem<-datos1[PE==1,]
  dpextrem<-datos[PE==1,c(1:ncllav)]
  pextori<-as.matrix(datos[PE==1,c((ncllav+1):(TAM[2]+ncllav))])
  tam<-dim(pextori)
  print(tam[1])
  if(tam[1]==1){
    estan<-as.numeric((pextrem-medias)/desvestan)
    idvarne<-names(datos[1,c((ncllav+1):(tam[2]+ncllav))])[abs(estan)>qnorm(nsig)]
  }
  if(tam[1]>1){
    estan<-matrix(rep(0,tam[1]*tam[2]),ncol=tam[2])
    idvarne<-matrix(rep(NA,tam[1]*tam[2]),ncol=tam[2])
    for(j in 1:tam[1]){
      estan[j,]<-as.numeric((pextrem[j,]-medias)/desvestan)
      NN<-sum(as.numeric(na.omit(abs(estan[j,])>qnorm(nsig))))
      if(NN>0)
        idvarne[j,1:NN]<-na.omit(names(datos[j,c((ncllav+1):(tam[2]+ncllav))])[abs(estan[j,])>qnorm(nsig)])
      else
        idvarne[j,]<-rep(NA,tam[2])
    }
  }
  return(list(Transforma=Trans,regID=dpextrem,valores=pextori,norm.estan=estan,vars.ext=idvarne))
}

###########################################################################
##  Funci?n para encontrar la mejor transformaci?n entre raiz cuad. y log #
##  para encontrar normalidad								  # 
###########################################################################
MejorTr<-function(variable){
	stopifnot(length(variable)>0)
	if(length(na.omit(variable[variable<0]))>0)
		variable<-variable - min(variable) + 0.25
	if(length(na.omit(variable[variable==0]))>0){
		varT<-sqrt(variable)
		PV<-shapiro.test(varT)$p.value
		return(list(vartra=varT,pvalor=PV,trans="Ra?z cuad."))}
	prue<-rep(0,3)
	prue[1]<-shapiro.test(variable)$p.value
	prue[2]<-shapiro.test(sqrt(variable))$p.value
	prue[3]<-shapiro.test(log(variable))$p.value
	o<-order(unlist(prue))
	tr<-c("ninguno","raiz_cuad","log")
	mejor<-tr[o][3]
	varT<-switch(mejor,ninguno=variable,raiz_cuad=sqrt(variable),log=log(variable))	
	PV<-max(unlist(prue))	
	return(list(vartra=varT,pvalor=PV,trans=mejor))}	


##############################################################################
##  Grafica la funci?n de densidad normal e identifica los valores de las
##  variables de los registros detectados como extremos con la funci?n PEconNA.
##  Par?metros de entrada: Resultados dela funci?n PEconNA,
##				   n?mero del punto extremo encontrado en el resultado 
##############################################################################
graNormPE<-function(rpe,npe){
  stopifnot(npe>0)
  numPE<-dim(rpe$norm.estan)[1]
  if(npe>numPE)
    stop("El punto extremo dado no existe")
  texto<-paste('Curva normal estandar con los valores del registro',
               as.vector(names(rpe$regID[1,]))[1],'=',rpe$regID[npe,1],
               as.vector(names(rpe$regID[1,]))[2],'=',rpe$regID[npe,2],sep=' ')
  numvars<-dim(rpe$norm.estan)[2]
  tintas<-c('black','red','blue','magenta','green4','brown','orange3','purple','cyan3') 
  plot(0,0,xlim=c(-4,4),ylim=c(0,.45),col="white",xlab='x',ylab='f(x)',
       main=texto)
  curve(dnorm(x,0,1),add=T,col='red',lwd=3)
  segments(c(-2,2),c(-1,-1),c(-2,2),c(.055,.055),col='black',lty=3,lwd=4)
  for(i in 1:numvars){
    abline(v=rpe$norm.estan[npe,i],col=tintas[i],lty=3,lwd=2)
    text(rpe$norm.estan[npe,i],0.05+0.03*i,col=tintas[i],
         paste(as.vector(names(rpe$valores[1,]))[i],'=',round(rpe$valores[npe,i],2)),cex=1.5)}
}

PonNA <- function(dat,n){ #values$a<-dat
  d<-dim(dat)
  dd<-length(n)
  for (i in 1:d[1])
  {
    for(j in 1:dd)
    if (dat[i,n[j]]==0)
      dat[i,n[j]]<-NA
  }
  return(dat)}
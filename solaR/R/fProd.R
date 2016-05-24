fProd<-function(inclin, 
                module=list(), 
                generator=list(), 
                inverter=list(),
                effSys=list()
                ){
  
  stopifnot(is.list(module),
            is.list(generator),
            is.list(inverter),
            is.list(effSys)
            )
###Preparo los datos	
  if (class(inclin)=='Gef') {
    indInclin=indexI(inclin)
    Gef=coredata(inclin@GefI$Gef)
    Ta=coredata(inclin@Ta)
  } else {
    if (class(inclin)=='zoo') {
      indInclin=index(inclin)
      Gef=coredata(inclin$Gef)
      Ta=coredata(inclin$Ta)
    } else {
      Gef=inclin$Gef
      Ta=inclin$Ta
    }
  }
  
###Listas de parámetros
  module.default=list(Vocn=57.6,Iscn=4.7,Vmn=46.08,Imn=4.35,Ncs=96,Ncp=1,CoefVT=0.0023, TONC=47)
  module=modifyList(module.default, module)
  
  generator.default=list(Nms=12,Nmp=11)
  generator=modifyList(generator.default, generator)

  inverter.default=list(Ki = c(0.01,0.025,0.05),Pinv=25000,Vmin=420, Vmax=750, Gumb=20)
  inverter=modifyList(inverter.default, inverter)

  effSys.default=list(ModQual=3,ModDisp=2,OhmDC=1.5,OhmAC=1.5,MPP=1,TrafoMT=1,Disp=0.5)
  effSys=modifyList(effSys.default, effSys)
  
###Constantes
  Gstc=1000
  Ct=(module$TONC-20)/800
  Vtn=0.025*(273+25)/300
  m=1.3
	
###Cálculo de Tc
  Tc=Ta+Ct*Gef

###Método de Ruiz para el cálculo de la tensión y corriente MPP de una célula en condiciones no estándar
###Calculos para UNA CÉLULA
  Rs=with(module,(Vocn/Ncs-Vmn/Ncs+m*Vtn*log(1-Imn/Iscn))/(Imn/Ncp))
  rs=with(module,Rs/((Vocn/Ncs)/(Iscn/Ncp)))
	
  Vt=0.025*(Tc+273)/300
  Voc_c=with(module,Vocn/Ncs-CoefVT*(Tc-25))
  Isc_c=with(module,Iscn/Ncp*Gef/Gstc)

  koc=Voc_c/(m*Vt)
  Dm0=(koc-1)/(koc-log(koc))
  Dm=Dm0+2*rs*Dm0^2
	
  Impp_c=Isc_c*(1-Dm/koc)
  Vmpp_c=Voc_c*(1-log(koc/Dm)/koc-rs*(1-Dm/koc))
	
###Corriente, Tensión y Potencia del generator
  Voc=with(module,Voc_c*Ncs*generator$Nms)
  Isc=with(module,Isc_c*Ncp*generator$Nmp)
  Impp=with(module,Impp_c*Ncp*generator$Nmp)
  Vmpp=with(module,Vmpp_c*Ncs*generator$Nms)
  Pmpp=Impp*Vmpp
	
  ##Cálculo de corriente para tensión diferente al MPP
  ##cuando el inverter limita por tensión fuera de rango
  f<-function(i,v,koc){
    vp=v+i*rs
    Is=1/(1-exp(-koc*(1-rs)))
    result=i-(1-Is*(exp(-koc*(1-vp))-exp(-koc*(1-rs))))}
  raiz<-function(M){
    v=M[1]
    koc=M[2]
    if (is.na(koc)) {
      result<-NA
    } else {
      result<-uniroot(f,c(0,1),v=v,koc=koc)$root}
  }
	
  Vdc=Vmpp
  Idc=Impp

  ##¿Está por debajo de la mínima tensión del inverter?
  if (any(Vmpp<inverter$Vmin,na.rm=1)){
    indMIN=which(Vmpp<inverter$Vmin)
    VocMIN=Voc[indMIN]
    kocMIN=koc[indMIN]
    vmin=inverter$Vmin/VocMIN
                                        #v debe estar entre 0 y 1
    vmin[vmin<0]=0
    vmin[vmin>1]=1
    imin=apply(cbind(vmin,kocMIN),1,FUN=raiz)
    IscMIN=Isc[indMIN]
    Idc[indMIN]=with(module,imin*IscMIN)
    Vdc[indMIN]=inverter$Vmin
    warning('Minimum MPP voltage of the inverter has been reached')}

  ##¿Está por encima de la máxima tensión del inverter?
  if (any(Vmpp>inverter$Vmax,na.rm=1)){
    indMAX=which(Vmpp>inverter$Vmax)
    VocMAX=Voc[indMAX]
    kocMAX=koc[indMAX]
    vmax=inverter$Vmax/VocMAX
                                        #v debe estar entre 0 y 1
    vmax[vmax<0]=0
    vmax[vmax>1]=1
    imax=apply(cbind(vmax,kocMAX),1,FUN=raiz)
    IscMAX=Isc[indMAX]
    Idc[indMAX]=with(module,imax*IscMAX)
    Vdc[indMAX]=inverter$Vmax
    warning('Maximum MPP voltage of the inverter has been reached')}
	
###Potencia DC normalizada al inverter
  PdcN=with(inverter,(Idc*Vdc)/Pinv*(1-effSys$ModQual/100)*(1-effSys$ModDisp/100)*(1-effSys$MPP/100)*(1-effSys$OhmDC/100)) 

###Potencia AC normalizada al inverter
  if (is.matrix(inverter$Ki)){ #Ki es una matriz de nueve coeficientes-->dependencia con tensión
    VP <- cbind(Vdc, PdcN)
    PacN <- apply(VP, 1, solvePac, inverter$Ki)
  } else { #Ki es un vector de tres coeficientes-->sin dependencia con la tensión
    PacN=with(inverter,{
      A=Ki[3]
      B=Ki[2]+1
      C=Ki[1]-(PdcN)
      result=(-B+sqrt(B^2-4*A*C))/(2*A)
    })
  }
  EffI <- PacN/PdcN
  pacNeg <- PacN <= 0
  PacN[pacNeg] <- PdcN[pacNeg] <- EffI[pacNeg] <- 0

	
###Potencia AC y DC sin la normalización
  Pac=with(inverter,PacN*Pinv*(Gef>Gumb)*(1-effSys$OhmAC/100)*(1-effSys$TrafoMT/100)*(1-effSys$Disp/100))
  Pdc=PdcN*inverter$Pinv*(Pac>0)
	
                                        
###Resultado
  Pg=generator$Nms*module$Vmn*generator$Nmp*module$Imn
  generator$Pg=Pg
  
  resProd <-data.frame(Tc,Voc,Isc,Vmpp,Impp,Vdc,Idc,Pac,Pdc,EffI)
  if (class(inclin) %in% c('Gef', 'zoo')) {
    result<-zoo(resProd, order.by=indInclin)
    attr(result, 'generator')=generator
    attr(result, 'module')=module
    attr(result, 'inverter')=inverter
    attr(result, 'effSys')=effSys
    return(result)
  } else {
    result<-cbind(inclin, resProd)
    return(result)
  }
}





prodPVPS<-function(lat, 
                   modeTrk='fixed', 
                   modeRad='prom', 
                   dataRad,
                   sample='hour',
                   keep.night=TRUE,
                   sunGeometry='michalsky',
                   corr, f,
                   betaLim=90, beta=abs(lat)-10, alfa=0,
                   iS=2, alb=0.2, horizBright=TRUE, HCPV=FALSE,
                   pump , H, 
                   Pg, converter= list(), #Pnom=Pg, Ki=c(0.01,0.025,0.05)),
                   effSys=list(),
                   ...){
    
  stopifnot(is.list(converter),
            is.list(effSys))

  if (modeRad!='prev'){                 #No utilizamos un cálculo prev

    radEf<-calcGef(lat=lat, modeTrk=modeTrk, modeRad=modeRad,
                   dataRad=dataRad,
                   sample=sample, keep.night=keep.night,
                   sunGeometry=sunGeometry,
                   corr=corr, f=f,
                   betaLim=betaLim, beta=beta, alfa=alfa,
                   iS=iS, alb=alb, horizBright=horizBright, HCPV=HCPV,
                   modeShd='', ...)
		
  } else { #Utilizamos un cálculo previo de calcG0, calcGef o prodSFCR
    stopifnot(class(dataRad) %in% c('G0', 'Gef', 'ProdPVPS'))
    radEf <- switch(class(dataRad),
                    G0=calcGef(lat=lat, 
                      modeTrk=modeTrk, modeRad='prev',
                      dataRad=dataRad,
                      betaLim=betaLim, beta=beta, alfa=alfa,
                      iS=iS, alb=alb, horizBright=horizBright, HCPV=HCPV,
                      modeShd='', ...),
                    Gef=dataRad,
                    ProdPVPS=as(dataRad, 'Gef')
                    )
  }


###Producción eléctrica
  converter.default=list(Ki = c(0.01,0.025,0.05), Pnom=Pg)
  converter=modifyList(converter.default, converter)

  effSys.default=list(ModQual=3,ModDisp=2,OhmDC=1.5,OhmAC=1.5,MPP=1,TrafoMT=1,Disp=0.5)
  effSys=modifyList(effSys.default, effSys)

  TONC=47
  Ct=(TONC-20)/800
  lambda=0.0045
  Gef=coredata(radEf@GefI$Gef)
  aman=coredata(radEf@solI$aman)
  Ta=coredata(radEf@Ta)
  
  Tc=Ta+Ct*Gef
  Pdc=Pg*Gef/1000*(1-lambda*(Tc-25))
  Pdc[is.na(Pdc)]=0    #Necesario para las funciones que entrega fPump
  PdcN=with(effSys,
    Pdc/converter$Pnom*(1-ModQual/100)*(1-ModDisp/100)*(1-OhmDC/100)
    )
  PacN=with(converter,{
    A=Ki[3]
    B=Ki[2]+1
    C=Ki[1]-(PdcN)
    ##Potencia AC normalizada al inversor
    result=(-B+sqrt(B^2-4*A*C))/(2*A)
  })
  PacN[PacN<0]<-0
	
  Pac=with(converter,
    PacN*Pnom*(1-effSys$OhmAC/100))
  Pdc=PdcN*converter$Pnom*(Pac>0)
	
		
###Bomba
  fun<-fPump(pump=pump, H=H)
  ##Limito la potencia al rango de funcionamiento de la bomba
  rango=with(fun,Pac>=lim[1] & Pac<=lim[2]) 
  Pac[!rango]<-0
  Pdc[!rango]<-0
  prodI=data.frame(Pac=Pac,Pdc=Pdc,Q=0,Pb=0,Ph=0,f=0)	
  prodI=within(prodI,{
    Q[rango]<-fun$fQ(Pac[rango])
    Pb[rango]<-fun$fPb(Pac[rango])
    Ph[rango]<-fun$fPh(Pac[rango])
    f[rango]<-fun$fFreq(Pac[rango])
    etam=Pb/Pac
    etab=Ph/Pb
  })
  ##Pdc[!aman]<-NA
  ##Pac[!aman]<-NA
  prodI[!aman,]<-NA
  prodI<-zoo(prodI, indexI(radEf))

###Cálculo de valores diarios, mensuales y anuales
  ##Cálculo de valores diarios, mensuales y anuales
  ##=======================================
  DayOfMonth=c(31,28,31,30,31,30,31,31,30,31,30,31) ###OJO
    
  if (radEf@type=='prom') {
    prodDm=aggregate(prodI[,c('Pac', 'Q')],
      by=as.yearmon, FUN=P2E, radEf@sample) 
    names(prodDm)=c('Eac', 'Qd')
    prodDm$Eac=prodDm$Eac/1000          #kWh
    prodDm$Yf=prodDm$Eac/(Pg/1000)      #kWh/kWp
    
    prodD=prodDm
    prodD$Eac=prodD$Eac*1000         #Wh
    index(prodD) <- indexD(radEf)    ##para que sea compatible con G0D

    prody=zoo(t(colSums(prodDm*DayOfMonth)),
      unique(year(index(prodDm))))
  } else {
    prodD=aggregate(prodI[,c('Pac', 'Q')],
      by=truncDay, FUN=P2E, radEf@sample) #Wh
    names(prodD)=c('Eac', 'Qd')
    prodD$Yf=prodD$Eac/Pg
    
    prodDm=aggregate(prodD, by=as.yearmon, mean, na.rm=1)
    prody=aggregate(prodD, by=year, sum, na.rm=1)

    prodDm$Eac=prodDm$Eac/1000
    prody$Eac=prody$Eac/1000

  }
  
  result <- new('ProdPVPS',
                radEf,                  #contains 'Gef'
                prodD=prodD,
                prodDm=prodDm,
                prody=prody,
                prodI=prodI,
                pump=pump,
                H=H,
                Pg=Pg,
                converter=converter,
                effSys=effSys
                )
}
	

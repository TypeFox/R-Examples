optimShd<-function(lat,
                   modeTrk='fixed', 
                   modeRad='prom', 
                   dataRad,
                   sample='hour',
                   keep.night=TRUE,
                   sunGeometry='michalsky',
                   betaLim=90, beta=abs(lat)-10, alfa=0,
                   iS=2, alb=0.2, HCPV=FALSE,
                   module=list(), 
                   generator=list(),
                   inverter=list(), 
                   effSys=list(), 
                   modeShd='',    
                   struct=list(), 
                   distances=data.frame(),
                   res=2,     #resolución, separación entre distancias
                   prog=TRUE){          #Pinto barra de progreso
		
  if (('bt' %in% modeShd) & (modeTrk!='horiz')) {
    modeShd[which(modeShd=='bt')]='area'
    warning('backtracking is only implemented for modeTrk=horiz')}

  ##Guardo argumentos de la función para utilizar después

  listArgs<-list(lat=lat, modeTrk=modeTrk, modeRad=modeRad,
                 dataRad=dataRad,
                 sample=sample, keep.night=keep.night,
                 sunGeometry=sunGeometry,
                 betaLim=betaLim, beta=beta, alfa=alfa,
                 iS=iS, alb=alb, HCPV=HCPV,
                 module=module, generator=generator,
                 inverter=inverter, effSys=effSys,
                 modeShd=modeShd, struct=struct, distances=data.frame(Lew=NA, Lns=NA, D=NA))
  
  
  ##Creo red en la que haré los cálculos
  Red=switch(modeTrk,
    horiz=with(distances,
      data.frame(Lew=seq(Lew[1],Lew[2],by=res),
                 H=0)),
    two=with(distances,
      expand.grid(Lew=seq(Lew[1],Lew[2],by=res),
                  Lns=seq(Lns[1],Lns[2],by=res),
                  H=0)),
    fixed=with(distances,
      data.frame(D=seq(D[1],D[2],by=res),
                 H=0))
    )
  
  casos<-dim(Red)[1]               #Número de posibilidades a estudiar

  ##Preparo la barra de progreso
  if (prog) {pb <- txtProgressBar(min = 0, max = casos+1, style = 3)
             setTxtProgressBar(pb, 0)}
	
###Cálculos	
  ##Referencia: Sin sombras	
  listArgs0 <- modifyList(listArgs,
                          list(modeShd='', struct=NULL, distances=NULL) )
  Prod0<-do.call(prodGCPV, listArgs0)
  YfAnual0=mean(Prod0@prody$Yf)   #Utilizo mean por si hay varios años
  if (prog) {setTxtProgressBar(pb, 1)}
	
  ##Empieza el bucle
  
  ##Creo un vector vacío de la misma longitud que los casos a estudiar
  YfAnual<-numeric(casos) 

  BT=('bt' %in% modeShd)
  if (BT) { ##Hay backtracking, luego debo partir de radiación horizontal
    RadBT <- as(Prod0, 'G0')
    for (i in seq_len(casos)){
      listArgsBT <- modifyList(listArgs,
                               list(modeRad='prev', dataRad=RadBT,
                                    distances=Red[i,]))
      prod.i <- do.call(prodGCPV, listArgsBT)
      YfAnual[i]=mean(prod.i@prody$Yf)
      if (prog) {setTxtProgressBar(pb, i+1)}
    }
  } else {
    prom=('prom' %in% modeShd)
    for (i in seq_len(casos)){
      Gef0=as(Prod0, 'Gef')
      GefShd=calcShd(Gef0, modeTrk=modeTrk, modeShd=modeShd,
        struct=struct, distances=Red[i,])
      listArgsShd <- modifyList(listArgs,
                                list(modeRad='prev', dataRad=GefShd)
                                )
      prod.i <- do.call(prodGCPV, listArgsShd)
      YfAnual[i]=mean(prod.i@prody$Yf)
      if (prog) {setTxtProgressBar(pb, i+1)}
    }
  }
  if (prog) {close(pb)}


###Entrego resultados
  FS=1-YfAnual/YfAnual0
  GRR=switch(modeTrk,
    two=with(Red,Lew*Lns)/with(struct,L*W),
    fixed=Red$D/struct$L,
    horiz=Red$Lew/struct$L)
  SombraDF=cbind(Red,GRR=GRR,FS=FS,Yf=YfAnual)
  FS.loess=switch(modeTrk,
    two=loess(FS~Lew*Lns,data=SombraDF),
    horiz=loess(FS~Lew,data=SombraDF),
    fixed=loess(FS~D,data=SombraDF))
  Yf.loess=switch(modeTrk,
    two=loess(Yf~Lew*Lns,data=SombraDF),
    horiz=loess(Yf~Lew,data=SombraDF),
    fixed=loess(Yf~D,data=SombraDF))
  result <- new('Shade',
                Prod0, ##contains ProdGCPV
                FS=FS,
                GRR=GRR,
                Yf=YfAnual,
                FS.loess=FS.loess,
                Yf.loess=Yf.loess,
                modeShd=modeShd,
                struct=struct,
                distances=Red,
                res=res
                )
  result
}

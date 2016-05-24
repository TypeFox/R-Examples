fTheta<-function(sol, beta, alfa=0, modeTrk='fixed', betaLim=90, 
                 BT=FALSE, struct, dist)
{

    stopifnot(modeTrk %in% c('two','horiz','fixed'))
    if (!missing(struct)) {stopifnot(is.list(struct))}
    if (!missing(dist)) {stopifnot(is.data.frame(dist))}

    betaLim=d2r(betaLim)
    lat=getLat(sol, 'rad')
    signLat=ifelse(sign(lat)==0, 1, sign(lat)) ##Cuando lat=0, sign(lat)=0. Lo cambio a sign(lat)=1

    solI<-as.data.frameI(sol, complete=TRUE, day=TRUE)
    AlS=solI$AlS
    AzS=solI$AzS
    decl=solI$decl
    w<-solI$w

    aman<-solI$aman

    Beta<-switch(modeTrk,
                 two = {Beta2x=pi/2-AlS
                                        #if (BT==TRUE) {Beta=	
                     Beta=Beta2x+(betaLim-Beta2x)*(Beta2x>betaLim)},
                 fixed = rep(d2r(beta), length(w)), 
                 horiz={BetaHoriz0=atan(abs(sin(AzS)/tan(AlS)))
                     if (BT){lew=dist$Lew/struct$L
                         Longitud=lew*cos(BetaHoriz0)
                         Cond=(Longitud>=1)
                         Longitud[Cond]=1
                                        #Cuando Cond==TRUE Longitud=1 
                                        #y por tanto asin(Longitud)=pi/2,
                                        #de forma que BetaHoriz=BetaHoriz0
                         BetaHoriz=BetaHoriz0+asin(Longitud)-pi/2
                                        #=ifelse(Cond, 
                                        #BetaHoriz0,#No hay sombra
                                        #BetaHoriz0+asin(Longitud)-pi/2)
                     } else {
                         BetaHoriz=BetaHoriz0
                         rm(BetaHoriz0)}
                     Beta=ifelse(BetaHoriz>betaLim,betaLim,BetaHoriz)}
                 )
    is.na(Beta) <- (!aman)

    Alfa<-switch(modeTrk,
                 two = AzS,
                 fixed = rep(d2r(alfa), length(w)),
                 horiz=pi/2*sign(AzS))
    is.na(Alfa) <- (!aman)

    cosTheta<-switch(modeTrk,
                     two=cos(Beta-(pi/2-AlS)),
                     horiz={
                         t1=sin(decl)*sin(lat)*cos(Beta)      
                         t2=cos(decl)*cos(w)*cos(lat)*cos(Beta)   
                         t3=cos(decl)*abs(sin(w))*sin(Beta)   
                         cosTheta=t1+t2+t3
                         rm(t1,t2,t3)
                         cosTheta
                     },
                     fixed={
                         t1=sin(decl)*sin(lat)*cos(Beta)      
                         t2=-signLat*sin(decl)*cos(lat)*sin(Beta)*cos(Alfa) 
                         t3=cos(decl)*cos(w)*cos(lat)*cos(Beta)   
                         t4=signLat*cos(decl)*cos(w)*sin(lat)*sin(Beta)*cos(Alfa) 
                         t5=cos(decl)*sin(w)*sin(Alfa)*sin(Beta)   
                         cosTheta=t1+t2+t3+t4+t5
                         rm(t1,t2,t3,t4,t5)
                         cosTheta
                     }
                     )
    is.na(cosTheta) <- (!aman)
    cosTheta=cosTheta*(cosTheta>0) #cuando cosTheta<0, Theta es mayor de 90º, y por tanto el Sol está detras del panel.
    
    result<-zoo(data.frame(Beta, Alfa, cosTheta), order.by=indexI(sol))   
}

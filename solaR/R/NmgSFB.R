## NmgPVPS: no visible binding for global variable ‘Pnom’
## NmgPVPS: no visible binding for global variable ‘group.value’

if(getRversion() >= "2.15.1") globalVariables(c('Pnom', 'group.value'))

NmgPVPS <- function(pump, Pg, H, Gd, Ta=30,
                    lambda=0.0045, TONC=47,
                    eta=0.95, Gmax=1200, t0=6, Nm=6,
                    title='', theme=custom.theme.2()){

    ##Construyo el dia tipo mediante procedimiento IEC
    t=seq(-t0,t0,l=2*t0*Nm);
    d=Gd/(Gmax*2*t0)
    s=(d*pi/2-1)/(1-pi/4)
    G=Gmax*cos(t/t0*pi/2)*(1+s*(1-cos(t/t0*pi/2)))
    G[G<0]<-0
    G=G/(sum(G,na.rm=1)/Nm)*Gd
    Red<-expand.grid(G=G,Pnom=Pg,H=H,Ta=Ta)
    Red<-within(Red,{Tcm<-Ta+G*(TONC-20)/800
                     Pdc=Pnom*G/1000*(1-lambda*(Tcm-25)) #Potencia DC disponible
                     Pac=Pdc*eta})      #Rendimiento del inversor

    res=cbind(Red,Q=0)

    for (i in seq_along(H)){
        fun=fPump(pump, H[i])
        Cond=res$H==H[i]
        x=res$Pac[Cond]
        z=res$Pdc[Cond]
        rango=with(fun,x>=lim[1] & x<=lim[2]) #Limito la potencia al rango de funcionamiento de la bomba
        x[!rango]<-0
        z[!rango]<-0
        y=res$Q[Cond]
        y[rango]<-fun$fQ(x[rango])
        res$Q[Cond]=y
        res$Pac[Cond]=x
        res$Pdc[Cond]=z
    }

    resumen<-aggregate(res[c("G","Pdc","Q")],res[c('Pnom','H')],
                       FUN=function(x)sum(x,na.rm=1)/Nm)

    param=list(pump=pump, Pg=Pg, H=H, Gd=Gd, Ta=Ta,
    lambda=lambda, TONC=TONC, eta=eta, Gmax=Gmax, t0=t0, Nm=Nm)


###Abaco con los ejes X comunes

    ##Compruebo si tengo disponible el paquete lattice, que debiera haber sido cargado en .First.lib
    lattice.disp<-("lattice" %in% .packages())
    latticeExtra.disp<-("latticeExtra" %in% .packages())
    if (lattice.disp && latticeExtra.disp){
        tema<-theme
        tema1 <- modifyList(tema, list(layout.width = list(panel=1,
                                       ylab = 2, axis.left=1.0,
                                       left.padding=1, ylab.axis.padding=1,
                                       axis.panel=1)))
        tema2 <- modifyList(tema, list(layout.width = list(panel=1,
                                       ylab = 2, axis.left=1.0, left.padding=1,
                                       ylab.axis.padding=1, axis.panel=1)))
        temaT <- modifyList(tema, list(layout.heights = list(panel = c(1, 1))))
        p1 <- xyplot(Q~Pdc, groups=H, data=resumen,
                     ylab="Qd (m\u00b3/d)",type=c('l','g'),
                     par.settings = tema1)

        p1lab<-p1+glayer(panel.text(x[1], y[1], group.value, pos=2, cex=0.7))

        ##Pinto la regresion lineal porque Pnom~Pdc depende de la altura
        p2 <- xyplot(Pnom~Pdc, groups=H, data=resumen,
                     ylab="Pg",type=c('l','g'), #type=c('smooth','g'),
                     par.settings = tema2)
        p2lab<-p2+glayer(panel.text(x[1], y[1], group.value, pos=2, cex=0.7))

        p<-update(c(p1lab, p2lab, x.same = TRUE),
                  main=paste(title, '\nSP', pump$Qn, 'A', pump$stages, ' ',
                  'Gd ', Gd/1000," kWh/m\u00b2",sep=''),
                  layout = c(1, 2),
                  scales=list(x=list(draw=FALSE)),
                  xlab='',              #no dibuja el eje X
                  ylab = list(c("Qd (m\u00b3/d)","Pg (Wp)"), y = c(1/4, 3/4)),
                  par.settings = temaT
                  )
        print(p)
        result<-list(I=res,D=resumen, plot=p, param=param)
    } else {
        warning('lattice, latticeExtra packages are not all available. Thus, the plot could not be created')
        result<-list(I=res,D=resumen, param=param)
    }
}

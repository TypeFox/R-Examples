fPump <- function(pump, H){

  w1=3000 ##frecuencia rpm sincronica
  wm=2870 #frecuencia rpm con deslizamiento al aplicar tensión a 50 Hz
  s=(w1-wm)/w1
  fen=50 ## Frecuencia electrica nominal
  fmin=sqrt(H/pump$a)
  fmax=with(pump, (-b*Qmax+sqrt(b^2*Qmax^2-4*a*(c*Qmax^2-H)))/(2*a))
  ## fb es frecuencia de giro (Hz) de la bomba, 
  ##fe es la frecuencia electrica aplicada al motor
  ##que le hace girar a una frecuencia fb (y por tanto tambien a la bomba).
  fb=seq(fmin,min(60,fmax),length=1000) # La frecuencia maxima es 60
  fe=fb/(1-s)

###Caudal
  Q=with(pump, (-b*fb-sqrt(b^2*fb^2-4*c*(a*fb^2-H)))/(2*c))
  Qmin=0.1*pump$Qn*fb/50
  Q=Q+(Qmin-Q)*(Q<Qmin)
###Potencia hidraúlica
  Ph=2.725*Q*H

###Potencia mecánica
  Q50=50*Q/fb
  H50=H*(50/fb)^2
  etab=with(pump, j*Q50^2+k*Q50+l)
  Pb50=2.725*H50*Q50/etab
  Pb=Pb50*(fb/50)^3

###Potencia eléctrica
  Pbc=Pb*50/fe
  etam=with(pump, g*(Pbc/Pmn)^2+h*(Pbc/Pmn)+i)
  Pmc=Pbc/etam
  Pm=Pmc*fe/50
  Pac=Pm
  ##Pdc=Pm/(etac*(1-cab))

###Construyo funciones para caudal, la frecuencia y las potencias
###para ajustar la potencia de alterna.
  fQ<-splinefun(Pac,Q)
  fFreq<-splinefun(Pac,fe)
  fPb<-splinefun(Pac,Pb)
  fPh<-splinefun(Pac,Ph)
  lim=c(min(Pac),max(Pac))
  ##lim marca el rango de funcionamiento de la bomba
  result<-list(lim=lim,fQ=fQ,fPb=fPb,fPh=fPh,fFreq=fFreq)
}





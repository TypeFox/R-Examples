
############  to illustrate the graphical representation of
############   an earthquake focal mechanism (for teaching, for example)
readline("To Continue Hit Enter Key\n")

print("Demonstration of vectors for Focal Mechanisms")

TEACHFOC(65, 32, -34, up=TRUE)

readline("To Continue Hit Enter Key\n")


#############  animation of faults for teaching:sfanc
if(FALSE)
  {

    ##############  these are animations...to be run outside of the demo
anim= seq(from=0, to=1, by=.1) 

 strikeslip.fault(anim=anim, Light=c(45,90) )
readline("To Continue Hit Enter Key\n")

   normal.fault(45, anim=anim, KAPPA=4, Light=c(-20, 80))
readline("To Continue Hit Enter Key\n")

thrust.fault(anim=anim, KAPPA=4, Light=c(-20, 80))
readline("To Continue Hit Enter Key\n")
}

###############  end animation part

print("Illustration of differnt kinds of faults")

###############   PLOT ALL three styles
   par(mfrow=c(3,1))

    anim=.5

############
###   graphics.off();  X11(w=8, h=10)

    strikeslip.fault(anim=anim, Light=c(45,90) )

    MFOC1 = SDRfoc(65,90,1, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol1 = foc.color(foc.icolor(MFOC1$rake1), pal=1)

    justfocXY( MFOC1, fcol = Fcol1, .5, .7 , size = c(.4,.4) )
    title("Strike-slip fault")

   
############
    normal.fault(45, anim=anim, KAPPA=4, Light=c(-20, 80))

    MFOC2 = SDRfoc(135,45,-90, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol2 = foc.color(foc.icolor(MFOC2$rake1), pal=1)

    justfocXY( MFOC2, fcol = Fcol2, .5, 1 , size = c(.45,.45) )
    title("Normal fault")

############

    thrust.fault(anim=anim, KAPPA=4, Light=c(-20, 80))

    MFOC3 = SDRfoc(135,45,90, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol3 = foc.color(foc.icolor(MFOC3$rake1), pal=1)

    justfocXY( MFOC3, fcol = Fcol3, 0, -1 , size = c(.45,.45) )

    title("Reverse (Thrust) fault")
readline("To Continue Hit Enter Key\n")

###################################
############################# ############################# ############################# 
############################# ############################# ############################# 
############################# ############################# ############################# 


######  load("/home/lees/Progs/R_PAX/RFOC/data/KAMCORN.RData")

print("Load data from Kamchatka-Aleutian corner")

data(KAMCORN)
  par(mfrow=c(1,1))


print("Plot Focal Mechanisms Geographically")

plot(KAMCORN$LON, KAMCORN$LAT, xlab="LON", ylab="LAT" , main="Kamchatka-Aleutian Inersection", asp=1)

####  set up these vectors for later use
Paz =vector()
Pdip =vector()
Taz =vector()
Tdip =vector()
h = vector()
v = vector()

IFcol = vector()
Fcol = vector()

for(i in 1:length(KAMCORN$LON))
  {
    Msdr = CONVERTSDR(KAMCORN$STRIKE[i], KAMCORN$DIP[i], KAMCORN$RAKE[i]   )
  MEC = MRake(Msdr$M)
  MEC$UP = FALSE 
  IFcol[i] = foc.icolor(MEC$rake1)
    Fcol[i] = foc.color(IFcol[i], 1)
    
      az1 = Msdr$M$az1
  dip1 = Msdr$M$d1
  az2 = Msdr$M$az2
  dip2 = Msdr$M$d2
  BBB = Bfocvec(az1, dip1,  az2,  dip2)
  V = ternfoc.point(BBB$Bdip, Msdr$M$pd, Msdr$M$td )
 Paz[i] = Msdr$M$paz
  Pdip[i] = Msdr$M$pd
  Taz[i] = Msdr$M$taz
  Tdip[i] = Msdr$M$td
  h[i] = V$h
  v[i] = V$v

     justfocXY( MEC, fcol = Fcol[i], KAMCORN$LON[i], KAMCORN$LAT[i] , size = c(.21,.21) )
  }


readline("To Continue Hit Enter Key\n")

print("Plot Focal Mechanisms showing only fault plane and slip vector point")

#############  plot only fault planes with slip vector:


plot(KAMCORN$LON, KAMCORN$LAT, xlab="LON", ylab="LAT" , main="Kamchatka-Aleutian Inersection", asp=1, type='n')

 for(i in 1:length(KAMCORN$LON) )
    {
      Msdr = CONVERTSDR(KAMCORN$STRIKE[i], KAMCORN$DIP[i], KAMCORN$RAKE[i]   )
      MEC = MRake(Msdr$M)
      MEC$UP = FALSE 
      ##  points(fxy$x, fxy$y)
      nipXY( MEC,  KAMCORN$LON[i], KAMCORN$LAT[i] , fcol = Fcol[i], nipcol=Fcol[i],  size = c(.21,.21) )

    }

readline("To Continue Hit Enter Key\n")

 print("show summary steronet plot of P and T-axes from Kamchatka Data")

####################################################
####################################################
####################################################  combined P-T axes
####################################################


  net()
  PZZ =  focpoint(Paz, Pdip, col='red',  pch=3, lab="", UP=FALSE)
  TZZ =  focpoint(Taz, Tdip, col='blue',  pch=4, lab="", UP=FALSE)
  text(0,1.04,labels="P&T-axes Projected", font=2, cex=1.2)
  legend("topright",c("P", "T") , col=c('red','blue') , pch=c(3,4))


readline("To Continue Hit Enter Key\n")

  net()
  PZZ =  focpoint(Paz, Pdip, col='red',  pch=3, lab="", UP=FALSE)

 ALPH = alpha95(Paz, Pdip)

addsmallcirc(ALPH$Dr, ALPH$Ir, ALPH$Alph95, BALL.radius = 1, N = 25,
add = TRUE, lwd=1, col='blue')
readline("To Continue Hit Enter Key\n")

 net()
TZZ =  focpoint(Taz, Tdip, col='blue',  pch=4, lab="", UP=FALSE)
 ALPH = alpha95(Taz, Tdip)

addsmallcirc(ALPH$Dr, ALPH$Ir, ALPH$Alph95, BALL.radius = 1, N = 25,
add = TRUE, lwd=1, col='red')




readline("To Continue Hit Enter Key\n")

 print("Contour P and T-axis poles")

library(MASS)

 KP = kde2d(PZZ$x, PZZ$y, n=50, lims=c(-1, 1, -1, 1))
  KT = kde2d(TZZ$x, TZZ$y, n=50, lims=c(-1, 1, -1, 1) )



 par(mfrow=c(1,3))
   #################  plot 1
    par(mai=c(0.2,0,.2,0))
    CC = PLTcirc(PLOT=FALSE, add=FALSE,  ndiv=36,  angs=c(-pi, pi))

    image(KP$x, KP$y, KP$z, add=TRUE, col=terrain.colors(100))

    antipolygon(CC$x,CC$y,col="white")

    net(add=1)
    focpoint(Paz, Pdip, col='red',  pch=3, lab="", UP=FALSE)
    text(0,1.04,labels="P-axes 2D Density", font=2, cex=1.2)

   #################  plot 2
    CC = PLTcirc(PLOT=FALSE, add=FALSE,  ndiv=36,  angs=c(-pi, pi))


    image(KT$x, KT$y, KT$z, add=TRUE, col=terrain.colors(100))
    
    antipolygon(CC$x,CC$y,col="white")
    net(add=1)
    focpoint(Taz, Tdip, col='blue',  pch=4, lab="", UP=FALSE)
    text(0,1.04,labels="T-axes 2D Density", font=2, cex=1.2)


 #################  plot 3
  CC = PLTcirc(PLOT=FALSE, add=FALSE,  ndiv=36,  angs=c(-pi, pi))

    image(KP$x, KP$y, KP$z, add=TRUE, col=terrain.colors(100))


    ##  focpoint(Paz, Pdip, col='red',  pch=3, lab="", UP=FALSE)

    net(add=1)


    contour(KT$x, KT$y, KT$z, add=TRUE, lwd=1.2)
    

    antipolygon(CC$x,CC$y,col="white")
    text(0,1.04,labels="Combined P-T 2D Density", font=2, cex=1.2)
  

readline("To Continue Hit Enter Key\n")

 print("Show summary ternary plots of P and T-axes from Kamchatka Data")

####################################################
####################################################
####################################################  Ternary Plots
####################################################
graphics.off()
par(mfrow=c(1,1))

  PlotTernfoc(h,v,x=0, y=0, siz=1, fcols=Fcol, add=FALSE, LAB=TRUE)

MFOC1 = SDRfoc(65,90,1, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol1 = foc.color(foc.icolor(MFOC1$rake1), pal=1)
 MFOC2 = SDRfoc(135,45,-90, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol2 = foc.color(foc.icolor(MFOC2$rake1), pal=1)
 MFOC3 = SDRfoc(135,45,90, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)
    Fcol3 = foc.color(foc.icolor(MFOC3$rake1), pal=1)

justfocXY( MFOC3, fcol = Fcol3, 1.2, -0.9, size = c(.1,.1) )
justfocXY( MFOC2, fcol = Fcol2, -1.2, -0.9, size = c(.1,.1) )
justfocXY( MFOC1, fcol = Fcol1, 0, 1.414443+.2, size = c(.1,.1) )

mtext("Ternary Plot of focal mecahnisms", side=1, line = 1, font=2)

readline("To Continue Hit Enter Key\n")

######################################
graphics.off()
plot(KAMCORN$LON, KAMCORN$LAT, xlab="LON", ylab="LAT" , main="Kamchatka-Aleutian Inersection", asp=1, type='n')

 u = par("usr")

    RI = RectDense(KAMCORN$LON , KAMCORN$LAT, icut=20, u=u, ndivs=6)

    
    points(KAMCORN$LON , KAMCORN$LAT, pch='.', cex=2)

    rect(RI$icorns[,1],RI$icorns[,2],RI$icorns[,3],RI$icorns[,4], col=NA, border='blue')


######################################


for(i in 1:length(RI$ipass))
      {
        flag = KAMCORN$LON>RI$icorns[i,1]& KAMCORN$LAT>RI$icorns[i,2] & KAMCORN$LON<RI$icorns[i,3] & KAMCORN$LAT<RI$icorns[i,4]
        jh =h[flag]
        jv= v[flag]
        PlotTernfoc(jh,jv,x=mean(RI$icorns[i,c(1,3)]), y=mean(RI$icorns[i,c(2,4)]), siz=.8, fcols=Fcol[flag], add=TRUE)
      }
    legend("topright",focleg(1:7), col=foc.color(1:7, pal=1), pch=16)


readline("To Continue Hit Enter Key\n")

####################################################
####################################################
####################################################  Ternary Plots
print("show geographics distribution of P-T axes around corner")



plot(KAMCORN$LON, KAMCORN$LAT, xlab="LON", ylab="LAT" , main="Kamchatka-Aleutian Inersection", asp=1, type='n')
points(KAMCORN$LON, KAMCORN$LAT, pch='.', cex=2)


KPspat =  matrix(NA, nrow=length(RI$ipass), ncol=10)
KTspat =  matrix(NA, nrow=length(RI$ipass), ncol=10)
colnames(KTspat) = c("x", "y", "n", "Ir",  "Dr", "R", "K", "S", "Alph95" , "Kappa")
colnames(KPspat) = c("x", "y", "n", "Ir",  "Dr", "R", "K", "S", "Alph95" , "Kappa")


 for(i in 1:length(RI$ipass))
      {
        flag = KAMCORN$LON>RI$icorns[i,1]& KAMCORN$LAT>RI$icorns[i,2] & KAMCORN$LON<RI$icorns[i,3] & KAMCORN$LAT<RI$icorns[i,4]
        paz=Paz[flag]
        pdip=Pdip[flag]
         taz=Taz[flag]
        tdip=Tdip[flag]
        x=mean(RI$icorns[i,c(1,3)])
        y=mean(RI$icorns[i,c(2,4)])
        siz=(RI$icorns[1,3]-RI$icorns[1,1])/2.5

        PlotPTsmooth(paz, pdip, x=x, y=y, siz=siz, border=NA, bcol='white' , LABS=FALSE, add=FALSE, IMAGE=TRUE, CONT=FALSE)
          PlotPTsmooth(taz, tdip, x=x, y=y, siz=siz, border=NA, bcol='white' , LABS=FALSE, add=TRUE, IMAGE=FALSE, CONT=TRUE)

######dev.set(2)
######pnet(MN, add=FALSE)
  
    ALPH = alpha95(paz, pdip)
        n = length( paz)
      KPspat[i,] =   c(x, y, n, ALPH$Ir,  ALPH$Dr, ALPH$R, ALPH$K, ALPH$S, ALPH$Alph95 , ALPH$Kappa)

        ALPH = alpha95(taz, tdip) 
      KTspat[i,] =    c(x, y, n, ALPH$Ir,  ALPH$Dr, ALPH$R, ALPH$K, ALPH$S, ALPH$Alph95 , ALPH$Kappa)

        

      }



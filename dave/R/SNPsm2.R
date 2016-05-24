SNPsm2 <-
function(trange=100,tsl=5.0,diff=0.001,r6=NULL) {
# Vegetation map of the initial state
# -----------------------------------
  imax<- 40
  jmax<- 30
  igmap<- rep(6,imax*jmax)
  igmap<- matrix(igmap,ncol=jmax)
  mmax<- max(igmap)
  igmap<- rep(6,imax*jmax)
  igmap<- matrix(igmap,ncol=jmax)

  mmax<- max(igmap)
  vegtypes<- c("1 Pinus","2 Carex","3 Festuca","4 Deschampsia","5 Trisetum","6 Aconitum")
  colors<- c("darkolivegreen4","lightgreen","gold","darkorange","red1","darkred")
# colors<- c("darkred","red1","orange","yellow","lightgreen","darkgreen")
  igmap[7,19]<-c(4)
  igmap[8,19]<-c(5)
  igmap[8,22:27]<- c(4,4,4,5,4,5)
  igmap[9,21:27]<- c(4,4,5,4,5,4,5)
  igmap[10,20:26]<- c(4,5,3,4,4,5,4)
  igmap[11,19:26]<- c(4,4,4,5,4,4,5,4)
  igmap[12,18:26]<- c(5,4,4,4,4,4,4,3,4)
  igmap[13,17:26]<- c(5,5,4,4,2,1,1,1,2,4)
  igmap[14,16:25]<- c(4,4,5,4,2,1,1,1,1,4)
  igmap[15,14:24]<- c(5,4,4,4,4,2,1,1,1,1,2)
  igmap[16,13:24]<- c(5,4,4,4,4,5,3,2,1,1,2,3)
  igmap[17,12:23]<- c(4,4,4,4,4,4,3,3,3,1,2,3)
  igmap[18,11:23]<- c(5,5,4,4,5,4,3,5,3,3,2,2,3)
  igmap[19,10:22]<- c(5,5,5,4,4,4,3,3,3,4,2,1,3)
  igmap[20,09:22]<- c(5,5,5,5,5,5,4,4,3,3,2,1,1,3)
  igmap[21,08:21]<- c(5,4,5,5,5,4,4,5,4,3,3,4,2,2)
  igmap[22,06:21]<- c(5,5,4,4,5,5,5,4,5,4,4,3,2,3,4,3)
  igmap[23,06:20]<- c(5,5,5,5,4,5,5,4,5,4,4,4,1,2,5)
  igmap[24,05:19]<- c(5,5,5,5,5,4,4,4,5,4,4,3,3,2,3) 
  igmap[25,05:19]<- c(5,5,5,5,4,4,5,5,5,4,4,4,2,1,3)
  igmap[26,07:10]<- c(5,5,4,5)
  igmap[26,13:19]<- c(5,4,4,4,3,4,4)
  igmap[27,12:19]<- c(5,4,4,4,4,3,4,1)
  igmap[28,12:19]<- c(5,4,4,4,4,4,5,4)
  igmap[29,13:19]<- c(5,4,4,4,3,3,4)
  igmap[30,13:19]<- c(5,5,4,4,4,3,4)
  igmap[31,13:19]<- c(5,5,4,4,4,4,4)
  igmap[32,15:18]<- c(4,4,4,4)
  igmap[33,15:18]<- c(5,4,4,4)
  igmap[34,16:18]<- c(5,4,4)
  igmap[35,17:18]<- c(5,4)
  igmap<- 7-igmap                             # reverse code of vegetation types
  frame<- igmap   ; frame[frame!=1]<-0        # frame delimits Stabelchod
#
# THE model SNP
# -------------
#   par(mfrow=c(4,6),mar=c(0,0,0,0),omi=c(0,0,0,0),lwd=0.4)
   testrange<- is.null(trange)
   if(testrange == TRUE) trange<- 400                 # Time range
   testtsl<- is.null(tsl)
   if(testtsl == TRUE)   tsl<- 1.0                    # Time step length
   nt<- (trange/tsl)                                  # No of time steps
   t<- seq(0,trange-tsl,tsl)                          # Time vector
   tmax<- length(t) 
   imax<- 40
   jmax<- 30
   vmax<- 6 
   t.map<- 0                            # current output map
   testdiff<- is.null(diff)
   if(testdiff == TRUE) diff<- 0.002      # Spatial diffusion coefficient
   x<- rep(0,imax*jmax*tmax*vmax)         # 4-dimensional state vector
   x<- array(x,c(imax,jmax,tmax,vmax))
   gandl<- rep(0,imax*jmax*vmax)          # gains and losses from spatial process
   gandl<- array(gandl,c(imax,jmax,vmax))
   vegdef<- rep(0,36)
   vegdef<- array(vegdef,c(6,6))          # difining composition of vegetation types
   vegdef[1,]<- c(50,  10, 7,  2,  1,  0)
   vegdef[2,]<- c(17.5,35, 15, 3,  1,  0)
   vegdef[3,]<- c(17.5,35, 35, 30, 10, 1)
   vegdef[4,]<- c(10,  15, 35, 42, 15, 1)
   vegdef[5,]<- c(5,   5,  6,  20, 65, 8)
   vegdef[6,]<- c(0,   0,  2,  3,  8,  90)
   colors<- c(gray(0.10),gray(0.25),gray(0.40),gray(0.55),gray(0.70),gray(0.85))
   tmap<-igmap                                 # tmap is used for plotting dominant type
#  r<- c(0.030,0.035,0.045,0.045,0.045,0.050)  # 6 growth rates, model 415
#  r<- c(0.022,0.026,0.020,0.022,0.022,0.018)  # 6 growth rates, model 585
   testr6<- is.null(r6)                        # default growth rates
   if(testr6 == TRUE) r<- c(0.011,0.013,0.010,0.011,0.011,0.009)
   if(testr6 == FALSE) r<- r6
   K<- 100                                     # Carrying capacity
   dx<- rep(0,6)                                     
   colnames(vegdef)<- c("Aconitum","Trisetum","Deschampsia","Festuca","Carex","Pinus")
#  Initial 6-state in the first time layer
   for (i in 1:imax) for (j in 1:jmax) x[i,j,1,]<- vegdef[1:6,igmap[i,j]] 
#  Plot discrete map of state at t=1, using correlation
   for (i in 1:imax) for (j in 1:jmax) {
       tmap[i,j]<- which.max(cor(vegdef,x[i,j,1,1:6]))
   }
#
# Simulation starts here
#
  cat("Simulation in progress . . .\n")
  for(t in 2:nt) {                                        # Main time loop
      for(i in 1:imax) for(j in 1:jmax) {                 # Processing spatial grid, time model
# Differential equations part 1: Logistic growth of vegetation types
       dx[1]<- r[1]*x[i,j,t-1,1]*((K- x[i,j,t-1,1])/K)              
       dx[2]<- r[2]*x[i,j,t-1,2]*((K- x[i,j,t-1,1]-x[i,j,t-1,2])/K)
       dx[3]<- r[3]*x[i,j,t-1,3]*((K- x[i,j,t-1,1]-x[i,j,t-1,2]-x[i,j,t-1,3])/K)
       dx[4]<- r[4]*x[i,j,t-1,4]*((K- x[i,j,t-1,1]-x[i,j,t-1,2]-x[i,j,t-1,3]-x[i,j,t-1,4])/K)
       dx[5]<- r[5]*x[i,j,t-1,5]*((K- x[i,j,t-1,1]-x[i,j,t-1,2]-x[i,j,t-1,3]-x[i,j,t-1,4]-x[i,j,t-1,5])/K)
       dx[6]<- r[6]*x[i,j,t-1,6]*((K- x[i,j,t-1,1]-x[i,j,t-1,2]-x[i,j,t-1,3]-x[i,j,t-1,4]-x[i,j,t-1,5]-x[i,j,t-1,6])/K)
# Differential equations part 2: Keeping gains and losses balanced
       if(dx[2]!=0) dx[2]<- dx[2] - (dx[1])*(x[i,j,t-1,2]/sum(x[i,j,t-1,2:6]))
       if(dx[3]!=0) dx[3]<- dx[3] - (dx[1]+dx[2])*(x[i,j,t-1,3]/sum(x[i,j,t-1,3:6]))
       if(dx[4]!=0) dx[4]<- dx[4] - (dx[1]+dx[2]+dx[3])*(x[i,j,t-1,4]/sum(x[i,j,t-1,4:6]))
       if(dx[5]!=0) dx[5]<- dx[5] - (dx[1]+dx[2]+dx[3]+dx[4])*(x[i,j,t-1,5]/sum(x[i,j,t-1,5:6]))
       if(dx[6]!=0) dx[6]<- dx[6] - (dx[1]+dx[2]+dx[3]+dx[4]+dx[5])*(x[i,j,t-1,6]/sum(x[i,j,t-1,6]))
# Numerical integration
       for(v in 1:6){
           x[i,j,t,v]<- x[i,j,t-1,v]+(dx[v]*tsl)
       }
      }
      for(i in 1:imax) for(j in 1:jmax) {                 # Processing spatial grid, space model
           gandl[i,j,1:6]<- x[max(1,i-1),j,t,1:6]+x[min(imax,i+1),j,t,1:6]
           gandl[i,j,1:6]<- gandl[i,j,1:6]+x[i,max(1,j-1),t,1:6]+x[i,min(jmax,j+1),t,1:6]
           gandl[i,j,1:6]<- gandl[i,j,1:6]-(4*x[i,j,t,1:6])
      }
      x[,,t,1:6]<- x[,,t,1:6]+(diff*tsl*gandl[,,1:6])

  }
  o.SNPsm<- list(n.time.steps=nt,imax=imax,jmax=jmax,time.step.length=tsl,veg.types=vegtypes,vegdef=vegdef,growth.rates=r,sim.data=x,tmap=tmap,igmap=igmap,frame=frame)
}

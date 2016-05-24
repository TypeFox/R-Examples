SNPtm2 <-
function(trange=100,tsl=1.0,x6=NULL,r6=NULL) { 
   testrange<- is.null(trange)
   if(testrange == TRUE) trange<- 400                 # Time range
   testtsl<- is.null(tsl)
   if(testtsl == TRUE)   tsl<- 1.0                    # Time step length
   x<- rep(0,trange*6/tsl)
   x<- matrix(x,ncol=6)
   testr<- is.null(r6)
   if(testr == TRUE) r<- c(0.025,0.040,0.045,0.045,0.045,0.045)  # 6 growth rates, model 410
   if(testr == FALSE) r<- r6
   testx<- is.null(x6)
   if(testx == TRUE) x[1,]<- c(0.50,0.70,1.70,6.0,13.0,78.0)      # 6 initial conditions, model 410
   if(testx == FALSE) x[1,]<- x6
   nt<- (trange/tsl)               # No of time steps
   tv<- seq(0,trange-tsl,tsl)      # Time vector
   symbols<- c(1,4,8,17,22,18)
   colors<- c("darkolivegreen4","lightgreen","gold","darkorange","red1","darkred")
   vegtypes<- c("Aconitum","Trisetum","Deschampsia","Festuca","Carex","Pinus")
   K<- 100                                     # Carrying capacity
   dx<- rep(0,6)                                     
 
   for(t in 2:nt) {
# Differential equations part 1: Logistic growth of vegetation types
       dx[1]<- r[1]*x[t-1,1]*((K- x[t-1,1])/K)              
       dx[2]<- r[2]*x[t-1,2]*((K- x[t-1,1]-x[t-1,2])/K)
       dx[3]<- r[3]*x[t-1,3]*((K- x[t-1,1]-x[t-1,2]-x[t-1,3])/K)
       dx[4]<- r[4]*x[t-1,4]*((K- x[t-1,1]-x[t-1,2]-x[t-1,3]-x[t-1,4])/K)
       dx[5]<- r[5]*x[t-1,5]*((K- x[t-1,1]-x[t-1,2]-x[t-1,3]-x[t-1,4]-x[t-1,5])/K)
       dx[6]<- r[6]*x[t-1,6]*((K- x[t-1,1]-x[t-1,2]-x[t-1,3]-x[t-1,4]-x[t-1,5]-x[t-1,6])/K)
# Differential equations part 2: Keeping gains and losses balanced
       dx[2]<- dx[2] - (dx[1])*(x[t-1,2]/sum(x[t-1,2:6]))
       dx[3]<- dx[3] - (dx[1]+dx[2])*(x[t-1,3]/sum(x[t-1,3:6]))
       dx[4]<- dx[4] - (dx[1]+dx[2]+dx[3])*(x[t-1,4]/sum(x[t-1,4:6]))
       dx[5]<- dx[5] - (dx[1]+dx[2]+dx[3]+dx[4])*(x[t-1,5]/sum(x[t-1,5:6]))
       dx[6]<- dx[6] - (dx[1]+dx[2]+dx[3]+dx[4]+dx[5])*(x[t-1,6]/sum(x[t-1,6]))
# Numerical integration
       for(v in 1:6){
           x[t,v]<- x[t-1,v]+(dx[v]*tsl)
       }
   }
   o.SNP<- list(n.time.steps=trange,time.step.length=tsl,time.vector=tv,veg.types=vegtypes,growth.rates=r,initial.cond=x[1,],sim.data=x)
}

ros <-
function(modeltype,w,s,delta,mx.dead,h,m,u,slope) {
  
  if (is.data.frame(m)==F) {
    m<-data.frame(t(m))}
  
  m <- m/100
  
  if (nrow(w)==1 || is.data.frame(w)==F) {
  w<-data.frame(t(matrix(rep(as.numeric(w[1:5]),nrow(m)),5,nrow(m))))}
  
  w <-w/10*0.2048
  
  if (nrow(s)==1 || is.data.frame(s)==F) {
  s <- data.frame(t(matrix(rep(as.numeric(s[1:5]),nrow(m)),5,nrow(m))))}
  
  s<-s/3.281
  
  if (nrow(h)==1 || is.data.frame(h)==F) {
    h <- data.frame(t(matrix(rep(as.numeric(h[1:5]),nrow(m)),5,nrow(m))))}
  
  h<-h*0.429922614
  
  if (is.data.frame(delta)==F) {
    delta <- data.frame(t(matrix(rep(delta,nrow(m)),1,nrow(m))))}
  
  delta <-(delta)*0.0328084
  
  if (is.data.frame(mx.dead)==F) {
    mx.dead <- data.frame(t(matrix(rep(mx.dead,nrow(m)),1,nrow(m))))}
  
  mx.dead <-(mx.dead)/100
  
  u<-(u)*54.6806649
  
  slope <-(slope)/100
  
 
  # set lower limit of .3 for moisture(herb)
  
  error.mh<-c()
  for (i in 1:nrow(m)) {
    if (length(which(m[i,4]>0 & m[i,4]<0.3))>0) 
    {
    if (m[i,4]>0 & m[i,4]<.3) error.mh[i]<-1
    }
  }
  
  if (length(error.mh)>0) {
  if (sum(na.omit(error.mh))>0) stop("Moisture of herbaceous fuels should be >= 30%")}
  
  # transfer partially cured herbaceous fuels to dead
    
  kt=rep(NA,nrow(m))
  if (modeltype=="D") {
    kt[which(m[,4]>=.3 & m[,4]<1.2)]<-(1.20-m[which(m[,4]>=.3 & m[,4]<1.2),4])/.9
    kt[which(is.na(kt))]<-0
    
    # weighting SAV from transferred cured herbaceous
    f1<-w[,1]*s[,1]/32
    f4<-(w[,4]*kt)*s[,4]/32
    s[,1]<-(f1*s[,1] +f4*s[,4])/(f1+f4)
    
    #       A<-f[,1]+f[,4]
    #       h[,1]<-(f[,1]/A)*h[,1]+(f[,4]/A)*h[,4]
    #       apparently not implemented in BehavePlus
    
    w[,1]=w[,1]+w[,4]*kt
    w[,4]=w[,4]-w[,4]*kt
    
  }
  
  rhop= 32 # 513*0.0624279606 Scott and Burgan (2005)
  st=0.0555
  se=0.01
  
  #area fractions and weights
  a<-cbind(s[,1]*w[,1],s[,2]*w[,2],s[,3]*w[,3],s[,4]*w[,4],s[,5]*w[,5])/rhop
  a.dead=a[,1]+a[,2]+a[,3]
  a.live=a[,4]+a[,5]
  a.tot=(a.dead+a.live)
  f<-cbind(a[,1]/a.dead,a[,2]/a.dead,a[,3]/a.dead,a[,4]/a.live,a[,5]/a.live)
  f[which(a.live==0),4]=0
  f[which(a.live==0),5]=0
  f[which(a.dead==0),1]=0
  f[which(a.dead==0),2]=0
  f[which(a.dead==0),3]=0
  
  f.dead=a.dead/a.tot
  f.live=a.live/a.tot
  
  #net (weighted) fuel loadings
  wn=w*(1-st) # Albini 1976
  wn.dead=f[,1]*wn[,1]+f[,2]*wn[,2]+f[,3]*wn[,3]
  #wn.live=sum(f[4:5]*wn[4:5])
  wn.live=wn[,4]+wn[,5] # corrects models w/ 2 live fuel classes  (undocumented)
  
  #weighted fuel moisture
  mf.dead=f[,1]*m[,1]+f[,2]*m[,2]+f[,3]*m[,3] 
  mf.live=f[,4]*m[,4]+f[,5]*m[,5]
  
  #weighted SAV ratio
  sigma.dead=f[,1]*s[,1]+f[,2]*s[,2]+f[,3]*s[,3] 
  sigma.live=f[,4]*s[,4]+f[,5]*s[,5]  
  sigma.tot=(f.dead*sigma.dead+f.live*sigma.live) #characteristic SAV
  
  #weighted heat content
  h.dead=f[,1]*h[,1]+f[,2]*h[,2]+f[,3]*h[,3] 
  h.live=f[,4]*h[,4]+f[,5]*h[,5]
  
  #mean packing ratio for fuel complex
  beta=1/delta*(w[,1]/rhop+w[,2]/rhop+w[,3]/rhop+w[,4]/rhop+w[,5]/rhop)
  
  #live fuel moisture of extinction
  
  W=(w[,1]*exp(-138/s[,1])+w[,2]*exp(-138/s[,2])+w[,3]*exp(-138/s[,3]))/
    (w[,4]*exp(-500/s[,4])+w[,5]*exp(-500/s[,5]))
  
  mfpd=(w[,1]*m[,1]*exp(-138/s[,1])+w[,2]*m[,2]*exp(-138/s[,2])+w[,3]*m[,3]*exp(-138/s[,3]))/(w[,1]*exp(-138/s[,1])+w[,2]*exp(-138/s[,2])+w[,3]*exp(-138/s[,3]))
  
  mx.live=2.9*W*(1-mfpd/mx.dead)-0.226
  
  if (length(which(mx.live[,1]<mx.dead[,1]))>0)
  {mx.live[which(mx.live[,1]<mx.dead[,1]),]<-mx.dead[which(mx.live[,1]<mx.dead[,1]),] }
  
  if (length(which(W==Inf))>0)
  {mx.live[which(W==Inf),1]<-mx.dead[which(W==Inf),1] }
  
  #damping coefficients
  ns=0.174*se[1]^(-0.19)
  
  nm.dead=1-2.59*(mf.dead/mx.dead)+5.11*(mf.dead/mx.dead)^2-3.52*(mf.dead/mx.dead)^3
  nm.live=1-2.59*(mf.live/mx.live)+5.11*(mf.live/mx.live)^2-3.52*(mf.live/mx.live)^3
  
  if (length(which(mf.dead>mx.dead[,1]))>0)
  {nm.dead[which(mf.dead>mx.dead[,1]),]=0} 
  
  if (length(which(mf.live>mx.live[,1]))>0)
  {nm.live[which(mf.live>mx.live[,1]),]=0} 
  # Andrews 2013 pag.E
  
  #optimum packing ratio
  beta.op=3.348*sigma.tot^(-0.8189)
  rpr=beta/beta.op #relative packing ratio
  
  #maximum reaction velocity
  gamma.max=(sigma.tot^1.5)/(495+0.0594*sigma.tot^(1.5))
  
  #reaction intensity
  sum.dead=(wn.dead*h.dead*nm.dead*ns)
  sum.live=(wn.live*h.live*nm.live*ns)
  
  #A=(6.7229*sigma.tot^0.1-7.27)
  A=133*sigma.tot^(-0.7913) #alternate formulation from Albini (1976)
  ir.dead=gamma.max*(rpr*exp(1-rpr))^A*sum.dead #*f.dead removed by Frandsen 73
  ir.live=gamma.max*(rpr*exp(1-rpr))^A*sum.live #*f.live removed by Frandsen 73
  ir=ir.dead+ir.live
  
  #propagating flux ratio
  xi=(192+0.2595*sigma.tot)^(-1)*exp((0.792+0.681*sigma.tot^0.5)*(beta+0.1))
  
  #wind coefficient
  C=7.47*exp(-0.133*sigma.tot^.55)
  B=0.02526*sigma.tot^.54
  E=0.715*exp(-3.59*10^(-4)*sigma.tot)
  fw=C*u^B*rpr^(-E)
  
  #slope coefficient
  fs=5.275*beta^(-0.3)*slope^2
  
  #heat sink (denominatore ROS)
  rhob=(1/delta)*(w[,1]+w[,2]+w[,3]+w[,4]+w[,5]) #ovendry bulk density
  qig=250+1116*(m)
  
  if (length(which(w[,1]==0))>0)  
{qig[which(w[,1]==0),1]<-0}
  
  if (length(which(w[,2]==0))>0)    
{qig[which(w[,2]==0),2]<-0}
  
  if (length(which(w[,3]==0))>0)  
{qig[which(w[,3]==0),3]<-0}
  
  if (length(which(w[,4]==0))>0)    
{qig[which(w[,4]==0),4]<-0}
  
  if (length(which(w[,5]==0))>0)  
{qig[which(w[,5]==0),5]<-0}
  
  eps=f.dead*(f[,1]*qig[,1]*exp(-138/s[,1])+f[,2]*qig[,2]*exp(-138/s[,2])+f[,3]*qig[,3]*exp(-138/s[,3]))+f.live*(f[,4]*qig[,4]*exp(-138/s[,4])+f[,5]*qig[,5]*exp(-138/s[,5]))
  
  #ROS
  r = (ir*xi*(1+fw+fs))/(rhob*eps)
  r = (0.3048*r)
  
  # return values for rothermel function
  
  output=list(mf.dead*100,mf.live*100,(mx.live*100)[,1],sigma.tot*3.281,rhob[,1]*16.0184634,beta[,1],rpr[,1],ir.dead[,1]* 0.1893,ir.live[,1]* 0.1893,ir[,1]* 0.1893,as.data.frame(fw)[,1],as.data.frame(fs)[,1],(ir*xi*(1+fw+fs))[,1]* 0.1893,(rhob*eps)[,1]*37.25894580781,r[,1])
  
  names(output)=c("Characteristic dead fuel moisture [%]","Characteristic live fuel moisture [%]","Live fuel moisture of extinction [%]","Characteristic SAV [m2/m3]","Bulk density [kg/m3]","Packing ratio [dimensionless]","Relative packing ratio [dimensionless]","Dead fuel Reaction intensity [kW/m2]","Live fuel Reaction intensity [kW/m2]","Reaction intensity [kW/m2]","Wind factor [0-100]","Slope factor [0-1]","Heat source [kW/m2]","Heat sink [kJ/m3]","ROS [m/min]")
  
  return(c(lapply(output[1:5],FUN=round,digits=2),
         lapply(output[6],FUN=round,digits=4),
         lapply(output[7:15],FUN=round,digits=2)))

  
  #missing: upper threshold for u > 0.9 ir
  
  # BehavePlus includes the option of not imposing
  # the wind limit, based on a reanalysis of the McArthur (1969)
  # data from which the wind limit function was derived and on
  # recent data that do not support the limit (Andrews et al. 2013).
  
  # missing: equations for fire spread in all directions 
  # (da codice GRASS o Jiménez et al. 2008)
  
}

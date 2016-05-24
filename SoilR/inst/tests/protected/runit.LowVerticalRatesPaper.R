#
# vim:set ff=unix expandtab ts=2 sw=2:
# This test checks the code published with Sierra et al. (2013, Biogeosciences 10: 3455)
test.LowVerticalRatesPaper=function(){
     require(RUnit)

     #Create new dataset with pre- and post-bomb data
     C14d=data.frame(AD=c(1950-IntCal09[1:(dim(IntCal09)[1]-1),1],C14Atm_NH[52:111,1]),D14C=c(IntCal09[1:(dim(IntCal09)[1]-1),4],C14Atm_NH[52:111,2]))
     
     #Create model as an object. See manuscript for model description
     ZafireModel=function(t, ks,  C0, In, a31,  a42, g, xi=1, FcAtm, F0, lambda=-0.0001209681, lag=0, solver=deSolve.lsoda.wrapper, pass=FALSE){
          t_start=min(t)
          t_stop=max(t)
          if(length(ks)!=4) stop("ks must be of length = 4")
          if(length(C0)!=4) stop("the vector with initial conditions must be of length = 4")
          
          if(length(In)==1) inputFluxes=new("TimeMap",
                                            t_start,
                                            t_stop,
                                            function(t){matrix(nrow=4,ncol=1,c(In*g,In*(1-g),0,0))}
          )
          
          if(length(In)==4) inputFluxes=new("TimeMap",
                                            t_start,
                                            t_stop,
                                            function(t){matrix(nrow=4,ncol=1,c(In[1],In[2],In[3],In[4]))}
          )
          
          if(class(In)=="data.frame"){
               x=In[,1]  
               y=In[,2]  
               inputFlux=function(t0){as.numeric(spline(x,y,xout=t0)[2])}
               inputFluxes=new("TimeMap",
                               t_start,
                               t_stop,
                               function(t){matrix(nrow=4,ncol=1,c(inputFlux(t),0,0,0))}
               )   
          }
          
          if(class(In)=="TimeMap") inputFluxes=In
          
          if(length(xi)==1) fX=function(t){xi}
          if(class(xi)=="data.frame"){
               X=xi[,1]
               Y=xi[,2]
               fX=function(t){as.numeric(spline(X,Y,xout=t)[2])}
          }
          
          A=-abs(diag(ks))
          A[3,1]=a31
          A[4,2]=a42
          
          At=new(Class="DecompositionOperator",
                 t_start,
                 t_stop,
                 function(t){
                      fX(t)*A
                 }
          ) 
          
          Fc=FcAtm.from.Dataframe(FcAtm,lag,format="Delta14C")
          ivF=SoilR.F0.new(F0,"Delta14C")
          mod=GeneralModel_14(t=t,A=At,ivList=C0,initialValF=ivF,inputFluxes=inputFluxes,Fc=Fc,di=lambda)
     }     
     
     #### Run model
     
     #Parameters for the terra-firme forest (ZAR-04)
     years=seq(-1000,2010,by=0.5)
     LIZ4=4.72 #Litter inputs for ZAR-04
     RI20=1.53 # Root inputs for ZAR-04 from 0 to 20 cm
     Y20=1-(0.962^20) #Proportion of root inputs up to 20cm
     RT=RI20/Y20 #Total root inputs in profile
     Y55=1-(0.962^55) #Proportion of root inputs up to 55cm
     RI55= RT*Y55 #Root inputs from 0 to 55cm
     RI2055=RI55-RI20 #Root inputs from 20 to 55 cm
     gamma1=0.1; gamma2=0.98
     IZ4=c((LIZ4+RI20)*gamma1,(LIZ4+RI20)*(1-gamma1),RI2055*gamma2,RI2055*(1-gamma2)) #Litter inputs
     ksZ4=c(k1=2, k2=0.095, k3=1, k4=1/5000) #Decomposition rates
     F0=C14d[which(C14d[,1]==min(years)),2] #Initial value of atmospheric fraction
     
     #Run the model and calcuate radiocarbon in C and Rh 
     Z4=ZafireModel(t=years,ks=ksZ4, C0=c(20,20,20,47), In=IZ4, g=0.1, a31=0.3*ksZ4[1], a42=0.0001*ksZ4[2], FcAtm=C14d,F0=F0)
     RtZ4=getReleaseFlux(Z4)
     C14Z4=getC14(Z4)
     C14pZ4=getF14(Z4)
     CtZ4=Z4["C"]
     R14mZ4=getReleaseFlux14(Z4)
     R14topZ4=Delta14C_from_AbsoluteFractionModern(rowSums(R14mZ4[,1:2])/rowSums(RtZ4[,1:2]))
     C14topZ4=Delta14C_from_AbsoluteFractionModern(rowSums(C14Z4[,1:2])/rowSums(CtZ4[,1:2]))
     R14subZ4=Delta14C_from_AbsoluteFractionModern(rowSums(R14mZ4[,3:4])/rowSums(RtZ4[,3:4]))
     C14subZ4=Delta14C_from_AbsoluteFractionModern(rowSums(C14Z4[,3:4])/rowSums(CtZ4[,3:4]))
     
     #The vertical transfers are calcuated as
     VtransferZ4=RtZ4[,1:2]%*%c(0.3, 0.0001)
     
     
     #Parameters for the varillal (ZAR-01)
     IZ1=2.48+2.98 #Litter inputs for ZAR-01
     ksZ1=c(k1=2, k2=0.11, k3=1, k4=1/5000)
     
     Z1=ZafireModel(t=years,ks=ksZ1, C0=c(20,20,20,47), In=IZ1, g=0.1, a31=0.017*ksZ1[1], a42=0.0025*ksZ1[2], FcAtm=C14d,F0=F0)
     RtZ1=getReleaseFlux(Z1)
     C14Z1=getC14(Z1)
     C14pZ1=getF14(Z1)
     CtZ1=Z1["C"]
     R14mZ1=getReleaseFlux14(Z1)
     R14topZ1=Delta14C_from_AbsoluteFractionModern(rowSums(R14mZ1[,1:2])/rowSums(RtZ1[,1:2]))
     C14topZ1=Delta14C_from_AbsoluteFractionModern(rowSums(C14Z1[,1:2])/rowSums(CtZ1[,1:2]))
     R14subZ1=Delta14C_from_AbsoluteFractionModern(rowSums(R14mZ1[,3:4])/rowSums(RtZ1[,3:4]))
     C14subZ1=Delta14C_from_AbsoluteFractionModern(rowSums(C14Z1[,3:4])/rowSums(CtZ1[,3:4]))
     
     #The vertical transfers are calcuated as
     VtransferZ1=RtZ1[,1:2]%*%c(0.017, 0.0025)
     
     ######## Figure 4
     par(mfrow=c(2,1), mar=c(5, 5, 1, 2))
     plot(C14d,type="l",xlab="Year",ylab=expression(paste(Delta^14,"C (\u2030)")),xlim=c(1950,2010),ylim=c(-200,1000),lty=2) 
     lines(years,C14topZ4)
     lines(years, R14topZ4,col=4)
     lines(years,R14subZ4,col=2)
     lines(years,C14subZ4,col=3)
     #points(2010,Prof[1,5],col=4,pch=20)
     #arrows(2010,Prof[1,5]+Prof[1,6],2010,Prof[1,5]-Prof[1,6],col=4,angle=90,code=3,length=0.1)
     #points(2010,Prof[1,7],pch=20)
     #arrows(2010,Prof[1,7]+Prof[1,8],2010,Prof[1,7]-Prof[1,8],angle=90,code=3,length=0.1)
     #points(2010,Prof[3,5],col=2,pch=20)
     #arrows(2010,Prof[3,5]+Prof[3,6],2010,Prof[3,5]-Prof[3,6],angle=90,code=3,length=0.1,col=2)
     #points(2010,Prof[3,7],col=3,pch=20)
     #arrows(2010,Prof[3,7]+Prof[3,8],2010,Prof[3,7]-Prof[3,8],angle=90,code=3,length=0.1,col=3)
     legend("topright",c(expression(paste(Delta^14,"C atmosphere")), expression(paste(Delta^14,"C-",C0[2], " from topsoil")), 
                         expression(paste(Delta^14,"C in SOM from topsoil")), 
                         expression(paste(Delta^14,"C-",C0[2], " from subsoil")), expression(paste(Delta^14,"C in SOM from subsoil"))),
            lty=c(2,1,1,1,1),col=c(1,4,1,2,3),bty="n")
     legend("topleft",legend="a",bty="n")
     
     plot(C14d,type="l",xlab="Year",ylab=expression(paste(Delta^14,"C (\u2030)")),xlim=c(1950,2010),ylim=c(-200,1000),lty=2) 
     lines(years,C14topZ1)
     lines(years, R14topZ1,col=4)
     lines(years,R14subZ1,col=2)
     lines(years,C14subZ1,col=3)
#      points(2010,Prof[5,5],col=4,pch=20)
#      arrows(2010,Prof[5,5]+Prof[5,6],2010,Prof[5,5]-Prof[5,6],angle=90,code=3,length=0.1,col=4)
#      points(2010,Prof[5,7],col=4,pch=20)
#      arrows(2010,Prof[5,7]+Prof[5,8],2010,Prof[5,7]-Prof[5,8],angle=90,code=3,length=0.1)
#      points(2010,Prof[6,5],col=2,pch=20)
#      arrows(2010,Prof[6,5]+Prof[6,6],2010,Prof[6,5]-Prof[6,6],angle=90,code=3,length=0.1,col=2)
#      points(2010,Prof[6,7],col=3,pch=20)
#      arrows(2010,Prof[6,7]+Prof[6,8],2010,Prof[6,7]-Prof[6,8],angle=90,code=3,length=0.1,col=3)
     #legend("topright",c(expression(paste(Delta^14,"C atmosphere")), expression(paste(Delta^14,"C-",C0[2], " from topsoil")), 
     #                   expression(paste(Delta^14,"C in SOM from topsoil")), 
     #                  expression(paste(Delta^14,"C-",C0[2], " from subsoil")), expression(paste(Delta^14,"C in SOM from subsoil"))),
     #    lty=c(2,1,1,1,1),col=c(1,4,1,2,3),bty="n")
     legend("topleft",legend="b",bty="n")
     par(mfrow=c(1,1))
     
     
}     

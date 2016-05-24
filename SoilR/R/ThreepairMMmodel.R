#
# vim:set ff=unix expandtab ts=2 sw=2:
ThreepairMMmodel<- structure(
     function # Implementation of a 6-pool Michaelis-Menten model
     ### This function implements a 6-pool Michaelis-Meneten model with pairs of microbial biomass and substrate pools.
     (t, ##<< vector of times to calculate a solution.
      ks, ##<< a vector of length 3 representing SOM decomposition rate (m3 d-1 (gCB)-1)
      kb, ##<< a vector of length 3 representing microbial decay rate (d-1)
      Km, ##<< a vector of length 3 representing the Michaelis constant (g m-3)
      r, ##<< a vector of length 3 representing the respired carbon fraction (unitless)
      Af=1, ##<< a scalar representing the Activity factor; i.e. a temperature and moisture modifier (unitless)
      ADD, ##<< a vector of length 3 representing the annual C input to the soil (g m-3 d-1)
      ival # a vector of length 6 with the initial values of the SOM pools and the microbial biomass pools (g m-3)
     )
     {
          t_start=min(t)
          t_end=max(t)
          nr=6
          if(length(ival)!=6) stop("The vector of initial values ival must be of length 6")
          
          
          #function with system of equations
          f=function(C,t){
               S1=C[1] #SOM pool 1
               B1=C[2] #Microbial biomass pool 1
               S2=C[3] #SOM pool 2
               B2=C[4] #Microbial biomass pool 2
               S3=C[5] #SOM pool 3
               B3=C[6] #Microbial biomass pool3
               O=matrix(byrow=TRUE,nrow=6,ncol=1,c((Af*ks[1])*B1*(S1/(Km[1]+S1)),
                                                   kb[1]*B1,
                                                   (Af*ks[2])*B2*(S2/(Km[2]+S2)),
                                                   kb[2]*B2,
                                                   (Af*ks[3])*B3*(S3/(Km[3]+S3)),
                                                   kb[3]*B3))
               return(O)
          }
          
          #List of transfer coefficients for T matrix
          alpha=list()
          alpha[["1_to_2"]]=function(C,t){
               1-r[1]
          }
          alpha[["2_to_1"]]=function(C,t){
               1
          }
          
          alpha[["3_to_4"]]=function(C,t){
               1-r[2]
          }
          alpha[["4_to_3"]]=function(C,t){
               1
          }
          alpha[["5_to_6"]]=function(C,t){
               1-r[3]
          }
          alpha[["6_to_5"]]=function(C,t){
               1
          }

          Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
          
          inputrates=BoundInFlux(
            function(t){
              matrix(
                nrow=nr,
                ncol=1,
                c(ADD[1],0,ADD[2],0,ADD[3],0)
              )
            },
            t_start,
            t_end
          )
          
          modnl=GeneralNlModel( t, Anl, ival, inputrates, deSolve.lsoda.wrapper)
          
          return(modnl)
          ### An object of class NlModel that can be further queried.
     }
     ,
     ex=function(){
          
          days=seq(0,1000)
          
          #Run the model with default parameter values
          MMmodel=ThreepairMMmodel(t=days,ival=rep(c(100,10),3),ks=c(0.1,0.05,0.01),
                                   kb=c(0.005,0.001,0.0005),Km=c(100,150,200),r=c(0.9,0.9,0.9),
                                   ADD=c(3,1,0.5))
          Cpools=getC(MMmodel)
          
          #Time solution
          matplot(days,Cpools,type="l",ylab="Concentrations",xlab="Days",lty=rep(1:2,3),
                  ylim=c(0,max(Cpools)*1.2),col=rep(1:3,each=2),main="Multi-substrate microbial model")
          legend("topright",c("Substrate 1", "Microbial biomass 1", 
                             "Substrate 2", "Microbial biomass 2",
                             "Substrate 3", "Microbial biomass 3"),lty=rep(1:2,3),col=rep(1:3,each=2),
                 bty="n")
          
          
          #State-space diagram
          plot(Cpools[,2],Cpools[,1],type="l",ylab="Substrate",xlab="Microbial biomass")
          lines(Cpools[,4],Cpools[,3],col=2)
          lines(Cpools[,6],Cpools[,5],col=3)
          legend("topright",c("Substrate-Enzyme pair 1","Substrate-Enzyme pair 2",
                              "Substrate-Enzyme pair 3"),col=1:3,lty=1,bty="n")
          
          #Microbial biomass over time
          plot(days,Cpools[,2],type="l",col=2,xlab="Days",ylab="Microbial biomass")
          
          
     }
)

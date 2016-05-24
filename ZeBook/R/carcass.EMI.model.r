################################ FUNCTIONS #####################################
#' @title The Carcass (growth of beef cattle) model with energy as input
#' @description \strong{Model description.}
#' This model is proposed by Hoch et. al (2004) to represent the growth of cattle and the relative body composition of diferent types of animals depending on nuritionnal conditions.
#' It simulates the dynamics of changes in the composition of the body fat and proteins according to nutrient intake. The system is represented by four state variables: the protein and fat in the carcass (resp. ProtC and LIPC) and other tissues (resp. ProtNC and LipNC) grouped under the name of the fifth district (again, gastrointestinal tract, skin .. .). These variables depend on time, the time step used is dt = 1 day.
#' The model is defined by 20 equations, with a total of 18 parameters for the described process.
#' @param protcmax : amounts of protein in the carcass of the adult animal (kg)
#' @param protncmax : amounts of protein in the 5th district of the adult animal (kg)
#' @param alphac : maximum protein synthesis rate in the frame (excluding basal metabolism) (j-1)
#' @param alphanc : maximum rate of protein synthesis in the 5th district (except basal metabolism) (j-1)
#' @param gammac : Maximum rate of protein degradation in the frame (excluding basal metabolism) (j-1)
#' @param gammanc : maximum rate of protein degradation in the 5th district (except basal metabolism) (j-1)
#' @param lip0 : maximum lipid concentration to the theoretical physiological age (percent)
#' @param lipc1 : increase coefficient of the maximum lipid concentration with the physiological age of the carcase (percent)
#' @param lipnc1 : increase coefficient of the highest lipid concentration with physiological age area in the 5th (percent)
#' @param beta : lipid synthesis rate (j-1)
#' @param delta : lipid degradation rate (d-1)
#' @param amW :
#' @param b0c : coefficient of the allometric equation linking mass and lipid-protein carcass
#' @param b1c :  exponent allometric equation linking mass and defatted protein carcass
#' @param b0nc : coefficient of the allometric equation linking mass and lipid-protein 5th district
#' @param b1nc : exponent allometric equation linking mass and lipid-protein 5th district
#' @param c0 : coefficient of the allometric equation between live weight and live weight empty
#' @param c1 : exponent allometric equation linking body weight and live weight empty
#' @param energie : Metabolizable energy available
#' @param PVi : initial liveweight
#' @param duration : duration of simulation
#' @return matrix with ProtC,LipC,ProtNC,LipNC,PV
#' @export
#Function
carcass.EMI.model <- function(protcmax,protncmax,alphac,alphanc,gammac,gammanc,lip0,lipc1,lipnc1,beta,delta,amW,b0c,b1c,b0nc,b1nc,c0,c1,energie,PVi,duration)
  {
   #Intialize variables
    # 5 states variables, as 5 vectors initialized to NA
      #ProtC : Proteins in carcass tissues
ProtC=rep(NA,1,duration)
      #LipC : Lipids in carcass tissues
LipC=rep(NA,1,duration)
      #ProtNC : Proteins in non-carcass tissues
ProtNC=rep(NA,1,duration)
      #LipNC : Lipids in non-carcass tissues
LipNC=rep(NA,1,duration)
      #PV : body weight
PV=rep(NA,1,duration)

   # Initialization of state variables

PVVi=(1/c0)*PVi^(1/c1)
PoidsCi=0.8*PVVi
PoidsNCi=PVVi-PoidsCi
MDCi=0.8*PoidsCi
ProtCi=(1/b0c)*MDCi^(1/b1c)
LipCi=PoidsCi-MDCi
MDNCi=0.7*PoidsNCi
ProtNCi=(1/b0nc)*MDNCi^(1/b1nc)
LipNCi=PoidsNCi-MDNCi


ProtC[1]=ProtCi
LipC[1]=LipCi
ProtNC[1]=ProtNCi
LipNC[1]=LipNCi

   # Simulation loop
     for (t in 1:duration)
       {

         #calcul of body weight PV

          MDC=b0c*ProtC[t]^b1c
          PoidsC=LipC[t]+MDC
          MDNC=b0nc*ProtNC[t]^b1nc
          PoidsNC=LipNC[t]+MDNC

          PVV=PoidsC+PoidsNC
          PV[t]=c0*PVV^c1

          EMI=energie$y[t]
          CPM=amW*PV[t]^0.75

         #Carcass tissues
           #Proteins
             #Synthesis
          ProtCSyn=alphac*ProtC[t]*log(protcmax/ProtC[t])*(EMI/(CPM+EMI))
             #Degradation
          ProtCDeg=gammac*ProtC[t]*log(protcmax/ProtC[t])

           #Lipids
          LipCmax=(lip0+lipc1*(ProtC[t]/protcmax))*PoidsC
             #Synthesis
          LipCSyn=beta*LipC[t]*log(LipCmax/LipC[t])*(EMI/(CPM+EMI))
             #Degradation
          LipCDeg=delta*LipC[t]*log(LipCmax/LipC[t])

        #Non-carcass tissues
          #Proteins
            #Synthesis
          ProtNCSyn=alphanc*ProtNC[t]*log(protncmax/ProtNC[t])*(EMI/(CPM+EMI))
            #Degradation
          ProtNCDeg=gammanc*ProtNC[t]*log(protncmax/ProtNC[t])
           #Lipids
          LipNCmax=(lip0+lipnc1*(ProtNC[t]/protncmax))*PoidsNC
             #Synthesis
          LipNCSyn=beta*LipNC[t]*log(LipNCmax/LipNC[t])*(EMI/(CPM+EMI))
             #Degradation
          LipNCDeg=delta*LipNC[t]*log(LipNCmax/LipNC[t])

         # Calculate rates of change of state variables (dPC,dLC,dPNC,dLNC)

          dPC=ProtCSyn-ProtCDeg
          dLC=LipCSyn-LipCDeg
          dPNC=ProtNCSyn-ProtNCDeg
          dLNC=LipNCSyn-LipNCDeg

         # Update state variables

          ProtC[t+1]=ProtC[t]+dPC
          LipC[t+1]=LipC[t]+dLC
          ProtNC[t+1]=ProtNC[t]+dPNC
          LipNC[t+1]=LipNC[t]+dLNC

        }
  # End of simulation loop

  results=data.frame(time=c(1:duration),ProtC=ProtC[1:duration],LipC=LipC[1:duration],ProtNC=ProtNC[1:duration],LipNC=LipNC[1:duration],PV=PV[1:duration])
  return(results)
  }
################################################################################
#' @title The Carcass model function for use with carcass.EMI.simule
#' @description see carcass.EMI.model for model description.
#' @param param : a vector of parameters
#' @param energie : Metabolizable energy available
#' @param PVi : initial liveweight
#' @param duration : duration of simulation
#' @return data.frame with PV,ProtC,ProtNC,LipC,LipNC
#' @export
carcass.EMI.model2=function(param,energie,PVi,duration)
{
 # use carcass.model.EMI function to run the model
return(carcass.EMI.model(param[1],param[2],param[3],param[4],param[5],param[6],param[7],param[8],param[9],param[10],param[11],param[12],param[13],param[14],param[15],param[16],param[17],param[18],energie,PVi,duration))
}

################################################################################
#' @title Wrapper function to the Carcass model for multiple sets of parameter values
#' @param X : parameter matrix
#' @param energy : Metabolizable energy available
#' @param PVi : initial liveweight
#' @param duration : duration of simulation
#' @return data.frame with PV,ProtC,ProtNC,LipC,LipNC
#' @export
carcass.EMI.simule = function(X, energy,PVi,duration){
# output : all
sim <- lapply(1:dim(X)[1], function(id) cbind(id,carcass.EMI.model2(X[id,],energy,PVi,duration)))
sim=do.call(rbind, sim)
return(sim)
}
################################################################################
#' @title Wrapper function to run Carcass model on several animals with different conditions
#' @param param : a vector of parameters
#' @param list_individuals : list of individuals
#' @param energy : Metabolizable energy available for all individuals
#' @param init_condition : initial condition for all individuals
#' @return data.frame with id, ration ,duration, day, PV,ProtC,ProtNC,LipC,LipNC
#' @export
carcass.EMI.multi <- function(param, list_individuals,energy,init_condition)
{
sim=data.frame()
  for (id in list_individuals)
  {
  energyI=subset(energy,energy$Individu==id)
  energyInd=energyI$energie/7*10
  day=energyI$time*7-6
  energy_interpol<-approx(day,energyInd,xout=1:max(day),method = "constant")
  duration=length(energy_interpol$y)

  model<-carcass.EMI.model2(param,energy_interpol,init_condition[init_condition$Individu==id,"Pvi"],duration)
  model=data.frame(cbind(individu=rep(id,max(day)),Ration=rep(energyI$Ration[1],max(day)),model))
  sim=rbind(sim,model)
  }
return(sim)
}
################################################################################
#' @title Define values of the parameters for the Carcass model
#' @param full : if TRUE, return the full description of distribution(default = FALSE)
#' @return matrix with parameter values (nominal, binf, bsup). A data.frame if full=TRUE
#' @examples carcass.define.param(full=TRUE)
#' @export
carcass.define.param = function(full=FALSE)
{
param <-
structure(c(115, 104, 126, 115, 11.5, 55, 50, 60, 55, 5, 0.016,
0.015, 0.016, 0.0155, 5e-04, 0.032, 0.031, 0.032, 0.0315, 5e-04,
0.0025, NA, NA, NA, NA, 0.009, NA, NA, NA, NA, 0.12, 0.11, 0.13,
0.12, 0.01, 0.13, 0.12, 0.14, 0.13, 0.01, 0.13, 0.12, 0.14, 0.13,
0.01, 0.1, 0.09, 0.11, 0.1, 0.01, 0.025, 0.0225, 0.0275, 0.025,
0.0025, 1.5, 1.2, 1.8, 1.5, 0.3, 5.665, NA, NA, NA, NA, 0.949,
NA, NA, NA, NA, 7.224, NA, NA, NA, NA, 0.85, NA, NA, NA, NA,
1.3022, NA, NA, NA, NA, 0.9766, NA, NA, NA, NA), .Dim = c(5L,
18L), .Dimnames = list(c("nominal", "binf", "bsup", "mean", "sd"
), c("protcmax", "protncmax", "alphac", "alphanc", "gammac",
"gammanc", "lip0", "lipc1", "lipnc1", "beta", "delta", "amW",
"b0c", "b1c", "b0nc", "b1nc", "c0", "c1")), distribution = c("normal",
"normal", "uniform", "uniform", "constant", "constant", "uniform",
"uniform", "uniform", "uniform", "uniform", "normal", "constant",
"constant", "constant", "constant", "constant", "constant"), remark = c("Robelin 1986",
"Robelin 1986", "Agabriel", "Agabriel", "not to vary (Agabriel)",
"not to vary (Agabriel)", "Hoch 2004 nb2", "Hoch 2004 nb2", "Hoch 2004 nb2",
"Danfaer (1991)", "Danfaer (1991)", "Hoch 2004", "ajusted (Robelin 1986, thesis)",
"ajusted (Robelin 1986, thesis)", "ajusted (Robelin 1986, thesis)",
"ajusted (Robelin 1986, thesis)", "ajusted (Garcia 2007)", "ajusted (Garcia 2007)"
))
if (full){
return(cbind(as.data.frame(t(param)), distribution=attr(param,"distribution"),remark=attr(param,"remark")))
} else {return(as.matrix(param))}
}

# End of file 

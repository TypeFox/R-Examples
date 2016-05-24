################################ FUNCTIONS #####################################
#' @title The Carcass (growth of beef cattle) model
#' @description \strong{Model description.}
#'
#'
#'
#' The model is defined by 20 equations, with a total of 19 parameters for the described process.
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
#' @param k : Parameter coefficient between the half-saturation of the Michaelis-Menten equation of the metabolic weight (MJ.kg^0.75)
#' @param b0c : coefficient of the allometric equation linking mass and lipid-protein carcass
#' @param b1c :  exponent allometric equation linking mass and defatted protein carcass
#' @param b0nc : coefficient of the allometric equation linking mass and lipid-protein 5th district
#' @param b1nc : exponent allometric equation linking mass and lipid-protein 5th district
#' @param c0 : coefficient of the allometric equation between live weight and live weight empty
#' @param c1 : exponent allometric equation linking body weight and live weight empty
#' @param cem :
#' @param duration : duration of simulation
#' @return data.frame with ProtC,LipC,ProtNC,LipNC,PV
#' @export
carcass.model <- function(protcmax,protncmax,alphac,alphanc,gammac,gammanc,lip0,lipc1,lipnc1,beta,delta,k,b0c,b1c,b0nc,b1nc,c0,c1,cem,duration)
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
ProtC[1]=30
LipC[1]=15
ProtNC[1]=15
LipNC[1]=8

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

          EMI=cem*(0.0157*(PV[t]^0.9)+3.3161)
          CPM=k*PV[t]^0.75
          
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
          
         # Uptade state variables
         
          ProtC[t+1]=ProtC[t]+dPC
          LipC[t+1]=LipC[t]+dLC
          ProtNC[t+1]=ProtNC[t]+dPNC
          LipNC[t+1]=LipNC[t]+dLNC

        }
  # End simulation loop
  
  results=data.frame(time=c(1:duration),ProtC=ProtC[1:duration],LipC=LipC[1:duration],ProtNC=ProtNC[1:duration],LipNC=LipNC[1:duration],PV=PV[1:duration])
  return(results)
  }
# End of file

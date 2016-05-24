################################ FUNCTIONS #####################################
# Contribution of Juliette Adrian, Master2 internship, january-jully 2013
#' @title The Lactation model with milking machine
#' @description \strong{Model description.}
#' This model is a model of lactating mammary glands of cattle described by Heather et al. (1983). This model was then inspired more complex models based on these principles.
#' This model simulates the dynamics of the production of cow's milk.
#' the system is represented by 6 state variables: change in hormone levels (H), the production and loss of milk secreting cells (CS), and removing the secretion of milk (M), the average quantity of milk contained in the animal (Mmean), the amount of milk removed (RM) and yield (Y).
#' The model has a time step dt = 0.001 for milking machines.
#' The model is defined by a few equations, with a total of twenty parameters for the described process.
#' @param cu : number of undifferentiated cells
#' @param kdiv : cell division rate, Michaelis-Menten constant
#' @param kdl : constant degradation of milk
#' @param kdh : rate of decomposition of the hormone
#' @param km :  constant secretion of milk
#' @param ksl : milk secretion rate, Michaelis-Menten constant
#' @param kr : average milk constant
#' @param ks : rate of degradation of the basal cells
#' @param ksm : constant rate of degradation of milk secreting cells
#' @param mh : parameter
#' @param mm : storage Capacity milk the animal
#' @param p : parameter
#' @param mum : setting the maximum rate of cell division
#' @param rma : parameter of milk m (t) function
#' @param t1 : parameter of milk m (t) function
#' @param t2 : parameter of milk m (t) function
#' @param t3 : parameter of milk m (t) function
#' @param t4 : parameter of milk m (t) function
#' @param t5 : parameter of milk m (t) function
#' @param t6 : parameter of milk m (t) function
#' @param duration : duration of simulation
#' @param dt : time step
#' @param CSi : initial Number of secretory cells
#' @param Mi : initial Quantity of milk in animal (kg)
#' @return matrix with CS,M,Mmoy,RM
#' @export
lactation.model.machine=function(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6,duration,dt,CSi,Mi)
{


 # Initialize variables
 # 6 states variables, as 6 vectors initialized to NA
     # H : Hormone effector of cell division (kg/m³)
CS=rep(NA,(duration-1)/dt)
    # CS : Number of secretory cells
H=rep(NA,(duration-1)/dt)
    # M : Quantity of milk in animal (kg)
M=rep(NA,(duration-1)/dt)
    # Mmoy : Time average of M (kg)
Mmoy=rep(NA,(duration-1)/dt)
    # RM : Rate of removal of milk
RM=rep(NA,(duration-1)/dt)

  # Initialization of state variables
CS[1]=CSi
H[1]=1.0
M[1]=Mi
Mmoy[1]=0.0
RM[1]=0

i=1
  # Simulation loop
for (t in seq(0, duration, by = dt))
  {
    # Calculate milking function mach
    if (t%%1 > t1 & t%%1< t2 | t%%1 > t3 & t%%1< t4 | t%%1 > t5 & t%%1< t6)   # 3 pulses, for 2 pulses t5=t6=0
       {
        mach = rma
	     } else {
        mach = 0
       }

    # Calculate rates of change of state variables (dH,dCS,dM,dMmoy)
dH = - kdh * H[i] * dt
dCS = ( mum * (H[i]/(kdiv+H[i]))*cu - (ks + ksm*((Mmoy[i]/mh)^p/(1+(Mmoy[i]/mh)^p)))*CS[i] ) * dt
dM = (km * CS[i] * ((mm-M[i])/(mm-M[i]+ksl))-(M[i]/(kdl+M[i]))*mach) * dt
dMmoy = kr*(M[i]-Mmoy[i]) * dt

  # Uptade state variables
    H[i+1]= H[i] +dH
    CS[i+1]= CS[i] + dCS
    M[i+1]= M[i] + dM
    Mmoy[i+1]= Mmoy[i] + dMmoy

  # Removal of milk
    RM[i+1]=(M[i+1]/(kdl+M[i+1]))*mach
    i = i +1
    }
   # End simulation loop
  # conversion day to week
   day=seq(dt,duration,by=dt)
   week=day%/%7

results1=data.frame(M=M[1:(duration/dt)],Mmoy=Mmoy[1:(duration/dt)],CS=CS[1:(duration/dt)],RM=RM[1:(duration/dt)],day=day,week=week)
# mean by week
result = by(results1[,c("week","M","Mmoy","CS","RM")],results1$week, function(x) apply(x,2,mean))

results2 = matrix(unlist(result),ncol=5, byrow=TRUE,dimnames=list(NULL, c("week","M","Mmoy","CS","RM")) )
return(results2)
}
# end of file

################################################################################
#' @title The Lactation model for use with lactation.machine.simule
#' @description see lactation.calf.model for model description.
#' @param param : a vector of parameters containning (cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rma,t1,t2,t3,t4,t5,t6)(see lactation.model.machine)
#' @param duration : duration of simulation
#' @param dt : time step
#' @param CSi : initial Number of secretory cells
#' @param Mi : initial Quantity of milk in animal (kg)
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @export
lactation.model.machine2=function(param,duration,dt,CSi,Mi)
{
 # use lactation.model.machine function to run the model
return(lactation.model.machine(param["cu"],param["kdiv"],param["kdl"],param["kdh"],param["km"],param["ksl"],param["kr"],param["ks"],param["ksm"],param["mh"],param["mm"],param["p"],param["mum"],param["rma"],param["t1"],param["t2"],param["t3"],param["t4"],param["t5"],param["t6"],duration,dt,CSi,Mi))
}
# End of file

################################ FUNCTIONS #####################################
# Contribution of Juliette Adrian, Master2 internship, january-jully 2013
#' @title The Lactation model
#' @description \strong{Model description.}
#' This model is a model of lactating mammary glands of cattle described by Heather et al. (1983). This model was then inspired more complex models based on these principles.
#' This model simulates the dynamics of the production of cow's milk.
#' the system is represented by 6 state variables: change in hormone levels (H), the production and loss of milk secreting cells (CS), and removing the secretion of milk (M), the average quantity of milk contained in the animal (Mmean), the amount of milk removed (RM) and yield (Y).
#' The model has a time step dt = 0.1 for regular consumption of milk by a calf.
#' The model is defined by a few equations, with a total of fourteen parameters for the described process.
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
#' @param rc : parameter of milk m (t) function
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @examples lactation.calf.model2(lactation.define.param()["nominal",],300,0.1)
#' @export
lactation.calf.model <- function(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rc,duration,dt)
{
 # Initialize variables
 # 5 states variables, as 5 vectors initialized to NA
    # H : Hormone effector of cell division (kg/m³)
H=rep(NA,(duration-1)/dt)
    # CS : Number of secretory cells
CS=rep(NA,(duration-1)/dt)
    # M : Quantity of milk in animal (kg)
M=rep(NA,(duration-1)/dt)
    # Mmoy : Time average of M (kg)
Mmoy=rep(NA,(duration-1)/dt)
    # RM : Rate of removal of milk
RM=rep(NA,(duration-1)/dt)

 # Initialization of state variables
H[1]=1.0
CS[1]=520
M[1]=0.0
Mmoy[1]=0.0

i=1
 # Simulation loop
for (t in seq(0, duration, by = dt))
  {
 # Calculate rates of change of state variables (dH,dCS,dM,dMmoy)
    dH = - kdh * H[i] * dt
    dCS = (mum * (H[i]/(kdiv+H[i]))*cu - (ks + ksm*((Mmoy[i]/mh)^p/(1+(Mmoy[i]/mh)^p)))*CS[i] ) * dt
    dM = (km * CS[i] * ((mm-M[i])/(mm-M[i]+ksl))-(M[i]/(kdl+M[i]))*rc	) * dt
    dMmoy = kr*(M[i]-Mmoy[i]) * dt

 # Uptade state variables
    H[i+1]= H[i] +dH
    CS[i+1]= CS[i] + dCS
    M[i+1]= M[i] + dM
    Mmoy[i+1]= Mmoy[i] + dMmoy

  # removal of milk
    RM[i]=(M[i]/(kdl+M[i]))*rc

  i=i+1
  }
  # End simulation loop
  # conversion day to week
   day=seq(dt,duration,by=dt)
   week=day%/%7

results1=data.frame(M=M[1:(duration/dt)],Mmoy=Mmoy[1:(duration/dt)],CS=CS[1:(duration/dt)],RM=RM[1:(duration/dt)],day=day,week=week)
# mean by week
#result = by(,results1$week, mean)
result = by(results1[,c("week","M","Mmoy","CS","RM")],results1$week,function(x) apply(x,2,mean))
results2 = matrix(unlist(result),ncol=5, byrow=TRUE,dimnames=list(NULL, c("week","M","Mmoy","CS","RM")) )
return(results2)

#results=data.frame(CS=CS[1:(duration/dt)],M=M[1:(duration/dt)],Mmoy=Mmoy[1:(duration/dt)],RM=RM[1:(duration/dt)],day=seq(0.1,duration,by=dt),week=seq(0.1/7,duration/7,by=dt/7))
#return(results)
}
################################################################################
#' @title The Lactation model for use with lactation.calf.simule
#' @description see lactation.calf.model for model description.
#' @param param : a vector of parameters
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with CS, M, Mmoy, RM, day, week
#' @examples sim=lactation.calf.model2(lactation.define.param()["nominal",],6+2*7, 0.1)
#' @export
lactation.calf.model2 <- function(param,duration,dt){
 # use lactation.calf.model function to run the model
  return(lactation.calf.model(param["cu"],param["kdiv"],param["kdl"],param["kdh"],param["km"],param["ksl"],param["kr"],param["ks"],param["ksm"],param["mh"],param["mm"],param["p"],param["mum"],param["rc"],duration,dt))
}
################################################################################
#' @title Wrapper function to run the Lactation model for multiple sets of parameter values
#' @param X : parameter matrix
#' @param duration : duration of simulation
#' @param dt : time step
#' @return data.frame with : number of paramter vector (line number from X), week, CS, M, Mmoy, RM, day, week
#' @export
lactation.calf.simule = function(X, duration, dt){
# output : all
#sim <- apply(X,1,function(v) lactation.calf.model2(v,duration, dt))
#sim=do.call(rbind, sim)
sim <- lapply(1:dim(X)[1], function(id) cbind(id,lactation.calf.model2(X[id,],duration, dt)))
sim=do.call(rbind, sim)
return(sim)
}
################################################################################
#' @title Define values of the parameters for the Lactation model
#' @description values from Heather et al. (1983) for different scenarios
#' @param type : for which model version ? "calf" or "machine"
#' @return matrix with parameter values (nominal, binf, bsup)
#' @examples lactation.define.param()
#' @export
lactation.define.param <- function(type="calf")
{
# nominal, binf, bsup
#cu : number of undifferentiated cells (Unit ?)
cu=c(1000, NA, NA)
#kdiv : cell division rate, Michaelis-Menten constant (Unit ?)
kdiv=c(0.2, NA, NA)
#kdl : constant degradation of milk (Unit ?)
kdl=c(4.43, NA, NA)
#kdh : rate of decomposition of the hormone (Unit ?)
kdh=c(0.01, NA, NA)
#km : constant secretion of milk ksl : milk secretion rate, Michaelis-Menten constant (Unit ?)
km=c(0.005, NA, NA)
#
ksl=c(3.0, NA, NA)
#kr : average milk constant (Unit ?)
kr=c(0.048, NA, NA)
#ks : rate of degradation of the basal cells (Unit ?)
ks=c(0.1, NA, NA)
#ksm : constant rate of degradation of milk secreting cells (Unit ?)
ksm=c(0.2, NA, NA)
#mh : parameter mm : storage Capacity milk the animal (Unit ?)
mh=c(27, NA, NA)
#
mm=c(30, NA, NA)
#p : parameter mum : setting the maximum rate of cell division (Unit ?)
p=c(10, NA, NA)
#
mum=c(1, NA, NA)
#for calf
#rc : parameter of milk m (t) function
rc=c(40, NA, NA)
#for machine
#rma : parameter of milk m (t) function (Unit ?)
rma=c(NA, NA, NA)

CSi=c(NA, NA, NA)
Mi=c(NA, NA, NA)
if (type=="calf"){param<-data.frame(cu,kdiv,kdl,kdh,km,ksl,kr,ks,ksm,mh,mm,p,mum,rc)}
else {if (type=="machine") {param<-data.frame("TODO")}}

row.names(param)<-c("nominal","binf","bsup")
return(as.matrix(param))
}
# end of file

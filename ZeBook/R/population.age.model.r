# ZeBook
# Francois Brun
# version 2013-03-07
################################ FUNCTIONS #####################################
#' @title The PopulationAge model (Population Dynamics with Age Classes)
#' @description Population Dynamics Model with Age Classes for an insect
#' @param rb : eggs laid per adult per unit area (day-1)
#' @param rE : eggs hatch (day-1)
#' @param r12 : relative rate L1->L2 (day-1)
#' @param r23 : relative rate L2->L3 (day-1)
#' @param r34 : relative rate L3->L4 (day-1)
#' @param r4P : relative rate L4->P (day-1)
#' @param rPA : relative rate P->A (day-1)
#' @param mE : relative mortality rate of egg (day-1)
#' @param m1 : relative mortality rate of larvae L1 (day-1)
#' @param m2 : relative mortality rate of larvae L2 (day-1)
#' @param m3 : relative mortality rate of larvae L3 (day-1)
#' @param m4 : relative mortality rate of larvae L4 (day-1)
#' @param mP : relative mortality rate of purpae (day-1)
#' @param mA : relative mortality rate of adult L1 (day-1)
#' @param iA : input rate of adult (unit.day-1)
#' @param duration : simulation duration
#' @param dt : time step for integration
#' @return data.frame with values for state variables for each time step.
#' @export
population.age.model = function(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=100,dt=1){
  # states variables, vectors initialized to NA
  # E : egg stage. homogenous population (density) (number per ha)
  # L1 : larvae1 stage. homogenous population (density) (number per ha)
  # L2 : larvae2 stage. homogenous population (density) (number per ha)
  # L3 : larvae3 stage. homogenous population (density) (number per ha)
  # L4 : larvae4 stage. homogenous population (density) (number per ha)
  # P : pupae stage. homogenous population (density) (number per ha)
  # A : adult stage. homogenous population (density) (number per ha)
  
  # V : matrix of state variable (one per column)
  V = matrix(NA,ncol=7,nrow=duration/dt+1, dimnames=list(NULL,c("E","L1","L2","L3","L4","P","A")))
  
  # initiation of state variable
  V[1,]<- c(5,0,0,0,0,0,0)
  
  # Simulation loop
  for (k in 1:(duration/dt)){
    # Calculate rates of change of state variables
    dE= (rb*V[k,"A"] - rE*V[k,"E"] - mE*V[k,"E"])*dt
    dL1= (rE*V[k,"E"] - r12*V[k,"L1"] - m1*V[k,"L1"])*dt
    dL2= (r12*V[k,"L1"] - r23*V[k,"L2"] - m2*V[k,"L2"])*dt
    dL3= (r23*V[k,"L2"] - r34*V[k,"L3"] - m3*V[k,"L3"])*dt
    dL4= (r34*V[k,"L3"] - r4P*V[k,"L4"] - m4*V[k,"L4"])*dt
    dP= (r4P*V[k,"L4"] - rPA*V[k,"P"] - mP*V[k,"P"])*dt
    dA= (rPA*V[k,"P"] - mA*V[k,"A"] + iA)*dt
    
    # vector of rates of change
    dV=c(dE,dL1,dL2,dL3,dL4,dP,dA)
    
    # Update state variables 
    V[k+1,]<- V[k,] + dV
  }
  # End simulation loop
  return(round(as.data.frame(cbind(time=(1:(duration/dt+1))*dt-dt, V)),10))    
}
################################################################################
#' @title The PopulationAge model (Population Dynamics with Age Classes) - matrix form
#' @description Population Dynamics Model with Age Classes for an insect
#' Exactly the same model as population.age.model, but written as a matrix computation. It's possible for this model. It's really more efficient and reduce computer time by 6!
#' 7 states variables
#' E : egg stage. homogenous population (density) (number per ha)
#' L1 : larvae1 stage. homogenous population (density) (number per ha)
#' L2 : larvae2 stage. homogenous population (density) (number per ha)
#' L3 : larvae3 stage. homogenous population (density) (number per ha)
#' L4 : larvae4 stage. homogenous population (density) (number per ha)
#' P : pupae stage. homogenous population (density) (number per ha)
#' A : adult stage. homogenous population (density) (number per ha)
#' @param rb : eggs laid per adult per unit area (day-1)
#' @param rE : eggs hatch (day-1)
#' @param r12 : relative rate L1->L2 (day-1)
#' @param r23 : relative rate L2->L3 (day-1)
#' @param r34 : relative rate L3->L4 (day-1)
#' @param r4P : relative rate L4->P (day-1)
#' @param rPA : relative rate P->A (day-1)
#' @param mE : relative mortality rate of egg (day-1)
#' @param m1 : relative mortality rate of larvae L1 (day-1)
#' @param m2 : relative mortality rate of larvae L2 (day-1)
#' @param m3 : relative mortality rate of larvae L3 (day-1)
#' @param m4 : relative mortality rate of larvae L4 (day-1)
#' @param mP : relative mortality rate of purpae (day-1)
#' @param mA : relative mortality rate of adult L1 (day-1)
#' @param iA : input rate of adult (unit.day-1)
#' @param duration : simulation duration
#' @param dt : time step for integration
#' @return data.frame with values for state variables for each time step.
#' @export
population.age.matrix.model = function(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=100,dt=1){

  # V : matrix of state variable (one per column), initialized to NA
  V = matrix(NA,ncol=7,nrow=duration/dt+1, dimnames=list(NULL,c("E","L1","L2","L3","L4","P","A")))
  
  # initiation of state variable
  V[1,]<- c(5,0,0,0,0,0,0)
  
  # defining matrix transition
  M= matrix(c(-rE-mE, 0,0,0,0,0,rb,
  rE, -r12-m1,0,0,0,0,0,
  0,r12,-r23-m2,0,0,0,0,
  0,0,r23,-r34-m3,0,0,0,
  0,0,0,r34,-r4P-m4,0,0,
  0,0,0,0,r4P,-rPA-mP,0,
  0,0,0,0,0,rPA,-mA),ncol=7,nrow=7, byrow=TRUE)
  # input/output rate matrix
  IN = matrix(c(0, 0,0,0,0,0,iA),ncol=1,nrow=7, byrow=FALSE) 
  # Simulation loop
  for (k in 1:(duration/dt)){
    # Calculate rates of change of state variables
    dV= (M %*% V[k,] + IN)*dt
    # Update state variables 
    V[k+1,]<- V[k,] + dV
  }
  # End simulation loop
  return(round(as.data.frame(cbind(time=(1:(duration/dt+1))*dt-dt, V)),10))    
}
################################################################################
#' @title The PopulationAge model (Population Dynamics with Age Classes) - ode form
#' @description Population Dynamics Model with Age Classes for an insect
#' Exactly the same model as population.age.model, but written as an ordinary differential equation system (ode) with deSolve package.
#' 7 states variables
#' E : egg stage. homogenous population (density) (number per ha)
#' L1 : larvae1 stage. homogenous population (density) (number per ha)
#' L2 : larvae2 stage. homogenous population (density) (number per ha)
#' L3 : larvae3 stage. homogenous population (density) (number per ha)
#' L4 : larvae4 stage. homogenous population (density) (number per ha)
#' P : pupae stage. homogenous population (density) (number per ha)
#' A : adult stage. homogenous population (density) (number per ha)
#' @param rb : eggs laid per adult per unit area (day-1)
#' @param rE : eggs hatch (day-1)
#' @param r12 : relative rate L1->L2 (day-1)
#' @param r23 : relative rate L2->L3 (day-1)
#' @param r34 : relative rate L3->L4 (day-1)
#' @param r4P : relative rate L4->P (day-1)
#' @param rPA : relative rate P->A (day-1)
#' @param mE : relative mortality rate of egg (day-1)
#' @param m1 : relative mortality rate of larvae L1 (day-1)
#' @param m2 : relative mortality rate of larvae L2 (day-1)
#' @param m3 : relative mortality rate of larvae L3 (day-1)
#' @param m4 : relative mortality rate of larvae L4 (day-1)
#' @param mP : relative mortality rate of purpae (day-1)
#' @param mA : relative mortality rate of adult L1 (day-1)
#' @param iA : input rate of adult (unit.day-1)
#' @param duration : simulation duration
#' @param dt : time step for integration
#' @param  method : integration method (euler, rk4,...)
#' @return data.frame with values for state variables for each time step.
#' @export
population.age.model.ode = function(rb=3.5,mE=0.017,rE=0.172,m1=0.060,r12=0.217,m2=0.032,r23=0.313,m3=0.022,r34=0.222,m4=0.020,r4P=0.135,mP=0.020,rPA=0.099,mA=0.027,iA=0,duration=100, dt=1, method="euler"){
  # states variables
  # E : egg stage. homogenous population (density) (number per ha)
  # L1 : larvae1 stage. homogenous population (density) (number per ha)
  # L2 : larvae2 stage. homogenous population (density) (number per ha)
  # L3 : larvae3 stage. homogenous population (density) (number per ha)
  # L4 : larvae4 stage. homogenous population (density) (number per ha)
  # P : pupae stage. homogenous population (density) (number per ha)
  # A : adult stage. homogenous population (density) (number per ha)
  E0=5
  L10=0
  L20=0
  L30=0
  L40=0
  P0=0
  A0=0
  # definiting the model as an ordinary differential equation system
  predator.prey.ode <- function(Time, State, Pars) {
      with(as.list(c(State, Pars)), {
          dE= (rb*A - rE*E - mE*E)*dt
          dL1= (rE*E - r12*L1 - m1*L1)*dt
          dL2= (r12*L1 - r23*L2 - m2*L2)*dt
          dL3= (r23*L2 - r34*L3 - m3*L3)*dt
          dL4= (r34*L3 - r4P*L4 - m4*L4)*dt
          dP= (r4P*L4 - rPA*P - mP*P)*dt
          dA= (rPA*P - mA*A + iA)*dt
          return(list(c(dE, dL1, dL2, dL3, dL4, dP, dA)))
      })
  }
  sim = ode(y=c(E=E0,L1=L10,L2=L20,L3=L30,L4=L40,P=P0,A=A0), times=seq(0,duration,by=dt), func=predator.prey.ode, parms=c(rb,mE,rE,m1, r12,m2,r23,m3,r34,m4,r4P,mP,rPA,mA,iA), method = rkMethod(method))
  return(as.data.frame(sim))
}
# end of file

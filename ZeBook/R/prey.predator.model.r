# Francois Brun, 2013-05-15
################################ FUNCTIONS #####################################
#' @title The PredatorPrey model (Predator-Prey Lotka-Volterra with logistic equation for prey)
#' @description Predator-Prey Lotka-Volterra model (with logistic prey)
#' @param grH : relative rate of prey population growth
#' @param kH : environment carrying capacity for prey (number per ha)
#' @param mrH : maximum predation rate (number per predator and per prey per day)
#' @param eff : efficiency, growth of predator population depending on predation (-) 
#' @param mrA : mortality of predator (-)
#' @param  H0 : size of population of prey, at time 0
#' @param  A0 : size of population of predator, at time 0
#' @param  duration : simulation duration
#' @param  dt : time step for integration
#' @param  method : integration method
#' @return data.frame with daily H and A
#' @export
predator.prey.model = function(grH=1, kH=10, mrH=0.2, eff=0.5, mrA=0.2, H0=1,A0=2,duration=200,dt=1, method="euler"){
# 2 states variables
# H : prey, Aphids, homogenous population (density) (number per ha)
# A : predators, Ladybeetles, homogenous population (density)  (number per ha)
# definiting the model as an ordinary differential equation system    
  predator.prey.ode <- function(Time, State, Pars) {
      with(as.list(c(State, Pars)), {
          dH <- grH*H*(1-H/kH) - mrH*H*A
          dA <- mrH*H*A*eff - mrA*A
          return(list(c(dH, dA)))
      })
  }
  sim = ode(y=c(H=H0,A=A0), times=seq(0,duration,by=dt), func=predator.prey.ode, parms=c(mrH,grH,mrA,eff,kH), method = rkMethod(method))
  return(as.data.frame(sim))
}
# end of file

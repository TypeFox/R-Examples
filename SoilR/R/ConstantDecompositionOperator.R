#
# vim:set ff=unix expandtab ts=2 sw=2:
correctnessOfDecompOp=function(object)
  A=object@mat
  
setClass(# constant decomposition operator 
    Class="ConstLinDecompOp",
    contains="DecompOp",
    slots=list( mat="matrix")
)
setMethod(
     f="initialize",
     signature="ConstLinDecompOp",
     definition=function #internal constructor 
     ### This mehtod is intendet to be used internally. It may change in the future.
     ### For user code it is recommended to use one of the generic constructor \code{ConstLinDecompOp} instead.
     (.Object,mat=matrix())
     {
        .Object@mat=mat
     return(.Object)
     }
)
############################constructors####################
setMethod(
      f="ConstLinDecompOp",
      ### 
      signature=c(mat="matrix"),
      definition=function # construct from matric
      ### This method creates a ConstLinDecompOp from a matrix
      ### The operator is assumed to act on the vector of carbon stocks
      ### by multiplication of the (time invariant) matrix from the left.
      (mat){
      return(new("ConstLinDecompOp",mat=mat))
     }
)

############################methods####################
setMethod(
    f="getFunctionDefinition",
    signature="ConstLinDecompOp",
    definition=function # creates a constant timedependent function and returns it
      ### The method creates a timedependent function from the existing matrix describing the operator 
    (object){
      return(function(t){object@mat})
    }
)
setMethod(
    f="getTimeRange",
    signature="ConstLinDecompOp",
    definition=function # return an (infinite) time range since the operator is constant
    ### some functions dealing with DecompOps in general rely on this
    ### so we have to implement it even though the timerange is always the same: (-inf,inf)
    (object)
    {
        return( c("t_min"=-Inf,"t_max"=Inf))
    }
)
setMethod(
      f="BoundLinDecompOp",
      signature=c(map="ConstLinDecompOp"),
      definition=function # convert a ConstLinDecompOp to a BoundLinDecompOp
      ### The method creates a BoundLinDecompOp consisting of a constant time dependent function 
      ### and the limits of its domain (starttime and endtime) set to -Inf and Inf respectively
      (map){
      f=getFunctionDefinition(map)
      return(BoundLinDecompOp(starttime=-Inf,endtime=Inf,map=f))
     }
)
setMethod(
  f= "getMeanTransitTime",
    signature= "ConstLinDecompOp",
    definition=function #compute the mean transit time 
      ### This method computes the mean transit time for the linear time invariant system 
      ### that can be constructed from the given operator and input distribution.
      ### 
      ### It relies on the mehtod \code{getTransitTimeDistributionDensity} using the same arguments.
      ##details<< To compute the mean transit time for the distribution we have to compute the integral
      ## \deqn{
      ## \bar{T} = \int_0^\infty T \cdot S_r\left( \frac{\vec{I}}{I},0,T\right)\; dT
      ## }
      ## for the numerically computed density.
      ## To avoid issues with numerical integration  we dont use \eqn{\infty}{\infty} as upper limit but cut off the integragion interval prematurely.
      ## For this purpose we calculate a maximum response time of the system as \cite{Lasaga}
      ## \deqn{
      ## \tau_{cycle} = \frac{1}{|\min(\lambda_i)|}
      ## }
      ##  where \eqn{\lambda_i}{\lambda_i} are non-zero eigenvalues of the matrix 
      ## \eqn{{\bf A}}{{\bf A}}. 
      ##references<< 
      ## Lasaga, A.: The kinetic treatment of geochemical cycles, Geochimica et
      ## Cosmochimica Acta, 44, 815 -- 828, doi{10.1016/0016-7037(80)90263-X}, 1980.
      (object,inputDistribution){

      f=getFunctionDefinition(object)
      g=function(t){spectralNorm(f(t))}
      t_max=function(t_end){
          t_step=t_end/10
          t=seq(0,t_end,t_step)
          norms=sapply(t,g)
          tm=100*max(norms)
	  print(paste("tm=",tm))
	  return(tm)
      } 
      t_end=20
      print(paste("t_end=",t_end,sep=""))
      t_end_new=t_max(t_end)
      print(paste("t_end_new_before while=",t_end_new))
      while(t_end_new>t_end){
          print(t_end)
	  t_end=t_end_new
	  t_end_new=t_max(t_end)
      }
      print("after while")
      longTailEstimate=t_end
      subd=10000
      t_step=t_end/subd
      t=seq(0,t_end,t_step)
      shortTailEstimate=min(sapply(t,g))
      
      ttdd=getTransitTimeDistributionDensity(object,inputDistribution,t)
      #print(paste("ttdd=",ttdd))
      #print(paste("t=",t))
      #print(paste("ttdd*t=",ttdd*t))
      meanTimeRiemann=sum(ttdd*t)*t_step
      int2=splinefun(t,ttdd*t)
      meanTimeIntegrate=integrate(int2,0,t_end,subdivisions=subd)[["value"]] 
      print(paste("meanTimeRiemann=",meanTimeRiemann))
      print(paste("meanTimeIntegrate=",meanTimeIntegrate))
      #meanTime_s=integrate(integrand,0,shortTailEstimate,subdivisions=subd)[["value"]] 
      # here we must first check if the two values differ significantly
      #return(meanTimeRiemann)
      return(meanTimeIntegrate)
   }
)
setMethod(
   f= "getTransitTimeDistributionDensity",
      signature= "ConstLinDecompOp",
      definition=function #compute the TransitTimeDistributionDensity 
      ### This mehtod computes the probability density of the transit time of the linear time invariant system
      ### that can be constructed from the given operator and input distribution.
      ##details<< In a forthcoming paper \cite{SoilRv1.2} 
      ## we derive the algorithm used in this implementation
      ## under the assumption of steady conditions having prevailed infinitely.
      ## We arrive at a formulation well known from the literature about 
      ## time invariant linlear systems, cited e.g. in \cite{ManzoniJGR}.\cr
      ## The somehow amazing result is that the weight of the transit time density \eqn{\psi(T)}{\psi(T)} 
      ## for a \emph{transit time} \eqn{T}{T} 
      ## for the steady state system is identical to the output \eqn{O(T)}{O(T)} 
      ## observed at time \eqn{T}{T} of a \emph{different} system which started with a normalized impulsive input 
      ## \eqn{\frac{\vec{I}}{I}}{\frac{\vec{I}}{I}} at time 
      ## \eqn{T=0}{T=0},
      ## where \eqn{I=\sum_{k=1}^m i_k} is the cumulative input flux to all pools.
      ## \cr
      ## This fact simpliefies the computation considerably.
      ## Translated into the language of an ode solver an impulsive input becomes a start vector \eqn{\frac{\vec{I}}{I}}{\frac{\vec{I}}{I}} 
      ## at time \eqn{T=0}{T=0} 
      ## and \eqn{O(T)}{O(T)} the respiration related to the solution of the initial value problem 
      ## observed at time \eqn{T}{T}. 
      ## \deqn{
      ## \psi(T)=S_r \left( \frac{\vec{I}}{I},0,T\right)
      ## }{}
      ## Note that from the perspective of the ode solver \eqn{S_r}{S_r} depends on the decomposition operator and the distribution of the input among the pools only.
      ## It is therefor possible to implement the transit time distribution as a function of the decomposition operator and the fixed input flux distribution.
      ## To insure steady state conditions the decomposition operator is not allowed to be a true function of time.
      ## We therefor implement the method only for the subclass 
      ## \code{ConstLinDecompOp} 
      ## \cr Remark:\cr
      ## The decision to implement this method for \code{transitTimeDensity} especially for 
      ## objects of class \code{ConstLinDecompOp}
      ## reflects the fact that an arbitrary  model in SoilR is by no means bound to end up in steady state. To insure this we would have to ignore the input part of a user created model which would be confusing. 
      ## \cr Remark:\cr
      ## In future versions of SoilR it will be possible to compute a dynamic, time dependent transit time distribution 
      ## for objects of class \code{ Model}
      ## with a time argument specifying for which time the distribution is sought. 
      ## The steady state results computed here could than be reproduced 
      ## with the user responsible for providing a model actually reaching it. 
      ##references<<
      ##  Manzoni, S., Katul, G.~G., and Porporato, A.: Analysis of soil carbon transit
      ##  times and age distributions using network theories, J. Geophys. Res., 114,
      ##  
      (object,inputDistribution,times){
      # we set the initial values to the value provided by the inputdistribution
      sVmat=inputDistribution
      n=length(inputDistribution)
      # we provide a zero inputflux
      inputFluxes=BoundInFlux(
        map=function(t0){matrix(nrow=n,ncol=1,0)},
        starttime= -Inf, 
        endtime=+Inf 
      ) 
      #we create a model 
      mod=GeneralModel(times,object,sVmat,inputFluxes)
      R=getReleaseFlux(mod)
      TTD=rowSums(R)
      return(TTD)
   }
)

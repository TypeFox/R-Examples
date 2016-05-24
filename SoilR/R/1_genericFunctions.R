#
# vim:set ff=unix expandtab ts=2 sw=2:
setGeneric(
    name="Delta14C",
    def=function( # Converts its argument to a Delta14C representation
    ### The function returns an object of the same type as its input,
    ### which can be of different type.
    ### Have a look at the methods for details.
    F ##<< an object that contains data and a format description.  So it can be converted into the AbsoluteFractionModern format if a conversion is implemented.
    ){
        standardGeneric("Delta14C")
    }
)
setGeneric(
    name="Delta14C_from_AbsoluteFractionModern",
    def=function( # Converts its argument from an Absolute Fraction Modern to a Delta14C representation
    ### The function returns an object of the same type as its input,
    ### which can be of different type.
    ### Have a look at the methods for details.
    AbsoluteFractionModern ##<< A numeric object
    ){
        standardGeneric("Delta14C_from_AbsoluteFractionModern")
    }
)
setGeneric(
    name="AbsoluteFractionModern",
    def=function( # Converts its argument to an Absolute Fraction Modern representation
    ### The function returns an object of the same type as its input,
    ### which can be of different type.
    ### Have a look at the methods for details.
    F ##<< An object that contains data and a format description.  So it can be converted into the AbsoluteFractionModern format if a conversion is implemented.
    ){
        standardGeneric("AbsoluteFractionModern")
    }
)
setGeneric(
    name="AbsoluteFractionModern_from_Delta14C",
    def=function( # Converts its argument to an Absolute Fraction Modern representation
    ### The function returns an object of the same type as its input,
    ### which can be a number or a matrix.
    ### Have a look at the methods for details.
    delta14C){
        standardGeneric("AbsoluteFractionModern_from_Delta14C")
    }
)
setGeneric(
    name="getFormat",
    def=function( # Extracts the format from an object that contains one
    ### The function returns a format string that describes the format of the given data.
    ### More detailed information is given in the methods		 
    object ##<< the type depends on the implementing method
    ){
        standardGeneric("getFormat")
    }
)
setGeneric(
    name="getValues",
    def=function( # Extracts numeric values from an object that contains additional information such as format or units
    ### The function returns the actual  number(s) 
    ### More detailed information is provided in the methods		 
    object ##<< The class of object depends on the method 
    ){
        standardGeneric("getValues")
    }
)
setMethod(
   f= "AbsoluteFractionModern_from_Delta14C",
      signature("numeric"),
      definition=function(# Converts from Delta14C to Absolute Fraction Modern
      ### Converts a number or vector containing Delta14C values to the appropriate Absolute Fraction Modern values.
      ### Have a look at the methods for details.
	delta14C ##<< A numeric object containing the values in Delta14C format
	){
	fprime=(delta14C/1000)+1
	return(fprime)
	}
)
setMethod(
   f= "Delta14C_from_AbsoluteFractionModern",
      signature("numeric"),
      definition=function(# Converts to Delta14C format
      ### This method produces Delta14C values from Absolute Fraction Modern
      ### Have a look at the methods for details.
      AbsoluteFractionModern ##<< A numeric object containing the values in Absolute Fraction Modern format
      ){
        D14C=(AbsoluteFractionModern-1)*1000
        return(D14C)
      }
)
setMethod(
   f= "AbsoluteFractionModern_from_Delta14C",
      signature("matrix"),
      definition=function #Converts from Delta14C to Absolute Fraction Modern
      ### This method produces a matrix of Delta14C values from a Matrix or number of  Absolute Fraction Modern
	(
   delta14C ##<< An object of class matrix containing the values in Delta14C format
	){
	fprime=matrix(
	    nrow=nrow(delta14C),
	    ncol=ncol(delta14C),
	    sapply(delta14C,AbsoluteFractionModern_from_Delta14C)
	)
	return(fprime)
	}
)
setMethod(
   f= "Delta14C_from_AbsoluteFractionModern",
      signature("matrix"),
      definition=function(# Converts Absolute Fraction Modern values to Delta14C
      ### This method produces a matrix of  Delta14C values from a matrix of values in Absolute Fraction Modern.
	AbsoluteFractionModern ##<< An object of class matrix containing the values in Absolute Fraction Modern format
	){
	D14C=matrix(
	    nrow=nrow(AbsoluteFractionModern),
	    ncol=ncol(AbsoluteFractionModern),
	    sapply(AbsoluteFractionModern,Delta14C_from_AbsoluteFractionModern)
	)
	return(D14C)
	}
)
setGeneric ( # This function 
   name= "getMeanTransitTime",
   def=function(# Access to the mean transit time
      ### This generic function assembles methods to 
      ### compute mean transit times.
      ### The nature of the results can change considerably 
      ### depending on the arguments of the function.
      ### For an argument of class \code{model} 
      ### it means something different than for an object of class 
      ### \code{DecompOp}
      ### To interpret them correctly refer also to the documentation
      ### of the methods.
      ##details<<
      ## The concept of mean transit time also known as mean 
      ## residence time can be used to describe 
      ## compartment models. 
      ## In particular it is very frequently used to characterize
      ## linear time invariant compartment models in steady state.  \cite{1, 2, 3, 4}. 
      ## It is very important to note that SoilR is \emph{not} limited to those. 
      ## To integrate the concept into the more general context therefor requires some care. 
      ## This starts with the definition.
      ## Assuming a time invariant system in steady state described by the 
      ## (constant) distribution of inputs to the pools and constant decomposition 
      ## and transfer rates one can define  
      ## the mean transit time as the average time a package of carbon 
      ## spends in the system from entry to exit.
      ## From SoilRs general perspective 
      ## this defintion is ambigous with regard to several points.
      ## \enumerate{
      ## \item
      ## It does not take into account that the mean transfer time may change 
      ## itself with time. For time invariant models in steady state this does not 
      ## matter since the mean transit time turns out to be time invariant also but 
      ## for a general model in SoilR 
      ## input fluxes and decomposition coefficients can be time dependent and the 
      ## system as a whole far from steady state.
      ## \item
      ## It does not specify the set of particles contributing to the mean value. 
      ## If the system is forever in a steady state it is possible to think of the 
      ## average transit time of all particles but if the system changes with time 
      ## such a definition would not be to usefull. To be able to compare 
      ## time dependent  models with real measurements  
      ## the set of particles leaving the system at a certain point in time is a more 
      ## natural choice
      ## }
      ## To incorporate the concept of transit times into SoilR we need to 
      ## address these ambiguities. We also would like the new 
      ## definition to agree with the old one in the special but often studied 
      ## case of linear systems in steady state.
      ## We suggest the following Definition: \cr
      ## Given a system described by 
      ## the complete history of inputs \eqn{\mathbf{I}(t)}{\mathbf{I}(t)} 
      ## for \eqn{t\in (t_{start},t_0)}{t\in (t_{start},t_0)} 
      ## to all pools until time \eqn{t_0}{t_0} 
      ## and the cumulative output \eqn{O(t_0)}{O(t_0)} 
      ## of all pools at time \eqn{t_0}{t_0}
      ## the mean transit time \eqn{\bar T_{t_0}}{\bar T_{t_0}} 
      ## \bold{of the system}
      ## \bold{at time} \eqn{t_0}{t_0} 
      ## is the average of the transit times of all particles leaving the system at time \eqn{t_0}{t_0}
      ## Remark:\cr
      ## For a system with several output channels one could define the mean transit time of particles leaving by this specific channel.
      ## Remark:\cr
      ## In future versions of SoilR it will be possible to compute a dynamic, time dependent mean transit time 
      ## for objects of class \code{ Model}
      ## There is also a method that constructs a time invariant mean transit time by creting a time invariant model in steady state from an input flux distribution and a constant decompostion operators.
      ## This emphasizes that different methods for this function really answer different questions.

      ##references<< Manzoni, S., G.G. Katul, and A. Porporato. 2009. Analysis of soil carbon transit times and age distributions using network theories.
      ## Journal of Geophysical Research-Biogeosciences 114, DOI: 10.1029/2009JG001070.
      ##
      ## Thompson, M.~V. and Randerson, J.~T.: Impulse response functions of terrestrial
      ## carbon cycle models: method and application, Global Change Biology, 5,
      ## 371--394, 10.1046/j.1365-2486.1999.00235.x, 1999.
      ##
      ## Bolin, B. and Rodhe, H.: A note on the concepts of age distribution and transit
      ## time in natural reservoirs, Tellus, 25, 58--62, 1973.
      ##
      ## Eriksson, E.: Compartment Models and Reservoir Theory, Annual Review of Ecology
      ## and Systematics, 2, 67--84, 1971.
     
      object,           ##<< a DecompOp Object. 
      inputDistribution ##<< a vector of length equal to the number of pools. The entries are weights, which must sum to 1.
      
	){standardGeneric("getMeanTransitTime")}
)
setGeneric ( # This function 
   name= "getTransitTimeDistributionDensity",
   def=function(# methods for transit time distributions 
      ### According to  \link{getMeanTransitTime} to we define the related density:
      ##details<< Given a system described by
      ## the complete history of  inputs \eqn{\mathbf{I}(t)}{\mathbf{I}(t)} 
      ## for \eqn{t\in (t_{start},t_0)}{t\in (t_{start},t_0)} 
      ## to all pools until time \eqn{t_0}{t_0} 
      ## and
      ## the cumulative output \eqn{O(t_0)}{O(t_0)} 
      ## of all pools at time \eqn{t_0}{t_0}
      ## the transit time density \eqn{\psi_{t_0}(T) }{\psi_{t_0}(T) }
      ## \bold{of the system}
      ## \bold{at time} \eqn{t_0}{t_0} is the probability density 
      ## with respect to \eqn{T}{T} implicitly defined by
      ##\deqn{\bar T_{t_0} = \int_0^{t-t_{start}} \psi_{t_0}(T) T \;dT}
      ##references<< Manzoni, S., G.G. Katul, and A. Porporato. 2009. Analysis of soil carbon transit times and age distributions using network theories.
                     ## Journal of Geophysical Research-Biogeosciences 114, DOI: 10.1029/2009JG001070.
                object, ##<< a protoDecompOp Object 
                inputDistribution, ##<< a vector of length equal to the number of pools. The entries are weights. That means that their sume must be equal to one!
                times ##<< the times for which the distribution density is sought
	){standardGeneric("getTransitTimeDistributionDensity")}
)
setGeneric (
   name= "getTimes",
   def=function(# Extracts the times argument
	### This functions extracts the times argument from an argument of class NlModel
   object){standardGeneric("getTimes")}
)
setGeneric (
   name= "getInitialValues",
   def=function(# Extracts the times argument
	### This functions extracts the times argument from an argument of class NlModel
   object){standardGeneric("getInitialValues")}
)

setGeneric ( # This function 
   name= "getOutputFluxes",
   def=function# Computes the output flux from all pools, including the proportion transferred to other pools.
      ### This functions computes the output flux for all pools. Note that for any given pool not all the output of the pool is released from the system because it migtht as well be fed into other pools. If you are interested what a pool releases from the system use the method \code{\link{getReleaseFlux}}, which internally makes use of this method but afterwards substracts all parts of the outputs  that are fed to other pools.
   (
	object ##<< An object of class Model or Model14 created by a call to \code{\link{GeneralModel}} or other model creating functions.
  ,as.closures=F ##<< if set to TRUE instead of a matrix a list of functions will be returned.  
	){standardGeneric("getOutputFluxes")
    ##value<< A matrix with m columns representing the number of pools, and n rows representing the time step as specified by the argument
    ##\code{t} in \code{\link{GeneralModel}} or other model creating function.
    ##details<< This function takes a Model object, which represents a system of ODEs of the form 
    ##\deqn{\frac{d \mathbf{C}(t)}{dt} = \mathbf{I}(t) + \mathbf{A}(t) \mathbf{C}(t)}{dC(t)/dt = I(t) + A(t)C(t)} 
    ##and solves the system for \eqn{\mathbf{C}(t)}{C(t)}. The numerical solver used can be specified in \code{\link{GeneralModel}}.
	  ##seealso<< See examples in \code{\link{GeneralModel}}, \code{\link{GeneralModel_14}}, \code{\link{TwopParallelModel}}, 
    ## \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}}, etc.
   }
   
)
setGeneric ( 
   name= "getC",
   def=function(# Calculates the C content of the pools 
    ### This function computes the carbon content of the pools as function of time.
    ### Have a look at the methods for details.
	object ##<< some model object, the actual class depends on the method used.
  ,as.closures=F ##<< if set to TRUE instead of a matrix a list of functions will be returned.  
	){standardGeneric("getC")
    ##value<< A matrix with m columns representing the number of pools, and n rows representing the times as specified by the argument
    ##\code{t} in \code{\link{GeneralModel}} or another model creating function.
    ##details<< This function takes a Model object, which represents a system of ODEs 
    ##and solves the system for \eqn{\mathbf{C}(t)}{C(t)}. The numerical solver used can be specified in the constructors of the Model classes
    ## e.g. \code{\link{Model}},\code{\link{Model_14}},\code{\link{GeneralModel}}.
	  ##seealso<< See examples in \code{\link{GeneralModel}}, \code{\link{GeneralModel_14}}, \code{\link{TwopParallelModel}}, 
    ## \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}}, etc.
   }
   
)
setGeneric( 
   name= "getParticleMonteCarloSimulator",
   def=function# creates an MCSimulator that can be run afterwards to get statistical estimates for mean transit times and ages either pool specific or for the entire system. 
	(object ##<< An object of class NlModel ,NlModel14 or their subclasses created by a call to \code{\link{NlModel}} or other model creating functions.
	){standardGeneric("getParticleMonteCarloSimulator")
	  ##value<< MCS a Monte Carlos Simulator object (basically a closure) that can be run with various input parameters
    }
)
setGeneric ( # This function 
   name= "getReleaseFlux",
   def=function# Calculates the release of C from each pool
   ### This function computes carbon release from each pool of the given model as funtion of time 
   ### Have a look at the methods for details.
   (
	object ##<< An model object (the actual class depends on the method e.g. Model or  Model14 
	    ##details<< This function takes a Model object, which represents a system of ODEs 
	    ## solves the system for \eqn{\mathbf{C}(t)}{C(t)}, calculates a diagonal matrix of release coefficients \eqn{\mathbf{R}(t)}{R(t)}, 
      ## and computes the release flux as \eqn{\mathbf{R}(t) \mathbf{C}(t)}{R(t) C(t)}.
      ## The numerical solver used can be specified in the model creating functions like e.g. \code{\link{Model}}.
  
	){standardGeneric("getReleaseFlux")
    ##value<< A matrix. Every column represents a pool and every row a point in time
	  ##seealso<< See examples in  \code{\link{Model}}, \code{\link{GeneralModel}}, \code{\link{GeneralModel_14}}, \code{\link{TwopParallelModel}}, 
	  ## \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}}, etc.
   }
)
setGeneric ( 
   name= "getAccumulatedRelease",
   def=function# Calculates the accumulated carbon release from the pools as a function of time
   ### This function computes the accumulated carbon release of the given model as funtion of time. 
	(object ##<< A Model object (e.g. of class Model or Model14)  
    ## Have a look at the methods for details.
	){standardGeneric("getAccumulatedRelease")
	  ##value<< A n x m matrix of cummulative release fluxes with m columns representing the number of pools, and n rows representing the times as specified by the argument
	  ##\code{t} in \code{\link{GeneralModel}} or other model creating functions.
	  ##details<< This function takes a Model object, calculates the release flux as specified by \code{\link{getReleaseFlux}}, 
    ##and integrates numerically the release flux up to each time in \code{t}.
	  ##seealso<< See examples in \code{\link{Model}}, \code{\link{GeneralModel}}, \code{\link{GeneralModel_14}}, \code{\link{TwopParallelModel}}, 
	  ## \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}}, etc.
	 }
)
setGeneric ( # This function 
   name= "getC14",
   def=function(# Calculates the mass of radiocarbon (14C fraction times C stock) in all pools
    ### This function computes the mass of 14C (14C fraction times C stock) for all pools as function of time.
    ### Have a look at the methods for details.
	object
	){standardGeneric("getC14")}
)
setGeneric ( # Computes the  \eqn{\frac{^{14}C}{C}}{14C/C} ratio 
   name= "getF14",
   def=function(# Calculates the 14C fraction of all pools
   ### This function computes the radiocarbon fraction for all pools as funtion of time.
   ### Have a look at the methods for details.
	object
	){standardGeneric("getF14")}
)
setGeneric ( # This function 
   name= "getReleaseFlux14",
   def=function(# Calculates the mass of radiocarbon in the release flux (14C fraction times release flux)
   ### This function computes the mass of radiocarbon in the release flux (14C fraction times release flux) as a function of time.
   ### Have a look at the methods for details.
	object
	){standardGeneric("getReleaseFlux14")}
)
setGeneric ( # This function 
  name= "getF14R",
  def=function(# Calculates the average radiocarbon fraction weighted by the amount of carbon release
    ### This function calculates the average radiocarbon fraction weighted by the amount of carbon release at each time step.
    ### Have a look at the methods for details.
    object
    ){standardGeneric("getF14R")}
  )
setGeneric ( # This function 
  name= "getF14C",
  def=function(# Calculates the average radiocarbon fraction weighted by the mass of carbon
    ### This function calculates the average radiocarbon fraction weighted by the mass of carbon at each time step 
    ### Have a look at the methods for details.
    object
    ){standardGeneric("getF14C")}
  )
setGeneric(
    name="getTimeRange",
    def=function(object){
    ### This function returns the time range of the given object. 
    ### Have a look at the methods for details.
        standardGeneric("getTimeRange")
    }
)
setGeneric(
    name="getFunctionDefinition",
    def=function(object){
    ### Extracts the function definition (the R-function) from the argument
    ### Have a look at the methods for details.
        standardGeneric("getFunctionDefinition")
    }
)
setGeneric(
    name="getNumberOfPools",
    def=function(object){
    ### gives the number of poosl from the argument
        standardGeneric("getNumberOfPools")
    }
)
setGeneric(
    name="getOutputReceivers",
    def=function(object,i){
    ### Extracts the pools that recieve output 
        standardGeneric("getOutputReceivers")
    }
)
setGeneric(
    name="getDecompOp",
    def=function(object){
    ### Extracts the Operator from a model object
        standardGeneric("getDecompOp")
    }
)
setGeneric(
    name="getInFluxes",
    def=function(object){
    ### Extracts the InputFluxes from a model object
        standardGeneric("getInFluxes")
    }
)
setGeneric(
    name="availableParticleProperties",
    def=function(object){
    ### Shows the variables available for every particle in every timestep
        standardGeneric("availableParticleProperties")
    }
)
setGeneric(
    name="availableParticleSets",
    def=function(object){
    ### Shows the particle sets available for computation
        standardGeneric("availableParticleSets")
    }
)
setGeneric(
    name="computeResults",
    def=function(object){
        standardGeneric("computeResults")
    }
)
setGeneric(
    name="getDotOut",
    def=function(object){
        standardGeneric("getDotOut")
    }
)
setGeneric(
    name="getTransferMatrix",
    def=function(object){
        standardGeneric("getTransferMatrix")
    }
)
setGeneric(
    name="getTransferCoefficients",
    def=function(object){
        standardGeneric("getTransferCoefficients")
    }
)
setGeneric(
    name="getTransferCoefficients",
    def=function(object,as.closures=F){
        standardGeneric("getTransferCoefficients")
    }
)
setGeneric(
    name="BoundFc",
    def=function # generic constructor
    ### create a BoundFc object from different sources
    (map,starttime,endtime,lag,format,interpolation)
    {
        standardGeneric("BoundFc")
    }
)
setGeneric(
    name="BoundInFlux",
    def=function # generic constructor
    ### create a BoundInFlux object from different sources
    (
      map,
      starttime,
      endtime,
      lag,
      interpolation
     )
    {
        standardGeneric("BoundInFlux")
    }
)
setGeneric(
    name="DecompOp",
    def=function # Generic constructor and converter 
    ### If the argument is already of a subclass of class DecompOp 
    ### the function returns the unchanged object.
    ### Ohterwise it creates an object of a subclass of DecompOp.
    ### Note that the actual class of the output depends on the argument.
    ### For examples, please look at the methods of this function.
    (object)
    {
    ### create a DecompositonOperator from different sources
        standardGeneric("DecompOp")
    }
)
setGeneric(
    name="InFlux",
    def=function # Generic constructor and converter 
    ### If the argument is already of a subclass of class InFlux
    ### the function returns the unchanged object.
    ### Otherwise, it creates an object of a subclass of InFlux.
    ### Note that the actual class of the output depends on the argument.
    ### For examples, please look at the methods of this function.
    (object)
    {
    ### Creates a DecompositonOperator object from different sources
        standardGeneric("InFlux")
    }
)
setGeneric(
    name="ConstLinDecompOp",
    def=function # Generic constructor
    ### Creates a ConstantDecompositonOperator object from different sources.
    ### Please look at the different methods to see what kind of input is supported. 
    (mat)
    {
        standardGeneric("ConstLinDecompOp")
    }
)
setGeneric(
    name="BoundLinDecompOp",
    def=function # Generic constructor
    ### Creates a LinearDecompositonOperator from different sources.
    ### Please look at the methods to see what kind of input is supported. 
    (map,starttime,endtime,lag)
    {
        standardGeneric("BoundLinDecompOp")
    }
)
setGeneric(
    name="GeneralModel",
    def=function # A general constructor 
    ### Creates a Model object from different sources
    ### Have a look at the methods for details.
    (t,A,ivList,inputFluxes,...){
        standardGeneric("GeneralModel")
    }
)
setGeneric(
    name="GeneralModel_14",
    def=function # A general constructor 
    ### Creates a Model14 object from different sources
    ### Have a look at the methods for details.
    (
      t,	
      A,	
      ivList,
      initialValF, 
      inputFluxes, 
      inputFc,
      Fc,
      di=-0.0001209681, 
      lambda=-0.0001209681,
      solverfunc=deSolve.lsoda.wrapper,		##<< The function used by to actually solve the ODE system. This can be \code{\link{deSolve.lsoda.wrapper}} or any other user provided function with the same interface. 
      pass=FALSE,  ##<< if TRUE Forces the constructor to create the model even if it is invalid 
      ...
    )
    {
        standardGeneric("GeneralModel_14")
    }
)
setGeneric(
    name="Model",
    def=function # A general constructor 
    ### Creates a Model object from different sources
    ### Have a look at the methods for details.
    (t,A,ivList,inputFluxes,...){
        standardGeneric("Model")
    }
)
setGeneric(
    name="Model_14",
    def=function # A general constructor 
    ### Creates a Model_14 object from different sources
    ### Have a look at the methods for details.
    (
      t, 
      A,
      ivList,
      initialValF,
      inputFluxes,
      inputFc,
      c14DecayRate,
      solverfunc,
      pass
    ){
        standardGeneric("Model_14")
    }
)

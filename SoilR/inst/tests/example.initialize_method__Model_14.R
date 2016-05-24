require(RUnit)
# We present three possible scenarios:
# 1.) create an object from valid input
# 2.) try to build an Model_14 object with unsound parameters and 
#     show the savety net in action.
# 3.) force an unsound model to be created that would be rejected by default

#1.)
	# create the ingredients and assemble them to a Model_14  in the final step 
	t_start=1960
	t_end=2010
	tn=220
	timestep=(t_end-t_start)/tn 
	## the time 
	t=seq(t_start,t_end,timestep) 
	## some Decomposition Operator
	n=3
	At=new(Class="BoundLinDecompOp",
	  t_start,
	  t_end,
	  function(t0){
	        matrix(nrow=n,ncol=n,byrow=TRUE,
	          c(-1,    0.1,    0, 
	             0.5  , -0.4,    0,   
	             0,    0.2,   -0.1)
	        )
	  }
	) 
	 
	c0=c(100, 100, 100)
	
	## Atmospheric C_14
	
	F0=ConstFc(c(0,10,10),"Delta14C")
	
	## constant inputrate
	inputFluxes=new(
	  "TimeMap",
	  t_start,
	  t_end,
	  function(t0){matrix(nrow=n,ncol=1,c(10,10,10))}
	) 
	# we have a dataframe representing the C_14 fraction 
	# note that the time unit is in years and the fraction is given in
	# the Absolute Fraction Modern format.
	# This means that all the other data provided are assumed to have the same value
	# This is especially true for the decay constants to be specified later
	Fc=BoundFc(C14Atm_NH,format="Delta14C")
	# add the C14 decay to the matrix which is done by a diagonal 
	# matrix which does not vary over time
	# we assume a half life th=5730 years
	th=5730
	k=log(0.5)/th #note that k is negative and has the unit y^-1
	
	mod=new("Model_14",t,At,c0,F0,inputFluxes,Fc,k)

#2.) provoke failure by demanding extrapolation to times where 
#    the model ingredienst are not specified (10 years later than the input e.g)
	t_toLong=seq(t_start,t_end+10,timestep) 
	checkException(new("Model_14",t_toLong,At,c0,F0,inputFluxes,Fc,k),
	"initialize must throw an exception because it is asked to build 
	an unvalid model"
	)
#3.) force an unsound model by settign pass to TRUE
	new("Model_14",t_toLong,At,c0,F0,inputFluxes,Fc,k,pass=TRUE)


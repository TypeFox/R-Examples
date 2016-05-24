# set up some variables 
t_start=0 
t_end=10 

# Create an explicit example of a TimeMap
c0=c(0.5, 0.5, 0.5)
inputFluxes=TimeMap.new(
  t_start,
  t_end,
  function(t0){matrix(nrow=n,ncol=1,c(0.0,0,0))}
) 
# now it you only have the inputFluxes
# you can ask it for which time interval it was specified

print(getTimeRange(inputFluxes))


# Construct a less explicit example of a TimeMap  with an object 
# of class # BoundLinDecompOp which is subclass of TimeMap
n=3
At=new("BoundLinDecompOp",
  t_start,
  t_end,
  function(t0){
        matrix(nrow=n,ncol=n,byrow=TRUE,
          c(-0.2,    0,    0, 
             0  , -0.3,    0,   
             0,      0,   -0.4/t0)
        )
  }
) 
print(getTimeRange(At))


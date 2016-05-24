require(RUnit)
# We present three possible scenarios:
# 1.) create an object from valid input
# 2.) try to build an Model object with unsound parameters and 
#     show the savety net in action.
# 3.) force an unsound model to be created that would be rejected by default
# 4.) show some other insensible models being rejected 
#     
#1.) we first create a sensible model
   t_start=0 
   t_end=10 
   tn=50
   timestep=(t_end-t_start)/tn 
   t=seq(t_start,t_end,timestep) 
   A=new("BoundLinDecompOp",
    t_start,
    t_end,
    function(times){matrix(nrow=3,ncol=3,byrow=TRUE,
        c(-1,    0,    0, 
         0.5,   -2,    0,   
           0,    1, -0.5)
   )    
   }
  )  
  I=TimeMap.new(
     t_start,
     t_end,
     function(times){
       matrix(nrow=3,ncol=1,byrow=TRUE,
           c(-1,    0,    0)
       )
     }
  )
  res=new("Model",t,mat=A,c(0,0,0),I)
  print(class(res))
#2.)
# Now we present some examples where the constructor protests
# test.correctnessOfModel.impossibleCoefficients
   t_start=0 
   t_end=10 
   tn=50
   timestep=(t_end-t_start)/tn 
   t=seq(t_start,t_end,timestep) 

   A=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=3,byrow=TRUE,
            c(-1,    0,    0, 
            1, -0.7,    0,   
            0,    1, -0.5)
        )
      }
   )   
   I=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=1,byrow=TRUE,
            c(-1,    0,    0)
        )
      }
    )
   checkException(
	new("Model",t,A,c(0,0,0),I), 
	"correctnessOfModel should have returned FALSE 
	because the matrix values indicate unbiological 
	behavior (ruwsum should be smaller than zero), 
	but has not"
   )	
#3.) 
# force it nevertheless 
	new("Model",t,A,c(0,0,0),I,pass=TRUE) 

#4.) further examples	
# test.correctnessOfModel.impossibleTimeRanges
   mess="correctnessOfModel should have returned FALSE, but has not"
   t_start=0 
   t_end=10 
   tdiff=t_end-t_start
   tn=50
   timestep=(tdiff)/tn 
   t=seq(t_start,t_end,timestep) 

   #we create an A(t) with sensible coeficients 
   #but where the time range begins to late 

   A=TimeMap.new(
      t_start+1/4*tdiff,
      t_end,
      function(times){
        matrix(nrow=3,ncol=3,byrow=TRUE,
            c(-1,    0,    0, 
            1, -0.7,    0,   
            0,    0.5, -0.5)
        )
      }
   )   
   I=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=1,byrow=TRUE,
            c(-1,    0,    0)
        )
      }
    )
   
   checkException(new("Model",t,A,c(0,0,0),I),mess)
   #now we do the same to the InFluxes(t) while A(t) is correct 
   A=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=3,byrow=TRUE,
            c(-1,    0,    0, 
            1, -0.7,    0,   
            0,  0.5, -0.5)
        )
      }
   )   
   I=TimeMap.new(
      t_start+1/4*tdiff,
      t_end,
      function(times){
        matrix(nrow=3,ncol=1,byrow=TRUE,
            c(-1,    0,    0)
        )
      }
    )
   checkException(new("Model",t,A,c(0,0,0),I),mess)

   #we create an A(t) with sensible coeficients 
   #but where the time range ends to early 

   A=TimeMap.new(
      t_start,
      t_end-1/4*tdiff,
      function(times){
        matrix(nrow=3,ncol=3,byrow=TRUE,
            c(-1,    0,    0, 
            1, -0.7,    0,   
            0,    0.5, -0.5)
        )
      }
   )   
   I=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=1,byrow=TRUE,
            c(-1,    0,    0)
        )
      }
    )
   checkException(new("Model",t,A,c(0,0,0),I),mess)
   #now we do the same to the InFluxes(t) while A(t) is correct 
   A=TimeMap.new(
      t_start,
      t_end,
      function(times){
        matrix(nrow=3,ncol=3,byrow=TRUE,
            c(-1,    0,    0, 
            1, -0.7,    0,   
            0,  0.5, -0.5)
        )
      }
   )   
   I=TimeMap.new(
      t_start,
      t_end-1/4*tdiff,
      function(times){
        matrix(nrow=3,ncol=1,byrow=TRUE,
            c(-1,    0,    0)
        )
      }
    )
   checkException(new("Model",t,A,c(0,0,0),I),mess)

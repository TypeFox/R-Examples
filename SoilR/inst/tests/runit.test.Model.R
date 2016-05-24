#
# vim:set ff=unix expandtab ts=2 sw=2:
require(RUnit)
test.OnepModel=function(){
   attr(OnepModel,"ex")()
}
test.TwopParallelModel=function(){
   attr(TwopParallelModel,"ex")()
}
test.TwopSeriesModel=function(){
   attr(TwopSeriesModel,"ex")()
}
test.TwopFeedbackModel=function(){
   attr(TwopFeedbackModel,"ex")()
}
test.ThreepParallelModel=function(){
   attr(ThreepParallelModel,"ex")()
}
test.ThreepFeedbackModel=function(){
   attr(ThreepFeedbackModel,"ex")()
}
test.ThreepSeriesModel=function(){
   attr(ThreepSeriesModel,"ex")()
}
test.ParallelModel=function(){
  attr(ParallelModel,"ex")()
}
#test.GeneralModel=function(){
#  attr(GeneralModel,"ex")()
#}
test.ICBMModel=function(){
  attr(ICBMModel,"ex")()
}
test.RothCModel=function(){
  attr(RothCModel,"ex")()
}
test.Yasso07Model=function(){
  attr(Yasso07Model,"ex")()
}
test.YassoModel=function(){
  attr(YassoModel,"ex")()
}
test.correctnessOfModel.impossibleCoefficients=function(){
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
   checkException(new("Model",t,A,c(0,0,0),I), "correctnessOfModel should have returned FALSE because the matrix values indicate unbiological behavior (ruwsum should be smaller than zero), but has not")
}
test.correctnessOfModel.impossibleTimeRanges=function(){
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
}
test.correctnessOfModel.correctModel=function(){
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
  res=new("Model",t,A,c(0,0,0),I)
  print(class(res))
  target="Model"
  #checkEquals(target, class(res), "correctnessOfModel should have returned TRUE because the model was correct, but has not")
}

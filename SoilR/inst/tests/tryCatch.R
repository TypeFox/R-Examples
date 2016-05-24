#
# vim:set ff=unix expandtab ts=2 sw=2:
#!/usr/bin/Rscript
source("prolog.R")
mimose=function(){
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
   #correctnessOfModel(t,A,c(0,0,0),I)
   mod=new("Model",t,A,c(0,0,0),I)
   return(mod)
}
   
mimose2=function(x){
    if (x>0) {
        e=simpleError("x larger than 0")
        stop(e)
    } 
    else {y <- x^2 ;return(y)}
}
print(mimose2(-2))
res=tryCatch(
             expr=mimose(),
             error=function(err){100000},#specifying the error handler which has to have the name error
             finally=print('exception handling finished')
             )
print(res)

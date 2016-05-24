gaRoth<- function(
  w_1h = range(SFM_metric[,-1][,1]),
  w_10h = range(SFM_metric[,-1][,2]),
  w_100h = range(SFM_metric[,-1][,3]), 
  w_Live_Herb = range(SFM_metric[,-1][,4]),
  w_Live_Woody = range(SFM_metric[,-1][,5]),
  s_1h = range(SFM_metric[,-1][,6]),
  s_10h = range(SFM_metric[,-1][,7]),
  s_100h = range(SFM_metric[,-1][,8]),
  s_Live_Herb = range(SFM_metric[,-1][,9]), 
  s_Live_Woody = range(SFM_metric[,-1][,10]), 
  delta = range(SFM_metric[,-1][,11]),
  mx.dead = range(SFM_metric[,-1][,12]),
  h_1h = range(SFM_metric[,-1][,13]),
  h_10h = range(SFM_metric[,-1][,14]),
  h_100h = range(SFM_metric[,-1][,15]), 
  h_Live_Herb = range(SFM_metric[,-1][,16]), 
  h_Live_Woody = range(SFM_metric[,-1][,17]),
  m, u, slope, modeltype, obs, 
  method="rmse", maxiter=50, popSize = 20, pcrossover = 0.8, 
  pmutation = 0.1, elitism = base::max(1, round(popSize * 0.05)),
  ...) {
  
SFM_metric<-get (data (SFM_metric,envir = environment ()))
 
ga.min = c(w_1h[1], w_10h[1], w_100h[1], w_Live_Herb[1], w_Live_Woody[1], s_1h[1], s_10h[1], s_100h[1], s_Live_Herb[1], s_Live_Woody[1], delta[1], mx.dead[1], h_1h[1], h_10h[1], h_100h[1], h_Live_Herb[1], h_Live_Woody[1] )

ga.max = c(w_1h[2], w_10h[2], w_100h[2], w_Live_Herb[2], w_Live_Woody[2], s_1h[2], s_10h[2], s_100h[2], s_Live_Herb[2], s_Live_Woody[2], delta[2], mx.dead[2], h_1h[2], h_10h[2], h_100h[2], h_Live_Herb[2], h_Live_Woody[2] )
  
true=obs
if(length(true) ==1) {stop("More than 1 ROS observation is needed")}

fitness <- function (PAR) -error (forecast=unlist(ros(
  w=PAR[1:5],s=PAR[6:10],delta=PAR[11],mx.dead=PAR[12],h=PAR[13:17],
    m=m, u=u, slope=slope,modeltype=modeltype)[15]), 
  true=true, method=method) 


ga (type = "real-valued", 
              fitness=fitness, 
              popSize=popSize, 
              pcrossover=pcrossover, 
              pmutation=pmutation, 
              elitism=elitism, 
              min=ga.min, 
              max=ga.max, 
              maxiter=maxiter,
              names = c("w_1h", "w_10h", "w_100h", "w_Live_Herb", "w_Live_Woody", "s_1h", "s_10h", "s_100h", "s_Live_Herb", "s_Live_Woody", "delta", "mx.dead", "h_1h", "h_10h", "h_100h", "h_Live_Herb", "h_Live_Woody"),
    ...)

}
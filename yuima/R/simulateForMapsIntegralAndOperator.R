# Method for Map
setMethod("simulate", "yuima.Output",
                    function(object, nsim=1, seed=NULL, xinit, true.parameter,
                             space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                             hurst, methodfGn="WoodChan",
                             sampling, subsampling,
                             #Initial = 0, Terminal = 1, n = 100, delta,
                             #	grid, random = FALSE, sdelta=as.numeric(NULL),
                             #	sgrid=as.numeric(NULL), interpolation="none"
                             ...){
                      res <- aux.simulatOutput(object, nsim = nsim, seed = seed,
                        xinit = xinit, true.parameter = true.parameter,
                        space.discretized = space.discretized, increment.W = increment.W,
                        increment.L = increment.L, method = method, hurst = hurst,
                        methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

                      return(res)
                    }
          )

#aux.simulatOutput<-function(yuima.output, param=list(), sampling){
 # param <- list(a=1,b=1,s1=0.1,s2=0.2,a1=0.1,a2=0.1)
aux.simulatOutput<-function(object, nsim, seed, xinit,
  true.parameter, space.discretized, increment.W, increment.L, method, hurst,
  methodfGn, sampling, subsampling){
  mod <- object@model

  if(missing(sampling)){
       sampling <- setSampling()
  }

  sim.Inputs <- simulate(mod, nsim, seed, xinit,
    true.parameter, space.discretized,
    increment.W, increment.L, method, hurst,
    methodfGn, sampling, subsampling)

  for(i in c(1:object@model@equation.number)){
    assign(object@Output@param@Input.var[[i]],
           get.zoo.data(sim.Inputs)[[i]])
  }
  assign(object@Output@param@time.var,
         sim.Inputs@sampling@grid[[1]])
  par <- unlist(true.parameter)
  nam <- names(par)
  for(i in c(1:length(nam))){
    assign(nam[i],par[i])
  }
  my.data<-eval(object@Output@formula[[1]])
  for(i in c(2:prod(object@Output@dimension))){
    my.data<-cbind(my.data,
      eval(object@Output@formula[[i]]))
  }
  names(my.data)<-object@Output@param@out.var

  data1 <- setData(my.data)
  object@data <- data1
  return(object)

}

ModelMaker <-
function(spp,env,n=99,fun, autosim,...){#reconstruct random variables
  if(!is.list(env))env<-list(env=env)
  rownames(spp)<-1:nrow(spp)
  if(identical(names(formals(MAT)),names(formals(fun)))){stop("There is no advantage to using ModelMaker with MAT")}
  if(missing(autosim))rnd<-matrix(runif(nrow(spp)*n),ncol=n)
  else rnd<-autosim
  models<-apply(cbind(as.data.frame(env),rnd),2,function(sim){
         m<-fun(spp, sim,...)
    })
  attr(models,"Nenv")<-length(env)  
  models
}


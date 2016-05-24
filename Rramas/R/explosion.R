explosion <-
function(rmas, bootsp=1000){
  x<-rmas # TODO: change x for rmas in all function
   if(class(x)[1]!="rmas") stop("explosion requires an rmas object
                              (i.e. a trajectory simulation from projectn)")
   if(length(names(x))>0){ # bifurcación provisional para separar rmas de projectn y de projectn2
         x<- x$vn
  }
   
  abundances <- sapply(x, function(rmas) apply(rmas,2,sum))
  abundances.min <-round(apply(abundances[-1,],2,max))
  abminbot <- as.list(1:bootsp)
   # bootstrap maximum abundances
  abminbot <- lapply(abminbot,function(x) x <-sample(abundances.min, replace=TRUE))
  thresholds <- sort(unique(abundances.min))

    # inner function to compute explosion probabilities
    decl.prob <- function(abundances.min, thresholds){
        cf <-NULL
        for(i in 1:length(thresholds)){
          # how many times has ocuured an abundance >= treshold[i]?
          cf <- c(cf, sum(abundances.min >= thresholds[i])) 
        }
        cf <- cf/length(abundances.min)
        cf <- data.frame(Threshold=thresholds, Probability=cf)
        return(cf)
    }

  cf.obs <- decl.prob(abundances.min, thresholds=thresholds)
  cf.boot <- lapply(abminbot, decl.prob, thresholds=thresholds)
  result<- list(cf.obs=cf.obs, cf.boot=cf.boot,abminbot=abminbot,
                 main="Explosion/Increase")
  class(result)<- c("rmas.risk", class(result))
  return(result)
}


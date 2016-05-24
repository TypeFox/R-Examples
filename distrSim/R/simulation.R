### Simulation-Class

#### changed for compatibility with stats 04-10-05 P.R.

### specific generation methods for #runs>1 from version 1.9 on

setMethod("simulate", signature(object = "Simulation"),
  function(object, nsim=-1, seed=-1, ...){
#   function(object, ...){
    if(!is.null(Data(object)))
      return(invisible())
    if(!(seed==-1))
        stop("Seed of an object of class Simulation is changed by seed(<object>,<value>)!")
    if(!(nsim==-1))
        stop("Sample size of an object of class Simulation is changed by samplesize(<object>,<value>)!")
    setRNG(seed(object))
###new:031006:
### later: something like: if(!is(distribution(object), "TSDistribution") {...}
####                       else 'call r(.) with two arguments...'
    Data0<-aperm(array(
             t(r(distribution(object))(object@runs*object@samplesize)),
               c(object@obsDim, object@samplesize, object@runs)),
                 perm=c(2,1,3))
    eval.parent(substitute(object@Data<-Data0))

###old:
### eval.parent(substitute(object@Data<-matrix(
###             r(distribution(object))(object@runs*object@samplesize),
###                        nrow=object@runs,ncol=object@samplesize)))
      return(invisible())
    })


### Cont-Simulation-Class

setMethod("simulate",signature(object="Contsimulation"),
          function(object, nsim=-1, seed=-1, ...){
            if(!is.null(Data(object)))
              return(invisible())
            if(!(seed==-1))
                stop("Seed of an object of class Simulation is changed by seed(<object>,<value>)!")
            if(!(nsim==-1))
                stop("Sample size of an object of class Simulation is changed by samplesize(<object>,<value>)!")

            setRNG(seed(object))

###new:031006:
            data.raw <- r(distribution.id(object))(object@runs*object@samplesize)

            Data.id.raw <- aperm(array( t(data.raw),
                                        c(object@obsDim, object@samplesize,
                                          object@runs)
                                   ), perm = c(2,1,3))

            data.raw <- r(distribution.c(object))(object@runs*object@samplesize)
            Data.c.raw <- aperm(array( t(data.raw),
                                   c(object@obsDim, object@samplesize,
                                     object@runs)
                                  ), perm = c(2,1,3))
            Ind.raw <- matrix(rbinom(object@runs*object@samplesize,1,
                                     object@rate),
                               object@samplesize,object@runs)

            Indx <- array(Ind.raw,c(samplesize(object),runs(object),
                                    obsDim(object))
                         )
            x.id <- aperm(aperm(Data.id.raw, perm = c(1,3,2))*(1-Indx),
                          perm = c(1,3,2))
            x.c  <- aperm(aperm(Data.c.raw,  perm = c(1,3,2))* Indx,
                          perm = c(1,3,2))

            Data.raw <- x.id + x.c
            eval.parent(substitute(object@Data.id <- Data.id.raw))
            eval.parent(substitute(object@Data.c <- Data.c.raw))
            eval.parent(substitute(object@ind <- Ind.raw))
            eval.parent(substitute(object@Data <- Data.raw))
### old:
###            eval.parent(substitute(object@Data.id <-
###                 matrix(r(distribution.id(object))(
###                          object@runs*object@samplesize),
###                          object@runs,object@samplesize)))
###            eval.parent(substitute(object@Data.c <-
###                 matrix(r(distribution.c(object))(
###                          object@runs*object@samplesize),
###                          object@runs,object@samplesize)))
###            eval.parent(substitute(object@ind <-
###                 matrix(rbinom(object@runs*object@samplesize,1,object@rate),
###                               object@runs,object@samplesize)))
###            eval.parent(substitute(object@Data <- (1-object@ind) *
###                       object@Data.id+object@ind*object@Data.c))

            return(invisible())
          })


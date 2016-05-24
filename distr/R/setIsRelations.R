
setIs("Exp", "Gammad", 
       coerce = function(obj){new("Gammad", shape = 1, scale = 1/rate(obj))},
       replace = function(obj, value) 
                     {new("Gammad", shape = value@shape, scale = value@scale)}
     )
   ## a Gamma distribution with shape = 1 and scale = 1/rate(obj)

setIs("Exp", "Weibull", 
       coerce = function(obj) {new("Weibull", shape = 1, scale = 1/rate(obj))},
       replace = function(obj, value) 
                     {new("Weibull", shape = value@shape, 
                                     scale = value@scale)}
     )
   ## a Weibull distribution with shape = 1 and scale = 1/rate(obj)

## if ncp == 0 a Gamma distribution with shape = df(obj)/2 and scale = 2
setIs("Chisq", "Gammad", test = function(obj) isTRUE(all.equal(ncp(obj), 0)),
       coerce = function(obj) {new("Gammad", shape = df(obj)/2, scale = 2)},
      replace = function(obj, value) 
                    {new("Gammad", shape = value@shape, scale = value@scale)}
      )


## if location==0 and scale==1 a T-distribution with df = 1
setIs("Cauchy", "Td", test = function(obj) 
                 {isTRUE(all.equal(location(obj),0)) && 
                  isTRUE(all.equal(scale(obj),1))}, 
       coerce = function(obj) {new("Td")},
       replace = function(obj, value)
                     {new("Td", df = value@df, ncp = value@ncp)}
       ) 


## if Min==0 and Max==1 a Beta Distribution with Param's shape1 = 1, shape2 = 2
setIs("Unif", "Beta", test = function(obj) 
             {isTRUE(all.equal(Min(obj),0)) && 
              isTRUE(all.equal(Max(obj),1))}, 
       coerce = function(obj) {new("Beta", shape1 = 1, shape2 = 1)},
       replace = function(obj, value) {new("Beta", shape1 = value@shape1, 
                          shape2 = value@shape2, ncp = value@ncp)}
       ) 

## if support is affine linear, a DiscreteDistribution is a LatticeDistribution
setAs("DiscreteDistribution", "LatticeDistribution",
      function(from){
        if(!.is.vector.lattice(from@support))
            return(from)
        else{ to <- new("LatticeDistribution")
              slotNames <- slotNames(from)
              lst <- sapply(slotNames, function(x) slot(from,x))
              names(lst) <- slotNames
              lst$lattice <- .make.lattice.es.vector(from@support)
              for (i in 1: length(lst))
                   slot(to, name = names(lst)[i]) <- lst[[i]]
              return(to)}
      })

## if support is affine linear, a DiscreteDistribution is a LatticeDistribution
setAs("AffLinDiscreteDistribution", "LatticeDistribution",
      function(from){
        if(!.is.vector.lattice(from@support))
            return(from)
        else{ to <- new("AffLinLatticeDistribution")
              slotNames <- slotNames(from)
              lst <- sapply(slotNames, function(x) slot(from,x))
              names(lst) <- slotNames
              lst$lattice <- .make.lattice.es.vector(from@support)
              for (i in 1: length(lst))
                   slot(to, name = names(lst)[i]) <- lst[[i]]
              return(to)}
      })

#setIs("DiscreteDistribution", "LatticeDistribution",
#      test = function(object) .is.vector.lattice(support(object)),
#      coerce = function(from) 
#               LatticeDistribution(DiscreteDistribution = from),
#      replace = function(from,value)
#                LatticeDistribution( 
#                r = value@r, d = value@d, q = value@q, p = value@p,
#                support = value@support, img =value@img, 
#                .withSim = value@.withSim, .withArith = value@.withArith,
#                lattice = value@lattice)
#      )


#setAs("LatticeDistribution", "DiscreteDistribution", 
#       function(from) 
#            new("DiscreteDistribution",  r = from@r, d = from@d, q = from@q, p = from@p,
#                support = from@support, img =from@img, 
#                .withSim = from@.withSim, .withArith = from@.withArith)
#      )



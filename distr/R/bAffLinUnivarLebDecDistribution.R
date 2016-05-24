############################## Accessor / Replacement functions

############################## Arithmetics

setMethod("*", c("AffLinUnivarLebDecDistribution","numeric"),
          function(e1, e2) {

          if (length(e2)>1) stop("length of operator must be 1")

          if (isTRUE(all.equal(e2,1))) return(e1)
          if (isTRUE(all.equal(e2,0)))
               return(new("Dirac", location = 0))
          
          if(.isEqual(e1@a*e2,1)&&.isEqual(e1@b,0)){
             obj <- e1@X0
             if(getdistrOption("simplifyD"))
                obj <- simplifyD(obj)
             return(obj)
          }   

          Distr <- UnivarLebDecDistribution(
                     discretePart = discretePart(e1)*e2,
                     acPart = acPart(e1)*e2,
                     discreteWeight = discreteWeight(e1),
                     acWeight = acWeight(e1))

          if(.isEqual(e1@a*e2,1)&&.isEqual(e1@b,0)){
             obj <- e1@X0
             if(getdistrOption("simplifyD"))
                obj <- simplifyD(obj)
             return(obj)
          }   

          Symmetry <- NoSymmetry()
          if(is(e1@Symmetry,"SphericalSymmetry"))
             Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry) + e2)   

          object <- new("AffLinUnivarLebDecDistribution",
                    r = Distr@r, d = Distr@d, p = Distr@p,
                    q = Distr@q, X0 = e1@X0, mixDistr = Distr@mixDistr,
                    mixCoeff = Distr@mixCoeff,
                    a = e1@a*e2, b = e1@b, .withSim  = e1@.withSim,
                    .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1),
                     gaps = gaps(Distr), support = support(Distr),
                     Symmetry = Symmetry
                     )
          object})

setMethod("+", c("AffLinUnivarLebDecDistribution","numeric"),
          function(e1, e2) {
          if (length(e2)>1) stop("length of operator must be 1")
          if (isTRUE(all.equal(e2,0))) return(e1)

          if(.isEqual(e1@a,1)&&.isEqual(e1@b+e2,0)){
             obj <- e1@X0
             if(getdistrOption("simplifyD"))
                obj <- simplifyD(obj)
             return(obj)
          }   
          
          Distr <- UnivarLebDecDistribution(
                     discretePart = discretePart(e1)+e2,
                     acPart = acPart(e1)+e2,
                     discreteWeight = discreteWeight(e1),
                     acWeight = acWeight(e1))

          if(.isEqual(e1@a,1)&&.isEqual(e1@b+e2,0)){
             obj <- e1@X0
             if(getdistrOption("simplifyD"))
                obj <- simplifyD(obj)
             return(obj)
          }   
          
          Symmetry <- NoSymmetry()
          if(is(e1@Symmetry,"SphericalSymmetry"))
             Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry) * e2)   

          object <- new("AffLinUnivarLebDecDistribution",
                    r = Distr@r, d = Distr@d, p = Distr@p,
                    q = Distr@q, X0 = e1@X0, mixDistr = Distr@mixDistr,
                    mixCoeff = Distr@mixCoeff,
                    a = e1@a, b = e1@b+e2, .withSim  = e1@.withSim,
                    .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1),
                     gaps = gaps(Distr), support = support(Distr),
                     Symmetry = Symmetry
                     )
          object})



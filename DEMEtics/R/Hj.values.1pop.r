# Function used within the function calc.r

  Calculate.Hj.values.1pop <- function(allelefrequency3.1poptable){

        Hj.one.population<-Hj(as.numeric(as.vector(allelefrequency3.1poptable$proportion)))
        # The Hj-value for the actual locus and population is calculated.
        Hj.one.population
  }
  

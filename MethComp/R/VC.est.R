VC.est <-
function( data,
           IxR = has.repl(data), linked = IxR,
           MxI = has.repl(data), matrix = MxI,
        corMxI = FALSE, # matrix effects are correlated within items
        varMxI = TRUE,  # variance of matrix effect varies across methods
          bias = TRUE,  # Estimate a bias between methods
         print = FALSE, # Print bias and variance?
 random.raters = FALSE, # Random methods/raters
    lmecontrol = lmeControl(msMaxIter=300) # Control options for lme
        )
{
# To prevnt annoying warnings from error in lme
oo <- options( warn=-1 )
on.exit( options(oo) )

# This function is just a simple wrapper that splits up
# the call depending on whether or not random raters
# are assumed
if( random.raters )
  return( MethComp:::VC.est.random( data = data,
                                     IxR = IxR,
                                     MxI = MxI,
                                  varMxI = varMxI,
                                    bias = FALSE,
                                   print = print,
                              lmecontrol = lmecontrol ) )
else
  return( MethComp:::VC.est.fixed ( data = data,
                                     IxR = IxR,
                                     MxI = MxI,
                                  corMxI = corMxI,
                                  varMxI = varMxI,
                                    bias = bias,
                                   print = print,
                              lmecontrol = lmecontrol) )
}

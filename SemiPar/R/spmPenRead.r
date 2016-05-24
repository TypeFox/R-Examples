 ########## S function: spmPenRead ##########

# For extracting information from a 
# penalized spline term.

# Last changed: 06 JUN 2002

spmPenRead <- function(term)
{

   # Extract argument list

   arg.list <- substring(term,3,(nchar(term)-1))

   # Extract variable name and values

   var.name <- break.string(arg.list,",")[1]
   var.val <- eval(parse(text=var.name))

   # Extract the knots.  
   # If no knots are specified, assign the default
 
   out <- arg.search(arg.list,"knots=")
   present <- out$present

   arg <- out$arg

   if (present==FALSE)
      knots <- default.knots(var.val)

   if (present==TRUE)
      knots <- spmArgRead(arg)$val

   # Extract the smoothing parameter.

   out <- arg.search(arg.list,"spar=")
   present <- out$present
   arg <- out$arg

   if (present==FALSE)
      spar <- NULL

   if (present==TRUE)
      spar <- spmArgRead(arg)$val

   # Extract the degrees of freedom (adf).  
 
   out <- arg.search(arg.list,"adf=")
   present <- out$present
   arg <- out$arg

   if (present==FALSE)
      adf <- "miss"

   if (present==TRUE)
      adf <- spmArgRead(arg)$val


   # Extract the spline basis (basis).  
 
   out <- arg.search(arg.list,"basis=")
   present <- out$present
   arg <- out$arg

   if (present==FALSE)
      basis <- "tps" 

   if (present==TRUE)
      basis <- spmArgRead(arg)$val


 # Extract the polynomial degree (degree).  
 
   out <- arg.search(arg.list,"degree=")
   present <- out$present
   arg <- out$arg

   if (present==FALSE)
      if (basis=="trunc.poly")
         degree <- 1
      else
         degree <- 3

   if (present==TRUE)
      degree <- spmArgRead(arg)$val


   return(list(name=var.name,var=var.val,knots=knots,spar=spar,
               adf=adf,degree=degree,basis=basis))
}

########## End of spmPenRead ##########

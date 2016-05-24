########## S function: spmKrigeRead ##########

# For extracting information from a 
# penalized spline term.

# Last changed: 28 FEB 2003 by MPW

  spmKrigeRead <- function(term)
  {
   # Extract argument list

   arg.list <- substring(term,3,(nchar(term)-1))

   # Extract variable name and values

   var.name <- break.string(arg.list,",")[1:2]
   var.val <- cbind(eval(parse(text=var.name[1])),
                    eval(parse(text=var.name[2])))

   # Extract the knots.  
   # If no knots are specified, assign the default
 
   out <- arg.search(arg.list,"knots=")
   present <- out$present
   arg <- out$arg

   if (present==FALSE)
   { 
      out <- arg.search(arg.list,"k=") 
      arg <- out$arg
      if (!is.null(arg)) 
      {
         num.knots <- spmArgRead(out$arg)$val
         knots <- default.knots.2D(var.val[,1],var.val[,2],num.knots)
      }   
      else
      {
         knots <- default.knots.2D(var.val[,1],var.val[,2])
         num.knots <- nrow(knots)
      }
   }

   if (present==TRUE)
   {
     knots <- spmArgRead(arg)$val
     num.knots <- nrow(knots)
   }

   # Extract the degree of the thin plate spline
   # radial basis functions.

   out <- arg.search(arg.list,"degree=")
   present <- out$present

   arg <- out$arg

   if (present==FALSE)
      degree <- 2
   if (present==TRUE)
      degree <- spmArgRead(arg)$val

   # Extract the smoothing parameter

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
   

   # Extract the boundary.

   out <- arg.search(arg.list,"bdry=")
   present <- out$present
   arg <- out$arg


   if (present==FALSE)
      bdry <- NULL

   if (present==TRUE)
      bdry <- spmArgRead(arg)$val   

  return(list(name=var.name,var=var.val,knots=knots,num.knots=num.knots,
               spar=spar,adf=adf,bdry=bdry,degree=degree))
}

########## End of spmKrigeRead ##########



########## S function: spmArgRead ##########

# For extracting the components of an argument
# assignment.


spmArgReadOS <- function(arg.assignment)
{

   out <- break.string(arg.assignment,"=")

   arg.name <- out[1]

   arg.val <- eval(parse(text=out[2]),envir=sys.frame(1))

   # Extract argument assignment information

   return(list(name=arg.name,val=arg.val))
}

########## End of spmArgRead ##########



gnames.fas <-
function(x = NULL)
{
   if(!inherits(x, "fasta")){
   stop("Make sure the data is a fasta object.")}
   if(is.null(x))
   {stop("You have to specify the input data.")}
   result = x[grepl(">", x)]
   result = gsub(">", "", result)
   return(result)
}


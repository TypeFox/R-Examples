as.fasta <-
function(x){
   if(inherits(x,"phy"))
   {result <- ConvFas(x, type = "phy")}
   else {
    if(inherits(x, "nxs")){
	  result <- ConvFas(x, type = "nxs")
	 }
     else{   
           if(any(regexpr(">", x[seq(1, length(x), by = 2)]) < 0)){
           xx <-  2 * which(regexpr(">", x[seq(1, length(x), by = 2)]) < 0)
		        if( length(xx) > 10){ 
		         xx <- xx[1:10]
	             }
			stop(paste("asfasta could not find \">\" in line: \n",
			     paste(xx,  collapse = ", "),"... \n",
	            "and can not convert", x, " to fasta format\n"))
	        }
	      result <- ConvFas(x, "fas")
		 }
	}
   return(result)
}


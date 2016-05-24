petals <- function(plot=TRUE, txt=TRUE) {

    #######  Don't Cheat  #######


    tmpstr <- "	  		  	 	 	   		   	  	  			  			 	   	 		 	  	    	  		   	  				 	 					    	 		 	 		  		   	 	 	   	 	  	  				 	 		  			  				 	 						  				  	 						  				  	 			 		  				  	 						  				  	 				 	  				  	 						  		 		 	 		  	  	 				    	 	   	 	 		 	 		"
tmpstr2 <- c(
             "
 O
   ","  O

O  ","O
 O
  O","O O

O O","O O
 O
O O","O O
O O
O O")
    ans <- eval(parse(text=rawToChar(packBits( unlist(strsplit(tmpstr,''))==' '))))

    resp <- TRUE
    while(resp) {
        roll <- unlist(dice(1,5,plot.it=plot))
        if(txt) {
            cat("\n---\n")
            cat(tmpstr2[roll], sep='\n---\n')
            cat("---\n")
        }
        petals <- ans(roll)
        resp <- readline('How many petals around the rose? ')
        if(nchar(resp)==0) {
            cat("There were", petals, "petals around the rose\n")
            resp <- FALSE
        } else {
            if(as.numeric(resp)==petals) {
                cat("correct, there were", petals,"petals around the rose\n", sep=' ')
            } else {
                cat("No, there were", petals, "petals around the rose\n", sep=' ')
            }
            resp <- TRUE
        }
    }

####### Don't Cheat   ################

}



## The following lines hid the source code from casual inspection in R 2.13
## but from 2.14 on this is no longer likely to work, see the R-help archives
## for a possible alternative.
#.onAttach <- function(...) {
#    petals <- petals
#    attr(petals,'source') <- "Don't Cheat!"
#    assign('petals',petals,'package:TeachingDemos')
#}

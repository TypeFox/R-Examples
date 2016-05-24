rename.fas <-
function(x, names = NULL)
{
   if(!inherits(x, "fasta")){
      stop("Make sure the data is a fasta object.")
	  }
   if(is.null(names)){
      return(x)
   }
   if(!length(gnames.fas(x)) == length(names)){
	    stop(paste("The number of provided names do not match to \n the number of sequences' names. There are",length(names),"newnames \n compared to",length(x),"original sequence names."))
	 }
      dnas <- x[!grepl(">", x)]
   if(!length(dnas)==length(names)){
	    stop(paste("The number of provided names do not match \n to the number of sequences. There are",length(names),"\n names for ",length(dnas),"sequences."))
	 }
	 namelinenumber <- grep(">", x)
	 result = x
	 for(i in namelinenumber) {
	         result[i] <- names[which(namelinenumber == i)]
           }
	 return(result)
}


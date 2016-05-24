

#File to print out box information in a standard format for reading into CARs

#The proposed format is:

#Box number
#var,lowerbound,upperbound

#with the space left blank otherwise if an upper or lower is not restricted
#and also with entirely unrestricted variables left off

###===VARIABLES====

#box.seq is a "box sequence" as returned by sdprim (possibly via sd.start)
#outfile is csv outfile, character, including extension
#whichboxes is vector of integers indicating which boxes in the _sequence_
#should be written out - note, these are not the numbers from the trajectory
#nobound gives the character to use when there is no lower or upper bound
#separate is the character for separating values - default is comma
#commentchar is the comment character for whatever language will be reading
#  in this csv file


csvboxes <- function(box.seq,outfile,whichboxes=1:length(box.seq),nobound="",separate=",", commentchar="//"){

#convert to vector of dim names
  
  cat("\n")

  sink(outfile) #Tells output that would otherwise go to screen to go to 
                #file instead
  
  for (i in whichboxes){
  
    cat("\n")
    cat("#","\n")
    cat("Box::",i,"\n",sep="") 
    
    box <- box.seq[[i]]
    
    dims <-length(box$dimlist$either)   #how many dims are there?
    
    dnames <- colnames(box$box)         
    
    dimsin <- c(1:dims)[box$dimlist$either]   #index of col numbers that have 
                                              #restricted dims
    
    for (j in dimsin){
    
      if(box$dimlist$lower[j]){  #If there is a lower bound 
        lb <- box$box[1,j]
      } else{ 
        lb <- nobound
      }
      
      if(box$dimlist$upper[j]){  #If there is a upper bound 
        ub <- box$box[2,j]
      } else{ 
        ub <- nobound 
      }
      
      cat(paste(dnames[j],lb,ub,sep=separate),"\n") 
    
    }
  
  cat(commentchar,"Supplemental box information:","\n",sep="")
  cat(commentchar,"Mean:",box$y.mean,"\n",sep="")
  cat(commentchar,"Marginal.coverage:",box$marcoverage,"\n",sep="")
  
  }
  
  sink() #End writing things to file
  
  cat("Output written.","\n")
  
}
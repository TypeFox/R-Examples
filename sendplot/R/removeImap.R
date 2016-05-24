#
# removes data mapping for a given figure (as set by Imap)
#

removeImap <- function(Splot,
                       figure,
                       subset=NA,
                       returnVl=TRUE,
                       saveFlag=FALSE,
                       saveName="Splot.RData"
                       ){

  # if removing default
  if(figure=="Default"){

    idx = which(names(Splot) == "Default")
    # if there is no default set cannot remove
    if(length(idx) == 0){
      cat("Note: There is no default set.\n     Cannot remove\n") 
    }else{
      cat("Removing Default\n")
      temp = Splot
      Splot = list()
      Splot = temp[-idx]
      class(Splot) = "Splot"
    }
    
  # if removing figure
  }else{

    # remove all interactivity for a figure
    if(is.na(subset[1])){
      cat(paste("Removing all maps for figure",figure, "\n")) 
      Splot$iMaps[[figure]] = list()
    # remove subsets iMaps for a figure
    }else{
      # check to make sure subset removing is within number of datasets for a figure
      len = length(Splot$iMaps[[figure]])
      notIn = which(subset > len)
      msg=""
      # if there are subsets outside the maximum number of already specified datasets
      # remove from subset list and print a message 
      if(length(notIn) != 0){
        msg = paste("Note: \n     There are only ", len, " Maps that can be removed for figure ", figure, "\n     Cannot remove map[s] ", subset[notIn[1]], sep="")
        if(length(notIn) > 1){
          for(i in 2:length(notIn)){
            msg = paste(msg, ",", subset[notIn[i]],"\n",sep= "")
          }
        }
        msg = paste(msg, "\n")
        cat(msg)
        subset = subset[-notIn]
      }
      subset = sort(subset)
      # after those that are outside limits are removed
      # if there are still data to remove continue
      if(length(subset) != 0){
        # print which sets are being removed for the figure
        msg = paste("Note: Removing map[s] ", subset[1],sep="")
        if(length(subset) > 1){
          for(i in 2:length(subset)){
            msg = paste(msg, ",", subset[i], "\n", sep="")
          }       
        }
        msg = paste(msg, "\n")
        cat(msg)
        # remove 
        temp = Splot$iMaps[[figure]]
        Splot$iMaps[[figure]] = list()
        Splot$iMaps[[figure]] = temp[-subset]
        # rename 
        if(length(Splot$iMaps[[figure]]) != 0){
          nms = paste("MapObj", c(1:length(Splot$iMaps[[figure]])), sep="")
          names(Splot$iMaps[[figure]]) = nms
        }
     
      }else{ # end if still objects to remove 
        cat("Note: no maps to remove based on given subset\n")
      }
    } # end else for if there is a subset of samples
  }# end if removing figure 

  # save and return
  if(saveFlag) save(Splot, file=saveName, compress=TRUE)
  if(returnVl) return(Splot)

}

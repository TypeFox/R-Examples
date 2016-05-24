cpg.length <-
function(indep,numpatients,covariates,chip.id) {
  theerror<-"Number of subjects per variable do not match up, check your variables\n"
covar<-covariates;
if(length(indep) != numpatients) {
      stop(theerror)
    }
if(!is.null(covar)) {
  if(length(indep)!= nrow(data.frame(covar)) | nrow(data.frame(covar))!=numpatients) {
    stop(theerror)
  }
  if(!is.null(chip.id)) {
    if(length(chip.id)!=nrow(data.frame(covar))) {
      stop(theerror)
  }   } 
     }
  
if(!is.null(chip.id)) {
  if(length(indep)!= length(chip.id) | length(chip.id)!=numpatients) {
     stop(theerror)
    }   }

        }

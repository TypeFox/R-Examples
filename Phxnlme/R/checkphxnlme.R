#' @export
checkphxnlme <- function(testchk=FALSE){
  lic.path <- "C:/Program Files (x86)/Pharsight/Phoenix/application/Plugins/DrugModelEffects/Executables/lservrc"
  lic.path2 <- "C:/Program Files/Pharsight/Phoenix/application/Plugins/DrugModelEffects/Executables/lservrc"
  
  if(testchk==FALSE){
    if (file.exists(lic.path)){
      return(1)
    }else{
      
      if(file.exists(lic.path2)){
        stop("Path")
      }else{
        stop("Phoenix NLME license file not found")  
      }      
    }
  }else{
    if (file.exists(lic.path)){
      return(1)
    }
  }   
}
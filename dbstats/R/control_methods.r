
 ########################
 #### control_method ####
 ########################

 ## description: Internal function. Check if the method is one we have defined.
 ##              If not stop the program.
 ##
 ##        Inputs:  - method ('character', containing the 'method' input by the 
 ##                  user). 
 ##                 - class_call (could be dblm, ldblm or ldbglm)
 ##        Outputs: - method (returned the method, if it's a defined method)
 ##
    
control_method<-function(method,class_call){

  # prospective methods for class dblm
  if (class_call=="dblm"){
    m <- c("OCV", "GCV", "AIC", "eff.rank", "BIC",
        "rel.gvar")   
  }

  if (class_call=="dbglm"){
    m <- c("GCV", "AIC", "eff.rank", "BIC",
        "rel.gvar")   
  }
 
  # prospective methods for class ldblm or ldbglm
  if (class_call=="ldblm"||class_call=="ldbglm"){
    m <- c("OCV", "GCV", "AIC", "eff.rank", "BIC",
         "user.h")
  }

  # prospective methods for class dbpls
  if (class_call=="dbplsr"){
    m <- c("OCV", "GCV", "AIC", "ncomp", "BIC")
  }

  # if character parameter: "method" contains only part of the word, finish
  # filling the word, in accordance with the methods m.
  methodaux <- m[pmatch(method[1],m)]
  if (is.na(methodaux) & class_call!="dbglm"){
    warning(gettextf("the method %s is not defined. Will apply the default method 'OCV'",method))
    method<-"OCV"
  }
   if (is.na(methodaux)&&class_call =="dbglm"){
    warning(gettextf("the method %s is not defined. Will apply the default method 'GCV'",method))
    method<-"GCV"
  }  else
    method<-methodaux
  
  return(method)
}

 ########################
 #### control_family ####
 ########################


 ##  Description: internal function. 
 ##               Get the function of the family "family".
 ##               It can be a 'character' or a 'function'
 ##
 ##      Iputs: family  (family parameter, input by the user.)
 ##      Outputs: family (the function of the family)
 ##


control_family<-function(family){


    # CASE 1: family is a character
    if (is.character(family)){
     
     # prospective families
     f <- c("binomial","gaussian","Gamma","inverse.gaussian","poisson","quasi",
         "quasibinomial","quasipoisson")
     
     # if character parameter: "family" contains only part of the word, 
     # finish filling the word, in accordance with the families "f".
     familyaux <- f[pmatch(family,f)]
     if (is.na(familyaux))
       stop(gettextf("the family %s is not defined",family))
     else
       family<-familyaux
     
     # get the function of the family.              
     family <- get(family, mode = "function", envir = parent.frame())
    }
   
    # CASE 2: family is a function
    if (is.function(family))
        family <- family()
   
    # CASE 3: family is not recognized
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    return(family)
}
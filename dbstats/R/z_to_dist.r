#######################
 #### formula_to_zy ####
 #######################


 ##  Description: internal function.
 ##               Get the z and y of the formula.
 ##
 ##      Iputs:   formula  (formula y~z, input by the user.)
 ##               data: optional data frame with z and y
 ##               mf: dblm call
 ##      Outputs: list with the explanatory variables z and the response y.
 ##

z_to_dist <- function(z,metric){
 
   # see if z is a distance object or have a explanatory variables
   # (in this case compute the distance matrix with daisy)
   if (class(z)[1]=="dist"||class(z)[1]=="dissimilarity"||class(z)[1]=="D2") {
     D<-z
     way<-"D2"
   }else{
     # z must be matrix or data.frame
     if (is.numeric(z))
       z<-as.matrix(z)
     if (is.factor(z))
       z<-data.frame(z)

     if (class(z)!="data.frame"&&class(z)!="matrix")
        stop("'z' must be matrix/data.frame or a numeric vector")
     else {

     D<-daisy(z,metric=metric)    # Computing the Distance matrix with the defined metric of daisy function
     way<-"Z"
      }
    }
    
    return (list(D=D,way=way))
}

 #######################
 #### formula_to_xy ####
 #######################


 ##  Description: internal function.
 ##               Get the x and y of the formula.
 ##
 ##      Iputs:   formula  (formula y~x, input by the user.)
 ##               data: optional data frame with x and y
 ##               mf: dblm call
 ##      Outputs: list with the explanatory variables x and the response y.
 ##

formula_to_zy<-function(formula,data,mf,class_mod,metric){
 
  if (missing(data))
    data <- environment(formula)
 ### recover x and y of the formula
  # number of components of formula, data and weights.
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  # eval the model.frame of the formula ( x, y and weights). save the terms.
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt,"intercept")<-0
  
  zini <-data.frame(mf[,2:ncol(mf)])
  names(zini)<-attr(mt,"term.labels")
  
  # the explanatory variables x (eval the terms labels of formula).
  #names(mf)=c("ncases","ncontrols","agregp")
  if (metric!="gower") 
    z<-model.matrix(mt,mf)
  else {
   if (any(attr(mt,"order")>1))
     z<-model.matrix(mt,mf)
   else{ 
     z <-data.frame(mf[,2:ncol(mf)])
     names(z)<-attr(mt,"term.labels")
   }
  }
  #x <- model.frame(mt,mf,data)
  #x <- x[,attr(mt,"term.labels")]

  # the reponse variable
  if (class_mod=="dblm")
    y <- model.response(mf, "numeric")
  
  if (class_mod=="dbglm")
    y <- model.response(mf, "any")
  
  return (list(z=z,y=y,zini=zini))

}
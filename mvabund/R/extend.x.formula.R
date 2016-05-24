################################################################################
## EXTEND.X.FORMULA: a method extend a formula description to it's extended   ##
## form, eg so that a*b is represented a a=b+a:b                              ##
## additionally formula terms that include several variables can be extended  ##
## and a vector that gives information whether a term is an interaction can   ##
## be returned                                                                ##
################################################################################

extend.x.formula <- function(formula, extend.term=TRUE, return.interaction=TRUE){

  terms.foo <- terms(formula)
  # foo.char  <- deparse(formula)
  term.labels <- attr(terms.foo,"term.labels")
  response <- attr(terms.foo, "response")
  intercept <- attr(terms.foo,"intercept")
  # if(response!=1) new.formula <- reformulate( new.labels ) else
  #    new.formula <- reformulate(new.labels, formula[[2]])

  length.terms <- length(term.labels)

   new.labels <- c()
   is.interact <- c()

   k <- 1   # the position in the new formula
   for(j in 1:length.terms){

      term.j <- terms.foo[j]
      
      # Get the number of influence factors for each term
      # If it's more than one, it's an interaction.
      if(ncol(attr(term.j,"factors"))>1) stop("formula has unknown properties")
      
      inflfac.num    <- sum(attr(term.j,"factors")!=0)
      
      if(inflfac.num == 1) {

          lab <- attr(term.j,"term.labels")   # =  term.labels[j]
          num.labs <- NCOL(eval(parse( text =lab )))
      
          if(num.labs == 1 | ! extend.term){
              new.labels[k] <- lab
              is.interact[k] <- FALSE
              k <- k+1
          } else {
            # Split the labels.
            for(i in 1:num.labs){
              new.labels[k] <- paste(lab,"[,",i,"]", sep="")
              is.interact[k] <- FALSE
              k <- k+1
            }
          }
          
      } else if(inflfac.num > 1) {
        new.labels[k] <- attr(term.j,"term.labels")
        is.interact[k] <- TRUE
        k <- k+1
      } # else leave out this part.


   }

   if(response!=1) {
    new.formula <- reformulate( new.labels )
     if(intercept==0)  new.formula <- update(new.formula,  ~ . + 0)
   } else {
      new.formula <- reformulate(new.labels, formula[[2]])
      if(intercept==0)  new.formula <- update(new.formula,  ~ . + 0)
   }
   
   if(return.interaction) {
      return(list(formula = new.formula, is.interaction = is.interact))
   } else  return(new.formula)
   
}


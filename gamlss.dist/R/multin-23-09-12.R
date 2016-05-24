## A set of functions to help with the multinomial logistic model
## MS and RAR
#-------------------------------------------------------------------------------------------
#   this seems OK
MULTIN <- function(type = "3")
  {
  if (type < 2 | type > 6) stop ("type should be 2, 3, 4 or 5") 
  type<-as.character(substitute(type))
  fam <- switch(type,
                  "2" = BI(),
                  "3" = MN3(),
                  "4" = MN4(),
                  "5" = MN5()
                  )
  fam
  }
#------------------------------------------------------------------------------------------
# MS this seems OK
fittedMN <- function(model)
  {
   #if (!is.gamlss(model))  stop(paste("This is not an gamlss object", "\n", ""))
   type <- as.character(length(model$parameters)+1)
   fit <- switch(type,
                  "2" = cbind(fitted(model), 1-fitted(model)),
                  "3" = cbind(fitted(model,"mu")/(1+fitted(model,"mu")+fitted(model,"sigma")), 
                        fitted(model,"sigma")/(1+fitted(model,"mu")+fitted(model,"sigma")), 
                        1/(1+fitted(model,"mu")+fitted(model,"sigma"))),
                  "4" = cbind(fitted(model,"mu")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")), 
                         fitted(model,"sigma")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")),
                         fitted(model,"nu")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")),          
                         1/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu"))),
                  "5" = cbind(fitted(model,"mu")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")+fitted(model, "tau")), 
                        fitted(model,"sigma")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")+fitted(model, "tau")), 
                        fitted(model,"nu")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")+fitted(model, "tau")), 
                        fitted(model,"tau")/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")+fitted(model, "tau")), 
                        1/(1+fitted(model,"mu")+fitted(model,"sigma")
                                  +fitted(model, "nu")+fitted(model, "tau")))
                  )
   fit
  }
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------  
#MULTINOM <- function(formula,  ... )#type=3,
#  {
#  if (type < 3 | type > 6) stop ("type should be 3. 4 or 5") 
#  type<-as.character(substitute(type))

#  mod <-gamlss(formula = formula, sigma.fo = formula, nu.fo = formula , tau.fo = formula, 
#               family = MULTIN(3), ...)
#  mod
#  }  


#------------------------------------------------------------------------------------------
#dMULTIN <- function(type = "3", ...)
#  {
#  if (type < 3 | type > 6) stop ("type should be 3, 4 or 5") 
#  type<-as.character(substitute(type))
#  fam <- switch(type,
#                  "3" = dMN3(...),
#                  "4" = dMN4(...),
#                  "5" = dMN5(...)
#                  )
#  fam
#  }  
  

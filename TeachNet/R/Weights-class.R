# structure of Weights: (alpha),(apha_1,...,alpha_H),(w_1,...,w_H),(Matrix with all w_ih i in 1,...,I, h in 1,...,H)

setClass(Class = "Weights",
         representation = representation(alpha = "numeric", ## Intercept from output layer
                                         alpha_h = "vector", ## Intercept from hidden layer
                                         w_h = "vector",  ## weights from hidden layer to output layer
                                         w_ih = "matrix"))   ## weights from input layer to hidden layer 
# define addition
setMethod("+", signature(e1 = "Weights", e2 = "Weights"), function (e1, e2) return(new("Weights", 
                                                                                       alpha = e1@alpha+e2@alpha, 
                                                                                       alpha_h = e1@alpha_h+e2@alpha_h, 
                                                                                       w_h = e1@w_h+e2@w_h, 
                                                                                       w_ih = e1@w_ih+e2@w_ih)))
# define subtraction
setMethod("-", signature(e1 = "Weights", e2 = "Weights"), function (e1, e2) return(new("Weights", 
                                                                                       alpha = e1@alpha-e2@alpha, 
                                                                                       alpha_h = e1@alpha_h-e2@alpha_h, 
                                                                                       w_h = e1@w_h-e2@w_h, 
                                                                                       w_ih = e1@w_ih-e2@w_ih)))
# define skalar multiplikation
setMethod("*", signature(e1 = "numeric", e2 = "Weights"), function (e1, e2) return(new("Weights", 
                                                                                       alpha = e1*e2@alpha, 
                                                                                       alpha_h = e1*e2@alpha_h, 
                                                                                       w_h = e1*e2@w_h, 
                                                                                       w_ih = e1*e2@w_ih)))
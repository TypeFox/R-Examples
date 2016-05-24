# structure of Weights2: (alpha),(apha_11,...,alpha_1M),(apha_21,...,alpha_2H),
#                       (w_1,...,w_H),(Matrix with all q_mh m in 1,...,M, h in 1,...,H),
#                       (Matrix with all w_im m in 1,...,I, m in 1,...,M)

setClass(Class = "Weights2",
         representation = representation(alpha = "numeric", ## Intercept from output layer
                                         alpha_1m = "vector", ## Intercept from first hidden layer
                                         alpha_2h = "vector", ## Intercept from second hidden layer
                                         w_h = "vector",  ## weights from second hidden layer to output layer
                                         q_mh = "matrix",   ## weights from first hidden layer to second hidden layer
                                         w_im = "matrix"))   ## weights from input layer to first hidden layer 
# define addition
setMethod("+", signature(e1 = "Weights2", e2 = "Weights2"), function (e1, e2) return(new("Weights2", 
                                                                                         alpha = e1@alpha+e2@alpha, 
                                                                                         alpha_1m = e1@alpha_1m+e2@alpha_1m,
                                                                                         alpha_2h = e1@alpha_2h+e2@alpha_2h,
                                                                                         w_h = e1@w_h+e2@w_h,
                                                                                         q_mh = e1@q_mh+e2@q_mh,  
                                                                                         w_im = e1@w_im+e2@w_im)))
# define subtraction
setMethod("-", signature(e1 = "Weights2", e2 = "Weights2"), function (e1, e2) return(new("Weights2", 
                                                                                         alpha = e1@alpha-e2@alpha, 
                                                                                         alpha_1m = e1@alpha_1m-e2@alpha_1m,
                                                                                         alpha_2h = e1@alpha_2h-e2@alpha_2h,
                                                                                         w_h = e1@w_h-e2@w_h,
                                                                                         q_mh = e1@q_mh-e2@q_mh,  
                                                                                         w_im = e1@w_im-e2@w_im)))
# define skalar multiplikation
setMethod("*", signature(e1 = "numeric", e2 = "Weights2"), function (e1, e2) return(new("Weights2", 
                                                                                        alpha = e1*e2@alpha, 
                                                                                        alpha_1m = e1*e2@alpha_1m,
                                                                                        alpha_2h = e1*e2@alpha_2h,
                                                                                        w_h = e1*e2@w_h,
                                                                                        q_mh = e1*e2@q_mh,  
                                                                                        w_im = e1*e2@w_im)))
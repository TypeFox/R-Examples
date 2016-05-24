gini.gb2 <- function(shape1,shape2,shape3){
G <- gb2.gini(shape1,shape2,shape3,tol=1e-08,maxiter=10000,debug=FALSE)$Gini
return(G)
}

gini.b2 <- function(shape2,shape3){
G <- beta(2*shape2,2*shape3-1)/(beta(shape2,shape3)^2)*(2/shape2)
return(G)
}

gini.dag <- function(shape1,shape2){
G <- ((gamma(shape2)*gamma(2*shape2+1/shape1))/(gamma(2*shape2)*gamma(shape2+1/shape1))-1)
return(G)
}

gini.sm <- function(shape1,shape3){
G <- (1-(gamma(shape3)*gamma(2*shape3-1/shape1))/(gamma(2*shape3)*gamma(shape3-1/shape1)))
return(G)
}




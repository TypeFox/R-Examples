# PPrates-Light-function.R

PP.Smith <- function(L, Pmax, alpha){
 PP <- Pmax*alpha*L/sqrt(Pmax^2 +(alpha*L)^2)
 return(PP)
}

PP.Thornley <- function(L, Pmax, alpha, ksi){
 a <- ksi; b <- alpha*L + Pmax; c <- alpha*L*Pmax
 PP <- (1/(2*a))* (b-sqrt(b^2-4*a*c))
 return(PP)
}

PP.Steele <- function(L, Pmax, Lopt){
 PP <-  Pmax*(L/Lopt)*exp(1-L/Lopt)
 return(PP)
}

PP.Eilers <- function(L, Pmax, alpha, Lopt){
 a <- 1/(alpha*Lopt^2); b <- 1/Pmax-2/(alpha*Lopt); c <- 1/alpha
 PP <- L/(a*L^2+b*L+c)
 return(PP)
}






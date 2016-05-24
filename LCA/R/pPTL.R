pPTL <- function(q,alpha,beta,gamma){
        evaluateCDF <- function(x,alpha,beta,gamma){
                out <- 0
                if(x < -2){
                        out <- 0
                }
                if(x<=0 & x>=-2){
                        out <- ((2*alpha*x)+(alpha*(x^2)/2)+(beta*exp(x/beta))+(4*gamma*x)-(gamma*(x^3)/3))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta)))) - ((2*alpha*(-2))+(alpha*((-2)^2)/2)+(beta*exp((-2)/beta))+(4*gamma*(-2))-(gamma*((-2)^3)/3))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta))))
                }
                if(x>0 & x<=2){
                        out <- 1 + ((2*alpha*(-2))+(alpha*((-2)^2)/2)+(beta*exp((-2)/beta))+(4*gamma*(-2))-(gamma*((-2)^3)/3))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta)))) - ((2*alpha*(-x))+(alpha*((-x)^2)/2)+(beta*exp((-x)/beta))+(4*gamma*(-x))-(gamma*((-x)^3)/3))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta))))
                }
                if(x > 2){
                        out <- 1
                }
                max(c(0,out))
        }
        (unlist(lapply(q,evaluateCDF,alpha=alpha,beta=beta,gamma))-evaluateCDF(-2,alpha=alpha,beta=beta,gamma=gamma))/evaluateCDF(2,alpha=alpha,beta=beta,gamma=gamma)
}


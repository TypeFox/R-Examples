dPTL <- function(x,alpha,beta,gamma){
	evaluatePDF <- function(x,alpha,beta,gamma){
                out <- 0
                if(x < -2){
                        out <- 0
                }
                if(x<0 & x>=-2){
                        out <- (exp(x/beta)+(alpha*(2+x))+(gamma*(4-(x^2))))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta))))

                }
                if(x>=0 & x<=2){
                        out <- (exp(-x/beta)+(alpha*(2-x))+(gamma*(4-(x^2))))/((4*alpha)+(32*gamma/3)+((2*beta)*(1-exp(-2/beta))))
                }
                if(x > 2){
                        out <- 0
                }
                max(c(0,out))
        }
        unlist(lapply(x,evaluatePDF,alpha=alpha,beta=beta,gamma))
}


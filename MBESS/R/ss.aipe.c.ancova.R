ss.aipe.c.ancova <- function(error.var.ancova=NULL, error.var.anova=NULL, rho=NULL, c.weights, width, conf.level=.95, 
assurance=NULL, certainty=NULL)
{
if (is.null(error.var.ancova)&is.null(error.var.anova)) stop("Please specify either the ANCOVA error variance, or both the ANOVA error variance and the correlation coefficient")
 
if(!is.null(assurance)& !is.null(certainty))
    {if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value")}
if(!is.null(certainty)) assurance<- certainty 
 
if(is.null(error.var.ancova))
    {if(is.null(error.var.anova)|is.null(rho)) stop("Please specify either the ANCOVA error variance, or both the ANOVA error variance and the correlation coefficient")
    error.var<- error.var.anova*(1-rho^2)
    }
if(!is.null(error.var.ancova))
    {if(!is.null(error.var.anova)|!is.null(rho)) stop("Since you input the ANCOVA error variance, do not input the ANOVA error variance and the correlation coefficient")
    error.var<- error.var.ancova
    }

if(sum(c.weights)!=0) stop("The sum of the coefficients must be zero")
if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")

alpha <- 1-conf.level
J <- length(c.weights)
sigma <- sqrt(error.var) 

n <- (sigma^2 * 4* (qnorm(1-alpha/2))^2* sum(c.weights^2) ) / width^2
tol<- 1e-6
dif<- tol+1

if(is.null(assurance))
    { while (dif > tol)
        {n.p <- n
        n<- (sigma^2 *4* (qt(1-alpha/2, n*J-J-1))^2 * sum(c.weights^2) ) / width^2
        dif <- abs(n-n.p)
        }
    n<- ceiling(n)
    return(n)
    }

if(!is.null(assurance))
    { while(dif > tol)
        {n.p<-n
        n<- ((sigma^2 *4* (qt(1-alpha/2, n*J-J-1))^2 * sum(c.weights^2) ) / width^2) * ( qchisq(assurance, n*J-J-1) / (n*J-J-1) )
        dif <- abs(n-n.p)
        }
    n<- ceiling(n)
    return(n)
    }

}

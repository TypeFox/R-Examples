ancova.random.data <- function(mu.y, mu.x, sigma.y, sigma.x, rho, J, n, randomized=TRUE)
## generate random data for simple ANCOVA, assuming random covariate, ie, x and y follow the multivariate normal distribution
{
if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

{if(length(mu.y)!=J) stop("'mu.y' should be a J by 1 vector that contains the means of the response variable in each group")

if(randomized==FALSE)
    {if (length(mu.x)!=J) stop ("'mu.x' should be a J by 1 vector that contains the means of the covariate in each group")}
if(randomized)
    {if (length(mu.x)!=1) stop ("In randomized design, the population means of the covariate in each group should be the same")
    mu.x<-array(mu.x, J)
    }
#####################################################
if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

x<-matrix(NA, n, J)
y<-matrix(NA, n, J)
cov.matrix<- matrix(c(sigma.y^2, rho*sigma.y*sigma.x, rho*sigma.y*sigma.x, sigma.x^2), 2) 

for(j in 1:J)
    {data<-MASS::mvrnorm(n, mu=c(mu.y[j], mu.x[j]), Sigma=cov.matrix)
    y[,j]<-data[,1]
    x[,j]<-data[,2]
    }

cbind(y,x)
}
}

#' VN Seasonality Test
#'
#' Seasonality test based on the coefficient of variance for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable Variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{VN.seas(x)}
#' @export


VN.seas <- function(variable){

  output <- vector("numeric", length(variable)/4)
  instances <- vector("numeric", length(variable)/4)

  for (i in 1:(length(variable)/4)){

    if (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])) <
          abs(sd(variable)/mean(variable))){

                            instances[i] <- i

                            output[i]<- (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])))

    }
    else{instances[i] <- 0
    output[i]<- 0
    }
    }
 if(sum(instances[instances>0])==0) {return(1)}

    if(length(instances[instances]>0)>0){
    plot(instances[instances>0],output[output>0],
         xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",
         ylim = c(0,2*abs(sd(variable)/mean(variable))),
         col=ifelse(output[output>0]==min(output[output>0]), "red", "black"), pch =ifelse(output[output>0]==min(output[output>0]), 19, 1))

         abline(h=abs(sd(variable)/mean(variable)), col="red",lty=5)
          text(length(instances[instances>0])/2,abs(sd(variable)/mean(variable)),adj=c(0,-.25),"Variable Coefficient of Variance",col='red')

        n<- rep(abs(sd(variable)/mean(variable)),length(instances[instances>0]))

      M<- matrix(c(instances[instances>0], output[output>0],n),
                     nrow=length(instances[instances>0]),  byrow= FALSE)


    colnames(M) <- c("Period","Coefficient of Variance","Variable Coefficient of Variance")
    print(M)
    }

    if(length(instances[instances>0])>0) {M[which.min(M[,2]),1]} else {1}


  }

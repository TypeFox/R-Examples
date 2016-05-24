#' @title Estimation of Regression Coefficients

#' @description An app to explore estimation of coefficients in simple regression.
#' 
#' @rdname RegEstimate
#' @usage RegEstimate(x=1:10)
#' @param x A numerical vector, specifying the fixed set of x-values. 
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) RegEstimate()
#' }
RegEstimate <-
function(x=1:10)  {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #x gives the set of fixed x-values in the model
  
  manipulate(
    a=slider(-1,1,step=0.1,initial=0,label="True Intercept a"),
    b=slider(-3,3,step=0.1,initial=1,label="True Slope b"),
    s=slider(0,4,step=0.1,initial=2,label="Likely Error Size"),
    coefs=checkbox(FALSE,"Print Coefficient Information"),
  {n <- length(x)
    y=a+b*x+rnorm(n,mean=0,sd=s)
   ymin <- min(a+b*x)-4*s
   ymax <- max(a+b*x)+4*s
   mod <- lm(y~x)
   ahat <- round(coef(mod)[1],2)
   bhat <- round(coef(mod)[2],2)
   plot(x,y,col="blue",pch=16,
        ylim=c(ymin,ymax))
   mtext(bquote(hat(a) == .(ahat)), line= 2)
   mtext(bquote(hat(b) == .(bhat)), line= 0.5)
   abline(a,b,col="red")
   abline(coef(mod),col="blue",lty=2)
   
   if (coefs)  {
     print(summary(mod)$coefficients)
    cat("\n \n")
   }
  
  }
    )
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("a","b","s","coefs"))

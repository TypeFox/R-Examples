#' @title Regression Line Too Shallow?

#' @description The regession line is not as steep as the SD Line (line through point of averages, with slope = sd(y)/sd(x)).  The difference
#' is especially noticeable when the scatterplot is the result of a sample from a bivariate normal distribution.  This app explains
#' why we use the regression line to predict y from x, even though the SD line appears to be a better linear summary of the
#' scatterplot.  Can be used as a starting-point for a discussion of "regression to the mean."
#' 
#' @rdname ShallowReg
#' @usage ShallowReg(n=900,rho=0.5)
#' @param n Number of points in the scatterplot.
#' @param rho Target correlation for the scatterplot.  Points are selected from a standardized bivariate normal distribtuion, with
#' correlation rho.  
#' @return Graphical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note Uses \code{manipulate}, available only in RStudio, and \code{mvrnorm} from package \code{MASS}.
#' @examples
#' \dontrun{
#' if (require(manipulate)) ShallowReg()
#' }
ShallowReg <-
function(n=900,rho=0.5)  {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #n is number of points in cloud
  #rho is the target correlation of cloud
  varcovar <- cbind(c(1,rho),c(rho,1))
  rpoints <- MASS::mvrnorm(n=n,mu=c(0,0),Sigma=varcovar)
  x <- rpoints[,1]
  y <- rpoints[,2]
  big <- max(abs(c(x,y)))*1.1
  
  mod <- lm(y~x)
  x.bounds <- min(x)+(0:10)/10*(max(x)-min(x))
  y.means <- numeric(10)
  for (i in 1:10)  {
    y.means[i] <- mean(y[x >= x.bounds[i] & x<=x.bounds[i+1]])
  }

  manipulate(
    
    slice=slider(1,10,step=1,initial=4,label="Slice of Cloud"),
    showslice=checkbox(FALSE,"Show Slice of Cloud"),
    showlines=checkbox(FALSE,"Show SD and Regression Lines"),
    showmeans=checkbox(FALSE,"Show All Means of Slices"),
  {plot(x,y,pch=16,cex=0.4,col=rgb(0,0,1,0.7),
          xlim=c(-big,big),ylim=c(-big,big))
   
  if (showlines)  {
    abline(coef(mod),col="blue") #Regression Line
    #Now for SD line.  This line also passes through
    #(mean(x),mean(y)), but its slope is sd(y)/sd(x)
    #(or - that if correlation is negative).  When the cloud
    #is result of random sampling from bivariate normal
    #distribtution, the SD line appears to describe the
    #cloud better than the regression line does:
    abline(mean(x)-mean(y)/coef(mod)[2],sign(coef(mod)[2])*sd(y)/sd(x),col="red") 
    if(rho >= 0)  {
      legend("topleft", c("Regression Line","SD Line"),
             fill = c("blue", "red"),cex=0.7)
    } else  {
      legend("topright", c("Regression Line","SD Line"),
             fill = c("blue", "red"),cex=0.7)
    }      
  }
  
  if (showslice)  {
    rect(x.bounds[slice],-big,x.bounds[slice+1],big,
         col=rgb(0,1,0,0.2))
    lines(x=c(x.bounds[slice],x.bounds[slice+1]),y=c(y.means[slice],y.means[slice]),lwd=2)
    x.val <- (x.bounds[slice]+x.bounds[slice+1])/2
    points(x.val,y.means[slice],pch=16,cex=1)
  }
   
  if(showmeans) {
    x.vals <- (x.bounds[1:10]+x.bounds[2:11])/2
    points(x.vals,y.means,cex=1,pch=16)
    #Note that even though it is "too shallow"
    #the regression line is the one to use for
    #predicting y from x.
  }
  
    
  }#end of manip body
    
    )#end manipulate  
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("showlines","showslice","slice","showmeans"))

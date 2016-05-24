#' @title A study of Influence

#' @description An app to explore the idea of influence.  Note how the influence of the blue point wanes as the number of
#' points in the central cloud increases, and also wanes as the correlation of the central cloud increases.
#' 
#' @rdname Points2Watch
#' @usage Points2Watch()  
#' @return Graphical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note Requires package \code{manipulate}, available only in RStudio.
#' Uses \code{mvrnorm} from package \code{MASS}.
#' @examples
#' \dontrun{
#' if (require(manipulate)) Points2Watch()
#' }
Points2Watch <-
function()  {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #In this app the red point (placed near middle of cloud)
  #is not influential, even when it is moved to be an oultier.
  #The blue point (outside) is more influential, but its influence
  #wanes as number of points in cloud increases, and as
  #correlation of cloud approaches 1 or -1.
  
  offloc <- 5
  middle <- c(0,0)
  offside <- c(offloc,0)
  
  big <- offloc
  xmax <- big+1
  xmin <- -xmax
  ymax<- big+1
  ymin <- -ymax
  
  #book-keeping:
  nprev <- -1
  rhoprev <- -2
  currentcol <- "red"

  
  manipulate(
    
    n=slider(10,200,step=10,initial=30,label="Number of points in cloud"),
    rho=slider(-1,1,step=0.01,initial=0,label="Target correlation for cloud"),
    ylevel=slider(ymin,ymax,initial=0,label="y-coordinate of selected point"),
    SelectPoint=picker("red","blue"),
    
    {
      
      varcovar <- cbind(c(1,rho),c(rho,1))
     
     #generate the n small black points, make sure they are in range:
     testval <- offloc+1
     if (!identical(n,nprev) || !identical(rho,rhoprev)) {
        while (testval>offloc)  {
          rpoints <- MASS::mvrnorm(n=n,mu=c(0,0),Sigma=varcovar)
          x <- rpoints[,1]
          y <- rpoints[,2]
          testval <- max(abs(c(x,y)))
      }
        if (!identical(n,nprev)) nprev <<- n
        if (!identical(rho,rhoprev)) rhoprev <<- rho
        
     }
     
     ord.mod <- lm(y~x)
     
     #plot:
     plot(rpoints,pch=16,cex=0.5,
          xlim=c(xmin,xmax),
          ylim=c(ymin,ymax),
          xlab="x",ylab="y")
     
     #The x and y axes:
     abline(0,0,lty=2,col="grey")
     lines(x=c(0,0),y=c(ymin,ymax),lty=2,col="grey")
     
     #Regression line for the n small black points:
     abline(coef(ord.mod))
     
     #The special points:
     if (SelectPoint=="red")  {
       if(SelectPoint==currentcol){
       middle <- c(0,ylevel)
       new.data <- rbind(rpoints,middle)
       offside <- c(5,5)
       } 
       else {
         new.data <- rbind(rpoints,middle)
         currentcol <<-SelectPoint
         offside <- c(5,5)
         middle <- c(0,0)
       }
     }
      
      if (SelectPoint=="blue") {
        if(SelectPoint==currentcol) {
          offside <- c(5,ylevel)
          new.data <- rbind(rpoints,offside)
          middle <- c(0,0)
        } else{
          middle <- c(0,0)
          offside <- c(5,5)
          new.data <- rbind(rpoints,offside)
          currentcol <<- SelectPoint
        }
        
      }
     
     points(x=middle[1],y=middle[2],col="red",pch=16,cex=1.5)
     points(x=offside[1],y=offside[2],col="blue",pch=16,cex=1.5)
     
     #add regression line for n points and the selected point:
     x.new <- new.data[,1]
     y.new <- new.data[,2]
     new.mod <- lm(y.new~x.new)
     abline(coef(new.mod),lwd=2,col=SelectPoint)
  
     legend("topleft", c("Cloud","Cloud + Select Point"),
            fill = c("black", SelectPoint), title="Regression Lines")
     
    }#end of manip body
    
    )#end of manip
  
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker", "row","SelectPoint","ylevel"))

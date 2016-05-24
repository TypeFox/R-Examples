#' @title Variability

#' @description An app to investigate how the variance and sample size 
#' affects the shape of a histogram and violin plot generated from normal data. Summary data (minimum, median,
#' mean, maximum, and quartiles) are displayed in the output for each random sample drawn.
#' 
#' @rdname Variability
#' @usage Variability()
#' @return Graphical and numerical output
#' @export
#' @author Rebekah Robinson \email{rebekah_robinson@@georgetowncollege.edu}
#' @note Uses \code{histogram} and \code{bwplot} from the \code{lattice} package.
#' @examples
#' \dontrun{
#' if (require(manipulate)) Variability()
#' }
Variability=function(){
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    n=slider(10,1000,initial=100,label="Sample Size n"),
    stdev=slider(1,5,initial=3,step=0.1,label="Standard deviation"),
{
  data=rnorm(n,mean=20,sd=stdev)
  print(summary(data))
  p1=lattice::histogram(~data,type="density",xlim=c(3,38),ylim=c(0,0.5),main="Histogram of Some Random Data",xlab=paste("Mean=20, SD=",stdev))
  p2=lattice::bwplot(~data,xlim=c(3,38),main="Violin Plot of the Same Data",xlab=paste("Mean=20, SD=",stdev), panel = function(..., box.ratio) {
    lattice::panel.violin(..., col = "bisque", box.ratio = box.ratio)
    lattice::panel.bwplot(..., box.ratio = 0.1)
  }) 
  print(p1,split=c(1,1,2,1),more=TRUE)
  print(p2,split=(c(2,1,2,1)))
}) 
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("stdev"))

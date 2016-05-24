# Construct the optimal continuous one factor completely randomised design - No overdispersion, A Optimality Only

# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @export
#' @importFrom stats optim
findOptimalApproxDesign<-function(means,silent=FALSE){

  # Find the betas
  b0<-mean(log(means))
  bi<-log(means) - b0

  # Calculate the optimal weights based on Poisson Regression model

  sumbi<- sum(bi)/2
  adjbi<-exp(sumbi - bi/2)
  weights<-adjbi/sum(adjbi)

  weights.table<-as.matrix(treatment=1:length(means),weights)
  # Print output if silent is FALSE
  if(silent==FALSE){
    cat("\n\nThe optimal design weights are:\n")
    print(weights.table)
#     cat(paste("The determinant of the information matrix is: ",round(current,digits=5),"\n",sep=""))
#     cat(paste("Progressvec is: ",progressvec,"\n",sep=""))
  }
  list("design"=weights)
}


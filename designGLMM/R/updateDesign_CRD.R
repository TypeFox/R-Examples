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

# Define an update algorithm
#' @export
#' @importFrom stats rmultinom

updateDesign_CRD <- function(des,ntmt,sige=0,means=c(1,1),probs=c(1)) {  # Generate new design
  numchange<-which(rmultinom(1,1,probs)==1)# Randomly generate the number of points to change
  nobj<-length(des)             # Number of points in the design
  newpos<-sample(1:nobj,size=numchange) # Choose a design point ot replace
  newtmt<-sample(1:ntmt,size=numchange) # Choose a new treatment to replace the old
  des[newpos]<-newtmt           # make replacement
  des                           # return new design
}

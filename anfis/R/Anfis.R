#' ANFIS S4 class implementation in R
#'  
#' Represent a concrete S4 class that represents an Adaptive Neuro Fuzzy 
#' Inference System in R, using type 3 Takagi and Sugeno's fuzzy if-then rule 
#' with multiple outputs.
#'
#' @slot premises list with the MembershipFunctions for each input.
#' @slot consequents numeric matrix with nrow= #rules, ncol= #outputs.
#' @slot rules matrix with the connectivity of the membership functions to the 
#'  rules.
#' @slot X input matrix with ncol=#inputs and nrow=#individuals.
#' @slot Y output matrix with ncol=#output and nrow=#individuals. 
#' @slot errors numeric vector with training errors.
#' @slot trainingType character describing the training algorithm used: 
#'  trainHybridJangOffLine, trainHybridOffLine or trainHybridJangOnLine.
#' @slot fitted.values numeric matrix with predicted values for training data X.
#' @slot residuals numeric matrix with residuals values for training data X.
#' @slot call call class object with training call.
#'
#' @section Features: 
#' \enumerate{
#'   \item Membership Functions (MF) flexible framework:
#'    \itemize{
#'      \item  Flexible user-defined membership functions(MF) extensible class.
#'      \item  Independent number of (MF) for each input.   
#'      \item  Different MF types, if required, for each input. 
#'    }
#'   \item Type 3 Takagi and Sugeno's fuzzy if-then rule 
#'   \item Full Rule combinations, e.g. 2 inputs 2 membership functions this 
#'     means that 4 fuzzy rules will be created.
#'   \item Different learning strategies:
#'    \describe{
#'      \item{trainHybridJangOffLine}{Hybrid learning, i.e. Descent 
#'        Gradient for precedents and Least Squares Estimation for consequents.}
#'      \item{trainHybridJangOnLine}{on-line version with hybrid learning.}
#'      \item{trainHybridOffLine}{Adaptive learning coefficient and 
#'        momentum term.}
#'     }  
#'   \item Multiple outputs support, i.e., the same input partition can be used
#'    to predict more than one output variable.
#' }
#' 
#' @section Functions:
#' ANFIS S4 class includes the following functions:
#' \describe{
#'  \item{initialize}{constructor of ANFIS Architecture to generate the rule 
#'    set and consequents} 
#'  \item{show/print}{generic output of the object}
#'  \item{getRules, getPremises, getConsequents, getErrors, getTrainingType}{
#'    return the respective ANFIS slots}
#'  \item{plotMF}{plot MembershipFunctions domain}
#'  \item{plotMFs}{plot all the MembershipFunctions for the input domain}
#'  \item{plot}{plot training error according with training Type}
#'  \item{LSE}{auxiliary function for Least Square Estimation to avoid singular 
#'    matrix system in off-line training}
#'  \item{trainHybridJangOffLine}{Jang's Hybrid off-line training}
#'  \item{trainHybridJangOnLine}{Jang's Hybrid on-line training}
#'  \item{trainHybridOffLine}{Hybrid off-line training with momentum and 
#'    adaptive learning rate}
#'  \item{summary, fitted, fitted.values, coef, coefficients, resid, residuals}{
#'    wrappers for traditional model functions}
#'  }
#'
#' @include MembershipFunction-evaluateMF.R
#' @import parallel
#' @name ANFIS-class
#' @rdname ANFIS-class
#' @exportClass ANFIS
#' @seealso \code{\link{BellMF-class}}, \code{\link{GaussianMF-class}} and 
#'  \code{\link{NormalizedGaussianMF-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
#' @examples
#' ##Set 2 cores using global options for parallel library
#' require("parallel")
#' if(.Platform$OS.type == "windows"){
#'  options(mc.cores=1)
#' }else{
#'  options(mc.cores=2) ##You could use all calling detectCores()  
#' }
#' 
#' ##Example domain for bidimentional sinc(x,y) function
#' x <- seq(-10, 10, length= 11)
#' trainingSet <- trainSet(x,x)
#' Z <- matrix(trainingSet[,"z"],ncol=length(x),nrow=length(x))
#' persp(x,x,Z,theta = 45, phi = 15, expand = 0.8, col = "lightblue",
#'  ticktype="detailed",main="sinc(x)*sinc(y)")
#' 
#' ##Training domain patterns
#' X <- trainingSet[,1:2]
#' Y <- trainingSet[,3,drop=FALSE]
#' 
#' ##Defining the required MembershipFunctions for the ANFIS
#' membershipFunction<-list(
#'  x=c(new(Class="NormalizedGaussianMF",parameters=c(mu=-10,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=-5,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=0,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=5,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=10,sigma=2))),
#'  y=c(new(Class="NormalizedGaussianMF",parameters=c(mu=-10,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=-5,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=0,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=5,sigma=2)),
#'    new(Class="NormalizedGaussianMF",parameters=c(mu=10,sigma=2))))
#' 
#' ##Creating the ANFIS network with 2 inputs and 4 MembershipFunctions in 
#' ##each input
#' anfis3 <- new(Class="ANFIS",X,Y,membershipFunction)
#' anfis3 
#' 
#' ##Check for epsilon-completeness in each input
#' plotMFs(anfis3)
#' 
#' ##Training the ANFIS network. 
#' \donttest{trainOutput <- trainHybridJangOffLine(anfis3, epochs=10)}
#' ##We will use instead an already trained object to reduce example time.
#' data(anfis3)
#' 
#' ##How the training went. You can keep on training as the training error
#' ##is still descending. 
#' plot(anfis3)
#' 
#' ##Test the fit, i. e., how the MembershipFunctions partition the input space 
#' plotMFs(anfis3)
#' 
#' ##Just to see if premises, consequents and errors were updated
#' getPremises(anfis3)[[input=1]][[mf=1]]
#' getConsequents(anfis3)[1:2,]
#' getErrors(anfis3) #Training errors 
#' getTrainingType(anfis3)
#' names(coef(anfis3))
#' ##An alternative to get premises and/or consequents ... 
#' coef(anfis3)$premises[[input=1]][[mf=1]]
#' coef(anfis3)$consequents[1:2,]
#' 
#' ##First five train pattern associated values for the training process
#' fitted(anfis3)[1:5,]
#' resid(anfis3)[1:5,]
#' summary(anfis3)
#' 
#' ##Surface comparison between the original training set and the predicted 
#' ##ANFIS network
#' y <- predict(anfis3,X)
#' z <- matrix(y[,1],ncol=length(x),nrow=length(x))
#' par(mfrow=c(1,2))
#' persp(x,x,Z,theta = 45, phi = 15, expand = 0.8, col = "lightblue",
#'  ticktype="detailed",main="Goal")
#' persp(x,x,z,theta = 45, phi = 15, expand = 0.8, col = "lightblue",
#'  ticktype="detailed",main="Fitted training Patterns", zlim=c(min(Z),max(Z)))
require(parallel)
ANFIS<-setClass(Class="ANFIS",
  slots=list(
    premises ="list",
    rules = "matrix",
    consequents = "matrix",
    X = "matrix",
    Y = "matrix",
    errors = "numeric",
    trainingType ="character",
    fitted.values = "matrix",
    residuals = "matrix",
    call = "call"),
  validity=function(object){
    ## Check input + MF 
    MF <- unlist(lapply(object@premises,length))
    if(length(MF) != ncol(object@X)){
      stop("Discrepancy between inputs and MembershipFunctions")
    }
    ## Check rules dimensions
    if(any(dim(object@rules) != c(prod(MF),length(MF)))){
      stop("Rules do not have the appropiate dimensions")
    }
    ## Check consequents
    if(any(dim(object@consequents) != c(nrow(object@rules)*
      (ncol(object@X)+1),ncol(object@Y)))){
      stop("Consequents do not have the appropiate dimensions")
    }
    return(TRUE)
  }
)

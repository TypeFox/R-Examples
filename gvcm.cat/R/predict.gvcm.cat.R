predict.gvcm.cat <-
function(
object,   
newdata, 
type="link", # "link", "response"
...
)
{
# check input
    if (("gvcm.cat" %in% is(object))==FALSE )
         stop ("object must be a 'gvcm.cat' object. \n")     
    if (("data.frame" %in% is(newdata))== FALSE)
         stop ("newdata must be a dataframe. \n")
    if ( (type %in% c("link", "response"))==FALSE )
         stop ("type must be one out of 'link', 'response'. \n")

# definitions out of given object
    formula <- object$formula
    control <- object$control

# na action
    if (missing(newdata))
        newdata <- environment(formula)
    data <- na.omit(newdata)

# model.matrix
    dsgn <- design(formula[c(1,3)],data)
    x <- dsgn$X
    x.reduced <- x %*% object$x.reduction
    
# standardize??
    if(control$center){
       x <- scale(x, center = object$centering, scale = FALSE)           
    }
    if(control$standardize){
       x <- scale(x, center = FALSE, scale = object$scaling)
    }

# response
    Y <- model.extract(dsgn$m, "response")
    if (is.factor(Y)==TRUE){Y <- as.numeric(Y)-1}
    if (!is.null(dim(Y)[2]) && object$family$family=="binomial") {
        Y <- Y[,1]/(Y[,1]+Y[,2])
        } 
    if (object$family$family=="binomial" && (sum(Y>1) || sum(Y<0))) 
        stop("No binomial response. \n") 
    if (object$family$family=="Gamma" && (sum(Y<=0))) 
        stop("No Gamma-distributed response. \n") 

# type
    if (type == "link") {link <- function(x){x}} else   # predictor level
                        {link <- object$family$linkinv} # response level

# return
return(list (
fit = link(x %*% object$coefficients), 
fit.refitted = link(x.reduced %*% object$coefficients.refitted),
fit.oml = link(x %*% object$coefficients.oml),
na.action = "na.omit"))  

}


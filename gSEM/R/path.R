##' Extract and display an equation of a pairwise path between two variables.
##'
##' Extract the "best" model between any two variables. The model name and the model equation are printed on screen. The model coefficients, as well as the model R object are also returned.
##' @title Extract Path Coefficients
##' @param x object of class "sgSEMp1", which is the return value of function sgSEMp1().
##' @param from character string. Name of the predictor.
##' @param to character string. Name of the response variable.
##' @param round a positive integer. The coefficients are rounded to this decimal place.
##' @return A list of the following items: 1) model: the best fitted model, 2) model.print: a character string of the model equation and 3) coefs: Model coefficients vector.
##' @export
##' @examples
##' ##' ## Load the sample acrylic data set
##' data(acrylic)
##' 
##' ## Run semi-gSEM principle one
##' ans <- sgSEMp1(acrylic, predictor = "IrradTot", response = "YI")
##'
##' ## Extract relations between IrradTot and IAD2
##' cf <- path(ans, from = "IrradTot", to = "IAD2")
##' print(cf)

path <- function(x, from, to, round = 3){
###------------------
### Model Names
###------------------
    ## IMPORTANT: These model names are taken from the main 'sgSEMp1' function.
    ## It needs to be updated if they are changed in 'sgSEMp1'
    modelNames <- c("SL", "Quad", "SQuad", "Exp", "Log", "nls", "CPSLSL")
    
###------------------
### Checking arguments
###------------------
    if(!inherits(x, "sgSEMp1"))
        stop("x is not of class 'sgSEMp1'.")
    vars <- colnames(x$bestModels)
    if(from %in% modelNames)
        stop("'from' variable is not found.")
    if(to %in% modelNames)
        stop("'to' variable is not found.")

###------------------
### Extract coefficients
###------------------
    if(is.na(x$bestModels[from, to])){
        cat("No relation was found from", from, "to", to, "\n")
        return(NA)
    }else{
        relName <- x$bestModels[from, to]
        theModel <- x$allModels[[from, to, relName]]
        cat("Model type:", relName, "\n")
        if(relName == "SL"){
            coefs <- round(theModel$coef, round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from)
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)
        } else if(relName == "Quad"){
            coefs <- round(theModel$coef, round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, cf(coefs[3]), " * ", from, "^2")
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)
        } else if(relName == "SQuad"){
            coefs <- round(theModel$coef, round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, "^2")
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)
        } else if(relName == "Exp"){
            coefs <- round(theModel$coef, round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * exp(", from, ")")
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)
        } else if(relName == "Log"){
            coefs <- round(theModel$coef, round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * log(", from, ")")
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)
        } else if(relName == "nls"){
            coefs <- round(summary(theModel)$coef[,1], round)
            ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * exp(", cf(coefs[3], TRUE), " * ", from, ")")
            cat("Model equation (round = ", round, "):\n", sep = "")
            print(ans.print)            
        } else if(relName == "CPSLSL"){
            ans.print <- NA
        } else {
            ans.print <- NA
        }
        invisible(list(model = theModel, model.print = ans.print, coefs = coefs))
    }
}

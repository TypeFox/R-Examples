#' Criterion-Related Profile Analysis
#'
#' Implements the criterion-related profile analysis described in Davison & Davenport (2002).
#'
#' The \code{cpa} function requires two arguments: criterion and predictors. The function returns the criterion-related 
#' profile analysis described in Davison & Davenport (2002). Missing data are presently handled by specifying 
#' \code{na.action = "na.omit"}, which performs listwise deletion and \code{na.action = "na.fail"}, the default, 
#' which causes the function to fail. The following S3 generic functions are available: \code{summary()},\code{anova()}, 
#' \code{print()}, and \code{plot()}. These functions provide a summary of the analysis (namely, R2 and the level a
#' nd pattern components); perform ANOVA of the R2 for the pattern, the level, and the overall model; provide 
#' output similar to \code{lm()}, and plots the pattern effect.
#' @export
#' @importFrom stats glm coef fitted cor pf
#' @param formula An object of class \code{\link{formula}} of the form \code{response ~ terms}.
#' @param data An optional data frame, list or environment containing the variables in the model.
#' @param k Corresponds to the scalar constant and must be greater than 0. Defaults to 100.
#' @param na.action How should missing data be handled? Function defaults to failing if missing data are present.
#' @param family A description of the error distribution and link function to be used in the model. See \code{\link{family}}.
#' @param weights An option vector of weights to be used in the fitting process.
#' @return An object of class \code{critpat} is returned, listing the following components: 
#' \itemize{
#'  \item \code{lvl.comp} - the level component
#'  \item \code{pat.comp} - the pattern component
#'  \item \code{b} - the unstandardized regression weights
#'  \item \code{bstar} - the mean centered regression weights
#'  \item \code{xc} - the scalar constant times bstar
#'  \item \code{k} - the scale constant
#'  \item \code{Covpc} - the pattern effect
#'  \item \code{Ypred} - the predicted values
#'  \item \code{r2} - the proportion of variability attributed to the different components
#'  \item \code{F.table} - the associated F-statistic table
#'  \item \code{F.statistic} - the F-statistics
#'  \item \code{df} - the df used in the test
#'  \item \code{pvalue} - the p-values for the test}
#' 
#' @examples
#' \dontrun{
#' data(IPMMc)
#' mod <- cpa(R ~ A + H + S + B, data = IPMMc)
#' print(mod)
#' summary(mod)
#' plot(mod)
#' anova(mod)
#' }
#'
#' @references Davison, M., & Davenport, E. (2002). Identifying criterion-related patterns of predictor scores using multiple regression. \emph{Psychological Methods, 7}(4), 468-484. DOI: 10.1037/1082-989X.7.4.468.
#' @seealso \code{\link{pcv}}
#' @keywords method

cpa <- function(formula, data, k=100, na.action = "na.fail", family = "gaussian", weights = NULL){
	if(is.null(weights))
    regweg <- glm(formula=formula,data=data,family = family,na.action = na.action)
    else regweg <- glm(formula=formula,data=data,family = family,na.action = na.action,weights=weights)
	  b <- coef(regweg)[-1]
    bstar <- b - mean(b)
    xc <- k*bstar
    
    if(is.null(weights))
    x <- regweg$model[,-1]
    else x <- regweg$model[,c(-1,-ncol(regweg$model))]
  	y <- regweg$model[,1]
    N <- nrow(x)
    v <- ncol(x)
    V <- 1/ncol(x)
    pat.comp <- x - apply(x,1,mean)
    Xp <- apply(x,1,mean)
    Covpc <- V*(as.matrix(pat.comp)%*%as.matrix(xc))
    ypred <- fitted(lm(y ~ 1 + Covpc + Xp, na.action = na.action))
    R2.f <- cor(ypred,y)^2
    R2.pat <- cor(y,as.vector(Covpc))^2 ## pattern effect
    R2.lvl <- cor(y,Xp)^2 ## level effect
    r2 <- rbind(R2.f,R2.pat,R2.lvl)
    colnames(r2) <- "R2"
    rownames(r2) <- c("Full Model","Pattern","Level")
    r2 <- round(r2,digits=6)
    
    full.df <- c(v,N-v-1)
    F.R2.full <- (R2.f*full.df[2])/((1-R2.f)*full.df[1])
    p.value.F.R2.full <- pf(F.R2.full,full.df[1],full.df[2],lower.tail=FALSE)
    
    pat.df <- c(v-1,N-v-1)
    F.R2.pat <- ((R2.f - R2.lvl)*pat.df[2])/((1-R2.f)*pat.df[1])
    p.value.F.R2.pat <- pf(F.R2.pat,pat.df[1],pat.df[2],lower.tail = FALSE)
    
    lvl.df<-c(1,N-v-1)
    F.R2.lvl <- ((R2.f - R2.pat)*lvl.df[2])/((1-R2.f))
    p.value.F.R2.lvl <- pf(F.R2.lvl,lvl.df[1],lvl.df[2],lower.tail=FALSE)
    
    F.R2.pat.only <- ((R2.pat)*pat.df[2])/((1-R2.pat)*pat.df[1])
    p.value.F.R2.pat.only <- pf(F.R2.pat.only,pat.df[1],pat.df[2],lower.tail = FALSE)
    
    F.R2.lvl.only <- ((R2.lvl)*lvl.df[2])/((1-R2.lvl)*lvl.df[1])
    p.value.F.R2.lvl.only <- pf(F.R2.lvl.only,lvl.df[1],lvl.df[2],lower.tail = FALSE)
    
    fvalue <- c(F.R2.full,F.R2.pat.only,F.R2.lvl.only,F.R2.pat,F.R2.lvl)
    df <- rbind(full.df,pat.df,lvl.df,pat.df,lvl.df)
    pvalue <- rbind(p.value.F.R2.full,p.value.F.R2.pat.only,p.value.F.R2.lvl.only,p.value.F.R2.pat,p.value.F.R2.lvl)
    ftable <- cbind(df,fvalue,pvalue)  
    
    rownames(ftable) <- c("R2.full = 0 ","R2.pat = 0","R2.lvl = 0","R2.full = R2.lvl","R2.full = R2.pat")
    colnames(ftable) <- c("df1", "df2", "F value", "Pr(>F)")
 
    call<- match.call()
    output <- list(call=call,lvl.comp=Xp,pat.comp=pat.comp,b=regweg,bstar=bstar, xc=xc, k=k, Covpc=Covpc, Ypred=ypred,r2=r2,ftable=ftable)
    
    class(output) <- "critpat"
    return(output)
  }
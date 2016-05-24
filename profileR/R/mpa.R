#' Moderated Profile Analysis
#'
#' Implements the moderated profile analysis approach developed by Davison & Davenport (unpublished) - EXPERIMENTAL SUPPORT
#'
#' The function returns the criterion-related moderated profile analysis described in Davison & Davenport (unpublished). Missing data are presently handled by specifying \code{na.action = "na.omit"}, which performs listwise deletion and \code{na.action = "na.fail"}, the default, which causes the function to fail. The following S3 generic functions are not yet available but will be in future implementations. \code{summary()},\code{anova()}, \code{print()}, and \code{plot()}. These functions provide a summary of the analysis (namely, R2 and the level and pattern components); perform ANOVA of the R2 for the pattern, the level, and the overall model; provide output similar to \code{lm()}, and plots the pattern effect.
#' @export
#' @importFrom stats lm anova coef update variable.names
#' @param formula An object of class \code{\link{formula}} of the form \code{response ~ terms}.
#' @param data An optional data frame, list or environment containing the variables in the model.
#' @param moderator Name of the moderator variable.
#' @param k Corresponds to the scalar constant and must be greater than 0. Defaults to 100.
#' @param na.action How should missing data be handled? Function defaults to failing if missing data are present.
#' @param stage2 Should stage 2 be executed regardless of stage 1 outcome? defaults to FALSE.
#' @return An object of class \code{critpat} is returned, listing the following components:
#' \itemize{
#'  \item \code{r2_full} - the proportion of variability attributed to the full model
#'  \item \code{r2_red} - the proportion of variability attributed to the reduced model
#'  \item \code{ftable} - the associated F-statistic table
#'  \item \code{F-statistic} - the F-statistics
#'  \item \code{df} - the dfs used in the test
#'  \item \code{pvalue} - the p-values for the test}
#'
#' @examples
#' \dontrun{
#' data(mod_data)
#' mod <- mpa(dv ~ pred1 * mod + pred2 * mod, moderator = "mod", data = mod_data)
#' }
#'
#' @references Davison, M., & Davenport, E. (unpublished). Comparing Criterion-Related Patterns of Predictor Variables across Populations Using Moderated Regression.
#' @seealso \code{\link{cpa}}
#' @keywords method

mpa <- function(formula, data, moderator, k=100, na.action = "na.fail", stage2 = FALSE){
  cat("# -------- Executing Stage 1 --------  #\n\n")
  stage1_mod <- lm(formula=formula,data=data,na.action = na.action)
  print(stage1_tab <- anova(stage1_mod))
  if(stage1_tab[nrow(stage1_tab) - 1, ncol(stage1_tab)] < .05 | stage2 == TRUE){
    cat("\n# ------- Executing Stage 2 --------  #\n\n")
    
    # Drop response and moderator variables
    loc <- which(names(stage1_mod$model) == moderator)
    resp <- stage1_mod$model[,1]
    
    # Identify number of predictors
    pred_num <- length(names(stage1_mod$model)[-c(1, loc)])
    x <- stage1_mod$model[,-c(1, loc)]
    z <- stage1_mod$model[,loc]
    z <- ifelse(z == levels(z)[1], 1, 0)
    
    # Find level effect
    Xp <- apply(x, 1, mean)
    
    # Pattern component by groups
    dat <- stage1_mod$model
    mod_data <- dat[,loc]
    reg_names <- levels(data[, moderator])
    ref <- which(mod_data == reg_names[1])
    foc <- which(mod_data == reg_names[2])
    x.ref <- dat[ref,-c(1, loc)]
    x.foc <- dat[foc,-c(1, loc)]
    pat.comp <- x - apply(x,1,mean)
    pat.comp.ref <- matrix(ncol = pred_num, rep(0,length(z)*pred_num))
    pat.comp.foc <- matrix(ncol = pred_num, rep(0,length(z)*pred_num))
    rownames(pat.comp.ref) <- c(1:dim(pat.comp.ref)[1])
    rownames(pat.comp.foc) <- c(1:dim(pat.comp.foc)[1])
    pat.comp.ref.t <- x.ref - apply(x.ref,1,mean)
    pat.comp.ref[which(mod_data == reg_names[1]),] <- as.matrix(pat.comp.ref.t)
    pat.comp.foc.t <- x.foc - apply(x.foc,1,mean)
    pat.comp.foc[which(mod_data == reg_names[2]),] <- as.matrix(pat.comp.foc.t)
    
    # Set up the regression weight contrast
    
    bref <- coef(stage1_mod)[-(grep(levels(data$mod)[2], names(coef(stage1_mod))))]
    bref <- bref[-1]
    brefsum <- sum(bref) 
    refstar <- (bref - mean(bref)) * k
    bfoc <- coef(stage1_mod)[-1]
    
    variable_names <- variable.names(stage1_mod$model)[-1]
    variable_names <- variable_names[-(grep(moderator, variable_names))]
    vlen <- length(variable_names)
    bfoc_tmp <- list()
    for(i in 1:vlen){
    bfoc_tmp[[i]] <- sum(coef(stage1_mod)[grep(variable_names[i], names(coef(stage1_mod)))])
    }
    bfoc_tmp <- unlist(bfoc_tmp)
    focstar <- (bfoc_tmp - mean(bfoc_tmp)) * k
    bfocsum <- sum(bfoc_tmp)
    
    V <- length(bref)
    Covpc.ref <- V*(as.matrix(pat.comp)%*%as.matrix(refstar))
    Covpc.foc <- V*(as.matrix(pat.comp)%*%as.matrix(focstar))
    
    model_data <- data.frame(resp, first = brefsum*Xp, sec = (bfocsum - brefsum)*z*Xp, thr = V*Covpc.ref, four =  (Covpc.foc - Covpc.ref)*V*z, z = z)
    
    lm_full <- lm(resp ~ first + sec + thr + four + z, model_data)
    lm_red <- update(lm_full, .~. -four)
    r2_full <- summary(lm_full)$r.squared
    r2_red <- summary(lm_red)$r.squared
    F <- (r2_full - r2_red)/(1 - r2_full) * (nrow(x) - 2*V - 2)/(V - 1)
    p.value <- 1 - pf(F, V- 1, nrow(x) - 2*V - 2)
    ftable <- cbind(F, df1 = V - 1, df2 = nrow(x - 2*V - 2), p.value)
    colnames(ftable) <- c("F-statistic", "df1", "df2", "p - value")
    rownames(ftable) <- c("R2full = R2red")
    cat("# -------- F-table --------  #\n\n")
    cat(paste("The F-statistic is:",F, "with", V-1, "and", nrow(x) - 2*V - 2, "degrees of freedom.\nThe p-value is:", p.value,"\n"))
    cat(paste("The R2 squared for the full model is:", r2_full, "and the R2 squared for the reduced model was:", r2_red,"\n"))
    call<- match.call()
    output <- list(call=call,r2_full, r2_red,ftable=ftable)
    #class(output) <- "critpat"
    #return(output)
   } else {
    cat("Interaction not significant. Not executing stage 2\n")
    break
   }
  }
#' Automated Multiple Regression Modelling
#'
#' @param outcome The dependent variable of the hierarchical model
#' @param block1 A character vector, with names of variables. The first block of
#'   independent variables.
#' @param ...  A character vector, with names of variables. Subsequent blocks of
#'   independent variables.
#' @param dataset A data frame containing variables refered to in
#'   \code{formulas}, passed to data argument of \code{lm}
#' @param type Family argument to pass to \code{glm}.  Specify "binomial" for
#'   binary logistic regression models.
#' @param assumptions.check Boolean, if TRUE, then assumption checks are run and
#'   output is produced. If FALSE, only model summary and coefficient tables are
#'   produced.
#' @param outliers.check Determines how many observations to display for
#'   outliers check.  Default is significant observations. "All" shows all
#'   residual and Cook's D values.
#' @param transform.outcome A boolean. If TRUE, a variable transformation of the
#'   outcome is substituted in the final model if outcome is non-normal. NOT
#'   IMPLEMENTED YET.
#' @details Calls other functions to generate model objects and test them, given
#'   specified model parameters and other options.  Formatted output is produced
#'   via \code{model_output}
#' @examples
#' run_model("y", c("lag.quarterly.revenue"), c("price.index", "income.level"),
#' dataset=freeny)
#' @export
run_model <- function(outcome, block1, ..., dataset, type="gaussian", assumptions.check = T, outliers.check = "significant", transform.outcome=F){
  # Main function that calls the others to generate model objects, and test and summarize those model objects
  if(type == "binomial"){
    forms <- create_formula_objects(outcome, block1, ...)
    models <- create_model_objects(forms, dataset, type="binomial")
    model_output_binomial(models, forms)
  }else{
  forms <- create_formula_objects(outcome, block1, ...)
  full_model <- lm(forms[[length(forms)]], dataset)
  full_model_frame <- model.frame(full_model)
  dataset <- na.omit(full_model_frame)
  models <- create_model_objects(forms, dataset)
  top_model <- models[[length(models)]]
  if(assumptions.check == T) {
    if(outliers.check == "significant") {
    checks <- assumptions_check(top_model)
    model_output(models, dataset, checks, forms, outliers = "significant")
    } else {
      checks <- assumptions_check(top_model)
      model_output(models, dataset, checks, forms, outliers = "all")
    }
  } else {
    model_output(models, formulas=forms)
  }
  }
}

#' Hierarchical Formula Generation
#'
#' @param outcome The dependent variable of the hierarchical model
#' @param block1 A character vector, with names of variables. The first block of
#'   independent variables.
#' @param ...  A character vector, with names of variables. Subsequent blocks of
#'   independent variables.
#' @return A list of \code{lm} formulas
#' @examples
#' create_formula_objects("y", c("lag.quarterly.revenue"), c("price.index"))
#' create_formula_objects("y", c("lag.quarterly.revenue"), c("price.index",
#' "income.level"))
#' @export
create_formula_objects <- function(outcome, block1, ...){
  # Creates formulas for hierarchical models from blocks of predictors
  # Pass character vectors with names of independent variables corresponding to each block
  blocks <- list(...)
  formula <- as.formula(paste(outcome, "~", paste(block1, collapse="+")))
  formulas <- list()
  if(length(blocks) != 0){
    formulas <- list(formula)
  for(block in blocks){
  formulas[[length(formulas) + 1]] <-  as.formula(paste(formulas[length(formulas)], "+", paste(block, collapse="+")))
  }
  }
  if(length(formulas) == 0){
    list(formula)
  } else {
    formulas
  }
}

#' Hierarchical Regression Model Generation
#'
#' @param formulas A set of \code{lm} formulas, created with
#'   create_formula_objects
#' @param dataset A data frame containing variables refered to in
#'   \code{formulas}, passed to data argument of \code{lm}
#' @param type Family argument to pass to \code{glm}.  Specify "binomial" for
#'   binary logistic regression models.
#' @return A list of \code{lm} model objects
#' @examples
#' create_model_objects(create_formula_objects("y", c("lag.quarterly.revenue")
#' , c("price.index")), dataset = freeny)
#' freeny_model_formulas <- create_formula_objects("y", c("lag.quarterly.revenue")
#' , c("price.index"))
#' create_model_objects(freeny_model_formulas, dataset = freeny)
#' @export
create_model_objects <- function(formulas, dataset, type="gaussian"){
  if (type == "binomial") {
    models <- lapply(X = formulas, data=dataset, glm, family="binomial")
  } else {
  # Creates all of the hierarchical models from a set of formulas created with create_formula_objects
  if (length(formulas) < 2) stop("Your model contains just one block and is not hierarchical, consider lm")
  models <- lapply(X = formulas, data=dataset, lm)
  }
}

#' Multiple Regression Assumption Checking
#'
#' @param model A \code{lm} model object.  \code{run_model} automatically calls
#'   this function for the model with all blocks of predictors included.
#' @details Creates objects related to multiple regression assumption checking.
#'   These objects are used by \code{model_output} to produce readable output.
#' @examples
#' freeny_model_formulas <- create_formula_objects("y", c("lag.quarterly.revenue")
#' , c("price.index"))
#' freeny_models <- create_model_objects(freeny_model_formulas,
#' dataset = freeny)
#' freeny_model <- freeny_models[[length(freeny_models)]]
#' assumptions_check(freeny_model)
#' @export
assumptions_check <- function(model){
  # Creates objects needed for assumption checking and output printing
  dw <- lmtest::dwtest(model)$statistic
  dwp <- lmtest::dwtest(model)$p.value
  cormat <- cor(data.frame(lapply(model.frame(model), as.numeric)), use="pairwise.complete.obs")
  vifs <- car::vif(model)
  standresids <- MASS::stdres(model)
  cdists <- cooks.distance(model)
  if(length(standresids) < 5000){
  normresids <- shapiro.test(standresids)$p.value
  } else {
  normresids = NULL
  }
  probDist <- pnorm(MASS::stdres(model))
  results <- list(Durbin.Watson=dw, DW.p.value=dwp, Correlation.Matrix=cormat, Var.Inf.Factor=vifs, Stand.Residuals  =standresids, Cooks.Dist=cdists, Normality.Resids=normresids, Prob.Dist=probDist)
  invisible(results)
}

#' Multiple Regression Output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @param data The dataframe from which the model was built.
#' @param checkList a list object created by \code{assumptions_check} used to
#'   create output.
#' @param formulas Formula list produced by \code{create_formula_objects}, used
#'   for summary table.
#' @param outliers Outlier option, select the number of observations to examine
#'   for outliers.
#' @details Creates plots and text output to summarize models and check
#'   assumptions via objects created by \code{assumptions_check}.  Uses full
#'   model with all predictors.
#' @examples
#' freeny_model_formulas <- create_formula_objects("y",
#' c("lag.quarterly.revenue"), c("price.index"))
#' freeny_models <- create_model_objects(freeny_model_formulas,
#' dataset = freeny)
#' freeny_model <- freeny_models[[length(freeny_models)]]
#' checks <- assumptions_check(freeny_model)
#' model_output(freeny_models, freeny, checks, freeny_model_formulas,
#' outliers = "significant")
#' @export
model_output <- function(models, data, checkList=NULL, formulas, outliers){
  # Produces plots and prints relevant messages and outputs.
  options(scipen = 100, digits = 4)
  if(!is.null(checkList)){
  model <- models[[length(models)]]
  cat("\n\nREGRESSION OUTPUT\n\n")
  cat("Durbin-Watson = ", checkList$Durbin.Watson, "p value = ", checkList$DW.p.value, "\n\n")
  cat("Partial Regression plots (all relationships should be linear):\n\n")
  car::avPlots(model)
  cat("Plot of studentized residuals: uniform distibution across predicted values required")
  plot(predict(model), MASS::studres(model), main="Residuals by Predicted value", xlab="Unstandardized Predicted Values", ylab="Studentized Residuals")
  cat("Correlation Matrix for model (correlation >.70 indicates severe multicollinearity)\n\n")
  print(checkList$Correlation.Matrix)
  cat("\nVariance inflation factor (<10 desired):\n\n")
  print(checkList$Var.Inf.Factor)
  cat("\nStandardized Residuals (observations > 3.00 problematic):\n\n")
  if(outliers == "significant") {
    res <- unlist(checkList$Stand.Residuals)
    res1 <- res[res >= 3.0]
    res2 <- res[res <= -3.0]
    res <- c(res1, res2)
  } else if (outliers == "all") {
    res <- sort(abs(checkList$Stand.Residuals), decreasing = T)
  }
  if(length(res) > 0) {
    print(res)
  } else {
    cat("No significant outliers\n")
  }
  cat("\nCook's distance (values >.2 problematic):\n\n")
  if(outliers == "significant") {
    cookd <- unlist(checkList$Cooks.Dist)
    cookd <- cookd[cookd >= 0.2]
  } else if (outliers == "all") {
    cookd <- sort(checkList$Cooks.Dist, decreasing = T)
  }
  if(length(cookd) > 0){
    print(cookd)
  } else {
    cat("No significant Cook's D values\n")
  }
  hist(checkList$Stand.Residuals, prob=T, breaks = 30, main = "Plot of Std Residuals", xlab="Std Residuals")
  distmean <- mean(MASS::stdres(model))
  distsd <- sd(MASS::stdres(model))
  #curve(dnorm(x, distmean, distsd, col="darkblue", lwd=2, add=TRUE, yaxt="n"))
  cat("\nNormality of standardized model residuals:", " Shapiro-Wilk (p-value): ", checkList$Normality.Resids, "\n\n")
  plot(ppoints(length(MASS::stdres(model))), sort(checkList$Prob.Dist), main = "PP Plot", xlab = "Observed Probability", ylab = "Expected Probability")
  abline(0,1)
  }
  # Output tables
  cat("Model change statistics\n\n")
  summary <- model_summary_table(models, formulas)
  cat("\nModel Coefficients\n\n")
  coefficients <- model_coefficient_table(models)
  results <- list(Summary=summary, Coefficients=coefficients, Checks=checkList, SummaryDF=as.data.frame(summary), CoefficientsDF=as.data.frame(coefficients))
  invisible(results)
}

#' Hierarchical regression: model summary output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @param formulas Formula list produced by \code{create_formula_objects}.
#' @details Creates table output to summarize model statistics for all models in
#'   a hierarchical regression analysis.
#' @examples
#' freeny_model_formulas <- create_formula_objects("y",
#' c("lag.quarterly.revenue"), c("price.index"))
#' freeny_models <- create_model_objects(freeny_model_formulas,
#' dataset = freeny)
#' model_summary_table(freeny_models, freeny_model_formulas)
#' @export
model_summary_table <- function(models, formulas) {
  model_para_list <- list()
  model_para_list[[1]] <- list(
  round(sqrt(summary(models[[1]])$r.squared), 4),
  round(summary(models[[1]])$r.squared, 4),
  round(summary(models[[1]])$adj.r.squared, 4),
  round(summary(models[[1]])$sigma, 4),
  round(summary(models[[1]])$r.squared, 4),
  round(as.numeric(summary(models[[1]])$fstatistic[1]), 4),
  as.numeric(summary(models[[1]])$fstatistic[2]),
  as.numeric(summary(models[[1]])$fstatistic[3]),
  round(1 - pf(
  as.numeric(summary(models[[1]])$fstatistic[1]),
  as.numeric(summary(models[[1]])$fstatistic[2]),
  as.numeric(summary(models[[1]])$fstatistic[3])
  ), 4),
  gtools::stars.pval(1 - pf(
    as.numeric(summary(models[[1]])$fstatistic[1]),
    as.numeric(summary(models[[1]])$fstatistic[2]),
    as.numeric(summary(models[[1]])$fstatistic[3])
  ))
  )
  for(i in 2:length(models)) {
    comparison <- modelCompareMod(models[[i-1]], models[[i]], printOutput = F)
    model_para_list[[(length(model_para_list) + 1)]] <- list(
    R_1 = round(sqrt(summary(models[[i]])$r.squared), 4),
    R2_1 = round(summary(models[[i]])$r.squared, 4),
    R2_Adj_1 = round(summary(models[[i]])$adj.r.squared, 4),
    SE_1 = round(summary(models[[i]])$sigma, 4),
    delta12 = round(comparison$DeltaR2, 4),
    F12 = round(comparison$Fstat, 4),
    DF1_2 = as.numeric(comparison$nDF),
    DF2_2 = as.numeric(comparison$dDF),
    F12p = round(comparison$p, 4),
    sig = gtools::stars.pval(comparison$p)
  )
  }
  output_matrix <- matrix(nrow=length(models), ncol=10)
  model_labels <- list()
  for(i in 1:length(models)){
    output_matrix[i,] <- unlist(model_para_list[[i]])
    model_labels[length(model_labels) + 1] <- as.character(paste("Model", i))
  }
  colnames(output_matrix) <- c("R", "R^2", "Adj R^2", "SE Est.", "Delta R^2", "F Change", "df1", "df2", "p Fch", "Sig")
  rownames(output_matrix) <- unlist(model_labels)
  output_matrix <- as.data.frame(output_matrix)
  #for(i in c(1:6,9)){
  #  round(output_matrix[,i])
  #}
  print(output_matrix, digits = 3)
  forms <- c()
  for(i in 1:length(formulas)){
    #forms <- c()
    cat("Model", i, ":", deparse(formulas[[i]], width.cutoff = 500), "\n")
    form <- paste(deparse(formulas[[i]], width.cutoff = 500))
    forms <- c(forms, form)
  }
  Results <- list(R=output_matrix$R, R2=output_matrix$`R^2`, AdjR2=output_matrix$`Adj R^2`, SE=output_matrix$`SE Est.`, DeltaR2=output_matrix$`Delta R^2`, Fch=output_matrix$`F Change`, DF1=output_matrix$df1, DF2=output_matrix$df2, SigFch=output_matrix$`p Fch`, pval=output_matrix$`Sig`, Formula=forms)
  attr(Results$R, "names") <- rownames(output_matrix)
  attr(Results$R2, "names") <- rownames(output_matrix)
  attr(Results$AdjR2, "names") <- rownames(output_matrix)
  attr(Results$SE, "names") <- rownames(output_matrix)
  attr(Results$DeltaR2, "names") <- rownames(output_matrix)
  attr(Results$Fch, "names") <- rownames(output_matrix)
  attr(Results$DF1, "names") <- rownames(output_matrix)
  attr(Results$DF2, "names") <- rownames(output_matrix)
  attr(Results$SigFch, "names") <- rownames(output_matrix)
  attr(Results$pval, "names") <- rownames(output_matrix)
  invisible(Results)
}

#' Hierarchical regression: Coefficient table output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @details Creates table output to summarize model coefficients for all models
#'   in a hierarchical regression analysis.
#' @examples
#' freeny_model_formulas <- create_formula_objects("y", c("lag.quarterly.revenue")
#' , c("price.index"))
#' freeny_models <- create_model_objects(freeny_model_formulas,
#' dataset = freeny)
#' model_coefficient_table(freeny_models)
#' @export
model_coefficient_table <- function(models){
  tables <- broom::tidy(models[[1]])
  tables <- cbind(Model = "Model 1", tables)
  for(i in 2:length(models)){
    tables <- rbind(tables, cbind(Model = paste("Model", i), broom::tidy(models[[i]])))
  }
  tables$p.value <- round(tables$p.value, digits = 4)
  sig.stars <- gtools::stars.pval(tables$p.value)
  tables$sig <- sig.stars
  print(tables, row.names = F, digits = 4)
  results <- as.list(tables)
  attr(results$estimate, "names") <- results$term
  attr(results$std.error, "names") <- results$term
  attr(results$statistic, "names") <- results$term
  attr(results$p.value, "names") <- results$term
  invisible(results)
}

#' Binary Logistic Regression: Model Output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @param formulas A list of model formulas, generated by \code{create_formula_objects}.
#' @details Creates output for results of hierarchical binary logisitic
#'   regression models.
#' @examples
#' formulas <- create_formula_objects("am", c("hp", "mpg"), c("disp"),
#' c("drat"))
#' mtcars_models <- create_model_objects(formulas, data=mtcars,
#' type="binomial")
#' model_output_binomial(mtcars_models, formulas)
#' @export
model_output_binomial <- function(models, formulas) {
  cat("\nModel Summary Table: Pseudo R^2\n\n")
  summary <- model_summary_table_binomial(models)
  cat("\nModel Coefficient Table\n\n")
  coefficients <- model_coefficient_table_binomial(models)
  full_model <- models[[length(models)]]
  response <- full_model$model[,1]
  cat("\nClassification Table\n\n")
  class_table <- classification_table(full_model, response)
  forms <- c()
  for(i in 1:length(formulas)){
    #forms <- c()
    cat("Model", i, ":", deparse(formulas[[i]], width.cutoff = 500), "\n")
    form <- paste(deparse(formulas[[i]], width.cutoff = 500))
    forms <- c(forms, form)
  }
  results <- list(Summary=summary, Coefficients=coefficients, Class_Table=class_table, Formulas=forms)
  invisible(results)
}

#' Binary Logistic Regression: Classification Table
#'
#' @param model A binary logistic regression model object.
#' @param response The dependent variable in \code{model}.
#' @details Creates classification table for binary logistic regresison model
#'   using optimal cut point for accuracy.
#' @examples
#' formulas <- create_formula_objects("am", c("hp", "mpg"), c("disp"),
#' c("drat"))
#' mtcars_models <- create_model_objects(formulas, data=mtcars,
#' type="binomial")
#' last_model <- mtcars_models[[length(mtcars_models)]]
#' classification_table(last_model, last_model$model[,1])
#' @export
classification_table <- function(model, response){
  frame <- dplyr::arrange(data.frame(predict.glm(model), response), dplyr::desc(predict.glm(model)))
  predictor <- ROCR::prediction(predict.glm(model), response)
  accuracy <- ROCR::performance(predictor, measure = "acc")
  cutpoint <- which.max(accuracy@y.values[[1]])
  cutpoint
  frame$cut <- NA
  for(i in 1:cutpoint){
    frame$cut[i] <- 1
  }
  for(i in (cutpoint+1):nrow(frame)){
    frame$cut[i] <- 0
  }
  predone <- frame[0:cutpoint,]
  predzero <- frame[(cutpoint + 1):nrow(frame),]
  totp <- nrow(predone)
  totn <- nrow(predzero)
  tp <- nrow(dplyr::filter(predone, response == 1))
  fp <- totp - tp
  tn <- nrow(dplyr::filter(predzero, response == 0))
  fn <- totn - tn
  specificity <- tp / totp
  sensitivity <- tn / totn
  tot_accuracy <- (tp + tn) / (totp + totn)
  tabs <- table(frame$cut, frame$response, dnn = c("Predict", "Actual"))
  print(tabs)
  cat("Specificity: ", specificity, "\n")
  cat("Sensitivity: ", sensitivity, "\n")
  cat("Total Accuracy: ", tot_accuracy, "\n")
  table <- as.data.frame(tabs)
  }

#' Binary Logistic Regression: Model Summary Output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @details Creates summary table with pseudo R^2 values for the binary logistic
#'   regression model objects
#' @examples
#' formulas <- create_formula_objects("am", c("hp", "mpg"), c("disp"),
#' c("drat"))
#' mtcars_models <- create_model_objects(formulas, data=mtcars,
#' type="binomial")
#' last_model <- mtcars_models[[length(mtcars_models)]]
#' @export
model_summary_table_binomial <- function(models){
  output_matrix <- matrix(nrow = length(models), ncol = 4)
  model_labels <- list()
  for(i in 1:length(models)){
  output_matrix[i,] <- BaylorEdPsych::PseudoR2(models[[i]])[1:4]
  model_labels[[i]] <- as.character(paste("Model", i))
  }
  colnames(output_matrix) <- c("McFadden's", "Adj McFadden's", "Cox-Snell", "Nagelkerke")
  rownames(output_matrix) <- unlist(model_labels)
  print(output_matrix)
}

#' Binary Logistic Regression: Coefficient Table Output
#'
#' @param models A list of \code{lm} model objects.  A set of model objects
#'   created by \code{create_model_object}.
#' @details Creates table output to summarize model coefficients for all models
#'   in a hierarchical binary logistic regression analysis.
#' @examples
#' formulas <- create_formula_objects("am", c("hp", "mpg"), c("disp"),
#' c("drat"))
#' mtcars_models <- create_model_objects(formulas, data=mtcars,
#' type="binomial")
#' \dontrun{model_summary_table_binomial(mtcars_models)}
#' @export
model_coefficient_table_binomial <- function(models){
  full_model <- models[[length(models)]]
  full_model_summary <- summary(full_model)
  full_model_table <- broom::tidy(full_model)
  #model_terms <- attr(full_model_summary$terms, "term.labels")
  #model_terms <- c("Constant", model_terms)
  model_terms <- attr(full_model$coefficients, "names")
  output_table <- data.frame(model_terms)
  for(i in 1:length(models)){
   model_summary <- summary(models[[i]])
   model_table <- broom::tidy(models[[i]])
   coefs <- round(model_table$estimate, digits = 4)
   SE <- round(model_table$std.error, digits = 4)
   wald <- c()
   p.value <- c()
   for(j in 1:length(SE)){
   wald.value <- aod::wald.test(Sigma = vcov(models[[i]]), b = coef(models[[i]]), Terms = j)
   wald <- c(wald, wald.value$result$chi2[1])
   p.value <- c(p.value, paste(round(wald.value$result$chi2[3], digits = 4), gtools::stars.pval(wald.value$result$chi2[3]), sep = " "))
   }
   wald <- round(wald, digits = 4)
   #p.value <- round(p.value, digits = 4)
   df <- data.frame(coefs, SE, wald, p.value)
   if(nrow(df) < nrow(full_model_table)){
   for(i in seq((nrow(df) + 1), (nrow(full_model_table)))){
     df <- rowr::insertRows(df, data.frame(list('--', '--', '--', '--')), i)
   }
   }
   output_table <- data.frame(output_table, df)
  }
  print(output_table, digits = 4)
  #print(df, digits = 4)
}

#' Modified modelCompare function from lmSupport package.  Modified to suppress
#' print output.
#'
#' @param ModelC A model \code{lm} object.
#' @param ModelA A model \code{lm} object.
#' @param printOutput Boolean parameter, if \code{TRUE}, print output is
#'   supressed from modelCompare function.
#' @details This is a modification of the modelCompare function that allows
#'   print output to be suppressed.
#' @export
modelCompareMod <- function(ModelC, ModelA, printOutput = T) {
    sseC = sum(residuals(ModelC)^2)
    sseA = sum(residuals(ModelA)^2)

    pC = length(coef(ModelC))
    pA = length(coef(ModelA))
    if (!(pA > pC))  stop('Invalid model comparison:  modelA does not have more parameters than modelC')

    #Added the next three lines to check whether the terms in model C are a subset of the terms in model A
    termsC <- attr(terms(ModelC), "term.labels")
    termsA <- attr(terms(ModelA), "term.labels")

    if (!all(termsC %in% termsA))  stop('Invalid model comparison:  modelC is not a subset of modelA')

    nC = ModelC$df.residual + pC
    nA = ModelA$df.residual + pA
    if (!(nC == nA))  stop('Invalid model comparison:  ModelA and ModelC have different N')

    nDF = pA - pC
    dDF = nA - pA
    FStat=   ((sseC - sseA) / (pA-pC)) / (sseA / (nA-pA))

    p = pf(FStat,nDF, dDF, lower.tail = FALSE)

    PRE = (sseC - sseA) / sseC
    DeltaR2 = summary(ModelA)$r.squared - summary(ModelC)$r.squared

    #print output
    if(printOutput == TRUE){
    cat('SSE (Compact) = ', sseC, '\n', sep=' ')
    cat('SSE (Augmented) = ', sseA,  '\n', sep=' ')
    cat('Delta R-Squared = ', DeltaR2,  '\n', sep=' ')
    cat('Partial Eta-Squared (PRE) = ', PRE,  '\n', sep=' ')
    cat('F(', nDF, ',', dDF, ') = ', FStat, ', ', 'p = ', p, '\n', sep='')
    }
    Results = list(sseC=sseC, sseA=sseA, pC=pC, pA=pA, nDF=nDF, dDF=dDF, Fstat=FStat, p=p,  PRE=PRE, DeltaR2=DeltaR2)
    invisible(Results)  #return but dont print list
}

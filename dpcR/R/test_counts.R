#' Test counts
#' 
#' The test for comparing counts from two or more digital PCR experiments.
#' 
#' @aliases test_counts
#' @param input object of class \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}}
#' with "nm" type.
#' @param model may have one of following values: \code{binomial}, \code{poisson},
#' \code{prop}, \code{ratio}. See Details.
#' @param conf.level confidence level of the intervals and groups.
#' @details 
#' \code{test_counts} incorporates two different approaches to models: GLM (General Linear 
#' Model) and multiple pair-wise tests. The GLM fits counts data from different 
#' digital PCR experiments using quasibinomial or quasipoisson \code{\link[stats]{family}}. 
#' Comparisons between single experiments utilize Tukey's contrast and multiple t-tests 
#' (as provided by function \code{\link{glht}}).
#' 
#' In case of pair-wise tests, (\code{\link[rateratio.test]{rateratio.test}} or 
#' \code{\link[stats]{prop.test}}) are used to compare all pairs of experiments. The 
#' p-values are adjusted using the Benjamini & Hochberg method (\code{\link[stats]{p.adjust}}).
#' Furthermore, confidence intervals are simultaneous.
#' 
#' @note Mean number of template molecules per partition and its confidence intervals will 
#' vary depending on input.
#' @export
#' @seealso
#' Functions used by \code{test_counts}: 
#' \itemize{
#' \item \code{\link[stats]{glm}},
#' \item \code{\link[multcomp]{glht}},
#' \item \code{\link[multcomp]{cld}},
#' \item \code{\link[stats]{prop.test}},
#' \item \code{\link[rateratio.test]{rateratio.test}}
#' }
#' 
#' GUI presenting capabilities of the test: \code{\link{test_counts_gui}}.
#' 
#' @return an object of class \code{\linkS4class{count_test}}.
#' @author Michal Burdukiewicz, Stefan Roediger, Piotr Sobczyk.
#' @references Bretz F, Hothorn T, Westfall P, \emph{Multiple comparisons using R}. 
#' Boca Raton, Florida, USA: Chapman & Hall/CRC Press (2010).
#' @examples
#' #be warned, the examples of test_counts are time-consuming
#' \dontrun{
#' adpcr1 <- sim_adpcr(m = 10, n = 765, times = 1000, pos_sums = FALSE, n_panels = 3)
#' adpcr2 <- sim_adpcr(m = 60, n = 550, times = 1000, pos_sums = FALSE, n_panels = 3)
#' adpcr2 <- rename_dpcr(adpcr2, exper = "Experiment2")
#' adpcr3 <- sim_adpcr(m = 10, n = 600, times = 1000, pos_sums = FALSE, n_panels = 3)
#' adpcr3 <- rename_dpcr(adpcr3, exper = "Experiment3")
#' 
#' #compare experiments using binomial regression
#' two_groups_bin <- test_counts(bind_dpcr(adpcr1, adpcr2), model = "binomial")
#' summary(two_groups_bin)
#' plot(two_groups_bin)
#' #plot aggregated results
#' plot(two_groups_bin, aggregate = TRUE)
#' #get coefficients
#' coef(two_groups_bin)
#' 
#' #this time use Poisson regression
#' two_groups_pois <- test_counts(bind_dpcr(adpcr1, adpcr2), model = "poisson")
#' summary(two_groups_pois)
#' plot(two_groups_pois)
#' 
#' #see how test behaves when results aren't significantly different
#' one_group <- test_counts(bind_dpcr(adpcr1, adpcr3))
#' summary(one_group)
#' plot(one_group)
#' }

test_counts <- function(input, model = "ratio", conf.level = 0.95) { 
  if(!(model %in% c("binomial", "poisson", "prop", "ratio")))
    stop("Must must have one of following values: 'binomial', 'poisson', 'ratio' or 'prop'.")
  
  
  if(model %in% c("prop", "ratio")) {
    test_function <- if(model == "prop") prop.test else rateratio.test
    
    #extract number of positive partitions
    positives <- if(slot(input, "type") == "tnp") {
      slot(input, ".Data")[1, ]
    } else {
      colSums(input > 0, na.rm = TRUE)
    }
    total <- slot(input, "n")
    #change sequence of rows to create output similar to glm
    test_ids <- combn(1L:length(total), 2)[c(2, 1), ]
    
    #for only two experiments
    if(!is.matrix(test_ids))
      test_ids <- as.matrix(test_ids)
    #perform one-against-one multiple proportion tests
    all_combns <- apply(test_ids, 2, function(i)
      test_function(positives[i], total[i]))
    
    #adjust p-values
    p_vals <- p.adjust(sapply(all_combns, function(i)
      i[["p.value"]]), method = "BH")
    
    #     #for only two experiments
    #     if(!is.matrix(p_vals))
    #       p_vals <- as.matrix(p_vals)
    
    #values of chi square statistic
    statistics <- sapply(all_combns, function(i)
      i[["statistic"]])

    #what should be in res, depends on the model
    test_res <-  if(model == "prop") {
      matrix(c(statistics, p_vals), ncol = 2,
             dimnames = list(apply(matrix(names(positives)[test_ids], ncol = 2, 
                                          byrow= TRUE),1, function(i) 
                                            paste(i, collapse = " - ")),
                             c("X_squared", "p_value")))
    } else {
      matrix(p_vals, ncol = 1,
             dimnames = list(apply(matrix(names(positives)[test_ids], ncol = 2, 
                                          byrow= TRUE),1, function(i) 
                                            paste(i, collapse = " - ")),
                             c("p_value")))
    }
    
    
    #split data in groups
    #calculate confidence intervals using Sidak's unequality
    group_vals <- fl(binom.confint(positives, 
                                   total, 
                                   conf.level = (conf.level)^(1/ncol(input)),
                                   "wilson")[, 4L:6])
    
    #only unsignif in reality
    only_signif <- test_ids[, p_vals > 1 - conf.level]
    if(!is.matrix(only_signif))
      only_signif <- as.matrix(only_signif)
    #every experiment belongs to different group when dim(only_signif)[2] == 0
    groups <- if(dim(only_signif)[2] == 0) {
      as.list(1L:length(total))
    } else {
      unique(lapply(1L:length(total), function(i)
        sort(unique(as.vector(only_signif[, as.logical(colSums(only_signif == i))])))))
    }
    
    #detect single experiment groups - they are empty before this
    groups_length <- sapply(groups, length)
    if(0 %in% groups_length) {
      groups <- groups[groups_length != 0]
      groups <- c(groups, as.list((1L:length(total))[!(1L:length(total) %in% 
                                                         unlist(groups))]))
    }
    
    
    group_matrix <- sapply(1L:length(total), function(experiment) 
      sapply(groups, function(single_group) experiment %in% single_group))
    #all experiments in one group
    if(is.null(dim(group_matrix)))
      group_matrix <- matrix(group_matrix, nrow = 1)
    #name groups using the abc convention and at the same time reorder them along to value
    dimnames(group_matrix) <- list((1L:length(groups))[order(sapply(groups, function(single_group) 
                                     mean(group_vals[single_group, 1])))], 
                                   names(positives))
    
    group_coef <- data.frame(apply(group_matrix, 2, function(i) 
      paste(sort(as.numeric(names(i[which(i)]))), collapse = ".")), 
      group_vals)
    colnames(group_coef) <- c("group", "lambda", "lambda.low", "lambda.up")
    
  } else {
    
    if(slot(input, "type") == "tnp")
      stop("GLM does not work with 'tnp' type.")
    
    #choose proper family
    if (model == "binomial") {
      if(!(slot(input, "type") %in% c("tp")))
        #binarize input which isn't already binary
        input <- binarize(input)
      #family for model
      fam <- quasibinomial(link = "log")
      #function transforming coefficients of model to lambdas
      trans_fun <- function(x) fl(exp(x))
    } 
    
    if (model == "poisson") {
      
      #family for model
      fam <- quasipoisson(link = "log")
      #function transforming coefficients of model to lambdas
      trans_fun <- function(x) exp(x)
    }
    
    #dpcr version of melt
    n_vector <- slot(input, "n")
    m_dpcr <- do.call(rbind, lapply(1L:length(n_vector), function(i) {
      vals <- input[1L:n_vector[i], i]
      data.frame(experiment = rep(colnames(input)[i], length(vals)), values = vals)
    }))
    
    
    #fit model
    fit <- glm(values ~ experiment + 0, data = m_dpcr, family = fam)
    
    #do multiple comparision
    multi_comp <- glht(fit, linfct = mcp(experiment = "Tukey"))
    
    coefs <- summary(fit)[["coefficients"]][, 1:2]
    lambdas <- trans_fun(matrix(c(coefs[, 1], 
                                  coefs[, 1] - coefs[, 2], 
                                  coefs[, 1] + coefs[, 2]), ncol = 3))
    
    summ_mc <- suppressWarnings(summary(multi_comp))
    groups <- suppressWarnings(cld(multi_comp)[["mcletters"]][["LetterMatrix"]])
    if(is.matrix(groups)) {
      #more than one group
      groups_vector <- apply(groups, 1, function(i)
        paste0(colnames(groups)[i], collapse = ""))
    } else {
      #only one group
      groups_vector <- rep("a", ncol(input))
      names(groups_vector) <- colnames(input)
    }
    
    group_coef <- data.frame(groups_vector, lambdas)
    colnames(group_coef) <- c("group", "lambda", "lambda.low", "lambda.up")
    rownames(group_coef) <- colnames(input)
    test_res <- cbind(t = summ_mc[["test"]][["tstat"]], 
                      p_value = summ_mc[["test"]][["pvalues"]])
  }
  
  summ <- summary(input, print = FALSE)[["summary"]]
  
  group_coef <- cbind(group_coef, summ[summ[["method"]] == "bhat", c("experiment", "replicate", "k", "n")])
  
  new("count_test", group_coef = group_coef, test_res = test_res, 
      model = model)
}
#' Define an ANOVA Design
#'
#' Constructs an "design.anova" object required by the \code{\link{boot.power.anova}} function.
#'
#' Based on the supplied within-subjects factors and between-subjects factors, this function
#' constructs all conditions of the ANOVA design and opens two dialog windows querying for
#' the expected correlation matrix and cell means (+ standard deviations) for all conditions.
#'
#' The first dialog window queries for the correlation matrix of the conditions. If you have
#' a pure between-subjects design, you may instantly close this window. Otherwise, enter the
#' expected correlations between all conditions that include within-subjects manipulations.
#' Using the "default_within_correlation" parameter, a default value can be set. You should
#' fill in only the lower triangle of the correlation matrix and only the values not containing
#' NAs.
#'
#' The second dialog window queries for the means and standard deviations expected for each
#' condition.
#'
#' Use the "save_input_as" parameter in order to define a file name prefix of the files where
#' the function saves your input values. This will populate the dialog windows with the saved
#' values on the next execution of this function. If the parameter is NULL, input values will
#' not be saved.
#'
#' @param between list, between-subjects factors including its levels
#' @param within list, within-subjects factors including its levels
#' @param default_within_correlation numeric, default within-subjects correlation the correlation
#'                                            matrix is populated with (for designs including
#'                                            within-subjects factors)
#' @param save_input_as character, file name prefix of the files the input values entered by
#'                                 the user are save to. File names are constructed as
#'                                 paste0(save_input_as,"_cor_matrix.csv") and
#'                                 paste0(save_input_as,"_means_and_sds.csv")
#' @param silent_load boolean, FALSE (default): always show input dialogs (even if data was successfully
#'                             loaded from a file); TRUE: show input dialogs only if file did not yet exist
#'                             and break with error if data from file does not match the design
#' @return object of type design.anova
#' @seealso \code{\link{boot.power.anova}}
#' @examples
#' \dontrun{
#' design <- design.anova(
#'    between = list(age = c("young","old"),
#'                   sex = c("male","female")),
#'    within = list(condition = c("cond1","cond2","cond3")),
#'    default_within_correlation = 0.7,
#'    save_input_as = "myexp1",
#'    silent_load = T
#' )
#' }
#'
#' @export
design.anova <- function(
  between = list(),
  within = list(),
  default_within_correlation = 0,
  save_input_as = NULL,
  silent_load = FALSE
) {
  if (! is.list(between)) stop("between must be a list")
  if (! is.list(within)) stop("within must be a list")
  if (length(between)==0 && length(within)==0) stop("between and/or within structure must be specified")
  if (default_within_correlation < 0 || default_within_correlation > 1) stop("default_within_correlation must be in range 0-1")

  # Replace spaces in factor levels by underscore
  between <- lapply(between, stringr::str_replace_all, " ", "_")
  within <- lapply(within, stringr::str_replace_all, " ", "_")

  # Replace dots in factor levels by underscore
  between <- lapply(between, stringr::str_replace_all, "\\.", "_")
  within <- lapply(within, stringr::str_replace_all, "\\.", "_")

  expanded <- expand.grid(c(within,between))
  factor_names <- names(expanded)
  vars <- apply(expanded,
                1,
                paste,collapse=".")

  if (! is.null(save_input_as)) {
    file_name_cor_matrix <- paste0(save_input_as,"_cor_matrix.csv")
    file_name_means_and_sds <- paste0(save_input_as,"_means_and_sds.csv")
  } else {
    file_name_cor_matrix <- NULL
    file_name_means_and_sds <- NULL
  }

  cor_matrix <- matrix()
  if (! is.null(file_name_cor_matrix)) tryCatch({cor_matrix <- as.matrix(utils::read.csv(file_name_cor_matrix,row.names=1))}, warning=function(w){}, error=function(e){})
  new_cor_matrix_created <- FALSE

  # Create new correlation matrix because none saved or saved matrix does not fit to design
  if (! (identical(rownames(cor_matrix),vars) && identical(colnames(cor_matrix),vars))) {
    if (silent_load && !is.na(cor_matrix)) { # !is.na(cor_matrix): file exists and was not empty but content does not match design
      stop(paste0("Content of ",file_name_cor_matrix," does not match the specified design. Disable silent_load or delete file to enter new data."))
    }
    new_cor_matrix_created <- TRUE

    cor_matrix <- diag(1,length(vars),length(vars))
    rownames(cor_matrix) <- vars
    colnames(cor_matrix) <- vars
    for (row in vars) {
      for (col in vars) {
        if (row != col) {
          # Cells including between manipulation: NA (replaced by 0 later on)
          if (! identical(stringr::str_split_fixed(row,"\\.",length(names(within))+1)[length(names(within))+1],
                          stringr::str_split_fixed(col,"\\.",length(names(within))+1)[length(names(within))+1]) ) {
            cor_matrix[row,col] <- NA
          } else { # Cells with pure within manipulations
            cor_matrix[row,col] <- default_within_correlation
          }
        }
      }
    }
    # Correlation Matrix is symmetric: User fills in lower triangle only (copied to upper triangle later on)
    gdata::upperTriangle(cor_matrix) <- NA
  }

  # Show input dialog if silent_load is disabled. Also show with silent_load if new data was created (no old file exists)
  if (! silent_load || new_cor_matrix_created) {
    cor_matrix <- utils::edit(cor_matrix)

    if (! is.null(file_name_cor_matrix)) utils::write.csv(cor_matrix,file_name_cor_matrix)
  }

  gdata::upperTriangle(cor_matrix) <- gdata::upperTriangle(t(cor_matrix))
  cor_matrix[is.na(cor_matrix)] <- 0

  if (! all(diag(cor_matrix)==1)) stop("All values on diagonal must be 1")

  means_sds <- data.frame()
  if (! is.null(file_name_means_and_sds)) tryCatch({means_sds <- utils::read.csv(file_name_means_and_sds,row.names=1)}, warning=function(w){}, error=function(e){})
  new_means_sds_created <- FALSE

  if (! (identical(rownames(means_sds),vars))) {
    if (silent_load && length(means_sds) > 0) { # length(means_sds) > 0: file exists and was not empty but content does not match design
      stop(paste0("Content of ",file_name_means_and_sds," does not match the specified design. Disable silent_load or delete file to enter new data."))
    }
    new_means_sds_created <- TRUE

    means_sds <- data.frame(mean=numeric(length(vars)),
                            sd  =numeric(length(vars)))
    rownames(means_sds) <- vars
  }

  # Show input dialog if silent_load is disabled. Also show with silent_load if new data was created (no old file exists)
  if (! silent_load || new_means_sds_created) {
    means_sds <- utils::edit(means_sds)

    if (! is.null(file_name_means_and_sds)) utils::write.csv(means_sds,file_name_means_and_sds)
  }

  if (sum(means_sds$sd <= 0) > 0) {
    stop("All SD must be >= 0")
  }

  diag(cor_matrix) <- means_sds$sd
  cov_matrix <- lme4::sdcor2cov(cor_matrix)

  num_between_conds <- dim(expand.grid(c(between)))[1]
  if (num_between_conds == 0) { # pure within design
    num_between_conds <- 1
  }

  out <- list(
    between            =  between,
    within             =  within,
    num_between_conds  =  num_between_conds,
    factor_names          =  factor_names,
    cov_matrix         =  cov_matrix,
    means              =  means_sds$mean
  )

  class(out) <- "design.anova"

  return( out )
}

#' Bootstrap the Power of an ANOVA Design
#'
#' This function bootstraps the power of each effect in an ANOVA design for a given range
#' of sample sizes. Power is computed by randomly drawing samples from a multivariate normal
#' distribution specified according to the values supplied by the \code{\link{design.anova}}
#' object. Power is defined as the proportion of bootstrap iterations the p-values of each
#' effect lie below the supplied alpha level. Note that this function runs many ANOVAs which
#' might be slow for large sample size ranges or bootstrap iterations (see Details below).
#' Further note that this function does not check for assumptions such as sphericity.
#'
#' Note that this function requires the computation of many ANOVAs and therefore becomes
#' slow with increasing sample size ranges and bootstrap iterations. It is therefore suggested
#' to first use a very low number of bootstrap iterations, such as 10, in order to determine
#' a sensible sample size range for the power of interest. Once done, use this small sample
#' size range and dramatically increase the bootstrap iterations, such as 3000, in order to
#' determine more reliable power estimates. Because the power-by-samplesize function is
#' monotonically increasing, a zigzag of power values with increasing sample sizes indicates
#' that the selected bootstrap iterations are too low.
#'
#' @param design object of type \code{\link{design.anova}}
#' @param n_from numeric, lower boundary of sample size range (inclusive)
#'                        ; Refers to N per between condition
#' @param n_to numeric, upper boundary of sample size range (inclusive)
#'                      ; Refers to N per between condition
#' @param num_iterations_bootstrap numeric, number of bootstrap iterations for each sample size
#' @param alpha numeric, alpha level
#' @return list containing power-by-samplesize data.frames for each effect
#' @seealso \code{\link{design.anova}}, \code{\link{plot.power_by_samplesize}}
#' @examples
#' \dontrun{
#' design <- design.anova(
#'    between = list(age = c("young","old"),
#'                   sex = c("male","female")),
#'    within = list(condition = c("cond1","cond2","cond3")),
#'    default_within_correlation = 0.7
#' )
#'
#' power_by_samplesize <- boot.power.anova(
#'    design,
#'    n_from = 40,
#'    n_to = 60,
#'    num_iterations_bootstrap = 1000
#' )
#'
#' plot(power_by_samplesize,
#'      crit_power = 0.9,
#'      plot_dir = "power_plots")
#' }
#'
#' @export
boot.power.anova <- function(
  design,
  n_from, # N per between condition
  n_to, # N per between condition
  num_iterations_bootstrap,
  alpha = 0.05
) {

  if (! class(design) == "design.anova") stop("design must be an object of class design.anova")
  if (n_from < 2) stop("n_from must be >= 2")
  if (n_to < n_from) stop("n_to must be >= n_from")
  if (num_iterations_bootstrap < 1) stop("num_iterations_bootstrap must be >= 1")
  if (alpha < 0 || alpha > 1) stop("alpha must be in range 0-1")

  if (num_iterations_bootstrap < 1000) warning("Please consider increasing the number of bootstrap iterations in order to get more reliable power estimates.")

  boot_length <- (n_to-n_from+1) * num_iterations_bootstrap
  boot_p <- data.frame(n_total = rep(n_from:n_to,each=num_iterations_bootstrap)*design$num_between_conds, iteration = rep(1:num_iterations_bootstrap,boot_length/num_iterations_bootstrap))

  progress_bar <- utils::txtProgressBar(min = 1, max = boot_length, initial = 1, style = 3)
  progress_i <- 1

  for (cur_n_per_between_cond in n_from:n_to) {
    for (cur_iteration in 1:num_iterations_bootstrap) {
      dat <- MASS::mvrnorm(cur_n_per_between_cond, mu = design$means, Sigma = design$cov_matrix)
      dat <- data.frame(dat)
      dat$participant <- 1:dim(dat)[1]

      dat <- reshape2::melt(dat,id.vars="participant")

      variable_split <- stringr::str_split_fixed(dat$variable,"\\.",length(design$factor_names))
      for (i in 1:length(design$factor_names)) {
        dat[,design$factor_names[i]] <- factor(variable_split[,i])
      }

      for (between_var in names(design$between)) {
        dat$participant <- paste0(dat$participant,"_",dat[,between_var])
      }

      dat$participant <- factor(dat$participant)

      aovCall <- paste0("summary(aov(value ~ ",paste(c(names(design$within),names(design$between)),collapse="*"))
      if (length(names(design$within)) > 0) {
        aovCall <- paste0(aovCall," + Error(participant / (",paste(names(design$within),collapse="*"),"))")
      }
      aovCall <- paste0(aovCall,", dat))")

      aov_summary <- eval(parse(text=aovCall))

      if (dim(boot_p)[2] == 2) { # First run -> get effect names and populate with NA
        effect_names <- c()
        for (res in aov_summary) {
          if (length(res) == 1) {
            res <- res[[1]]
          }
          for (effect_name in stringr::str_trim(rownames(res))) {
            if (effect_name != "Residuals") {
              effect_names <- c(effect_names,effect_name)
            }
          }
        }
        effect_names <- effect_names[order(stringr::str_count(effect_names,':'))]

        for (effect_name in effect_names) {
          boot_p[,effect_name] <- NA
        }
      }

      cur_n_total <- length(unique(dat$participant))

      for (res in aov_summary) {
        if (length(res) == 1) {
          res <- res[[1]]
        }
        for (effect_name in rownames(res)) { # str_trim below instead of on rownames because value extraction from res works only with trailing whitespaces
          if (stringr::str_trim(effect_name) != "Residuals") {
            boot_p[boot_p$n_total==cur_n_total & boot_p$iteration==cur_iteration,stringr::str_trim(effect_name)] <- res[effect_name,"Pr(>F)"]
          }
        }
      }

      utils::setTxtProgressBar(progress_bar, progress_i)
      progress_i <- progress_i + 1
    }
  }

  # Power for each effect
  out <- list()

  n_total <- cur_effect <- NULL # for R CMD check

  for (effect_name in names(boot_p)[3:length(boot_p)]) {
    boot_p$cur_effect <- boot_p[,effect_name]
    eff <- plyr::ddply(boot_p,plyr::.(n_total),plyr::here(plyr::summarize), power = sum(cur_effect < alpha)/length(cur_effect))
    boot_p$cur_effect <- NULL

    out[[effect_name]] <- eff
  }

  close(progress_bar)

  class(out) <- "power_by_samplesize"

  return( out )
}

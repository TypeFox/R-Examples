#' @name aggSEM
#' @aliases aggSEM
#' @title Group-level structural equation model search.
#' @description Concatenates all individual-level data files and fits a group model to the data.
#' @usage
#' aggSEM(data   = "",
#'        out    = "",
#'        sep    = "",
#'        header = ,
#'        ar     = FALSE,
#'        plot   = TRUE,
#'        paths  = NULL)
#' @param data The path to the directory where the data files are located. Each file must contain one matrix for
#' each individual containing a T (time) by p (number of variables) matrix where the columns represent variables
#' and the rows represent time.
#' @param out The path to the directory where the results will be stored. This directory must be generated
#' by the user prior to running the function.
#' @param sep The spacing of the data files. "" indicates space-delimited, "/t" indicates tab-delimited, ","
#' indicates comma delimited.
#' @param header Logical. Indicate TRUE for data files with a header.
#' @param ar Logical. If TRUE, begins search for group model with autoregressive (AR) paths open. Defaults
#' to FALSE.
#' @param plot Logical. If TRUE, graphs depicting relations among variables of interest will automatically
#' be created. For aggregate-level plot, red paths represent positive weights and blue paths represent negative weights. Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which to begin model estimation. That is, Y~X indicates that Y is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation) by the column number. If a
#' header is used, variables should be referred to using variable names. To reference lag variables, "lag"
#' should be added to the end of the variable name with no separation. Defaults to NULL.
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{allBetas}} Matrix. Contains estimates for each path in the aggregate-level model. The row variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{allStdErrors}} Matrix. Contains standard errors for each path in the aggregate-level model. The row variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{allPathEstimates}} {Contains estimate, standard error, p-value, and z-value for each path for the concatenated data.}
#'  \item{\strong{summaryFit}} {Contains model fit information for the aggregate-level model.}
#'  \item{\strong{summaryPathsPlot}} Contains aggregate-level plot. Red paths represent positive weights and blue paths represent negative weights.
#' }
#' @author Stephanie Lane
#' @examples
#' data(ts1,ts2,ts3,ts4,ts5)
#' input.path  <- file.path(tempdir(),"input")
#' output.path <- file.path(tempdir(),"output")
#' dir.create(input.path)
#' dir.create(output.path)
#' write.table(ts1,file.path(input.path,"ts1.txt"),col.names=FALSE,row.names=FALSE)
#' write.table(ts2,file.path(input.path,"ts2.txt"),col.names=FALSE,row.names=FALSE)
#' write.table(ts3,file.path(input.path,"ts3.txt"),col.names=FALSE,row.names=FALSE)
#' write.table(ts4,file.path(input.path,"ts4.txt"),col.names=FALSE,row.names=FALSE)
#' write.table(ts5,file.path(input.path,"ts5.txt"),col.names=FALSE,row.names=FALSE)
#' aggSEM(data   = input.path,
#'        out    = output.path,
#'        sep    = "",
#'        header = FALSE,
#'        ar     = TRUE,
#'        plot   = TRUE,
#'        paths  = NULL)
#' @export
aggSEM <- function(data,
                   out,
                   sep,
                   header,
                   ar    = FALSE,
                   plot  = TRUE,
                   paths = NULL){

  agg.internal <- function(setup.out){

  subjects = setup.out$subjects
  varnames = setup.out$varnames
  syntax   = setup.out$syntax
  plot     = setup.out$plot
  files    = list.files(setup.out$data,full.names=TRUE)
  header   = setup.out$header
  sep      = setup.out$sep
  subgroup = setup.out$subgroup
  agg      = setup.out$agg
  plot     = setup.out$plot

  data.all <- data.frame()
  for (k in 1:subjects){
    data.file <- read.data(files[k],
                           header = header,
                           sep    = sep)
    data.all  <- rbind(data.all,data.file)
  }
  colnames(data.all) <- c(varnames)

  addind.out <- addind(done            = 0,
                       evaluate        = 1,
                       syntax          = syntax,
                       data.file       = data.all,
                       setup.out       = setup.out)

  evalind.out <- evalind(addind.out = addind.out,
                         setup.out  = setup.out,
                         data.file  = data.all)

  fixfitind.out <- fixfitind(setup.out   = setup.out,
                             evalind.out = evalind.out,
                             data.file   = data.all)

  final.fit.out <- final.fit(setup.out     = setup.out,
                             fixfitind.out = fixfitind.out,
                             data.file     = data.all,
                             k             = 1)

  all.elements <- final.fit.out$ind.elements
  all.fit      <- as.matrix(final.fit.out$ind.fit)
  all.fit[,1]  <- "all"
  colnames(all.fit) <- c("subject", "chisq", "df", "pval",
                         "rmsea", "srmr", "nnfi", "cfi", "status")
  row.names(all.elements) <- NULL
  list <- list("all.elements" = all.elements,
               "all.fit"      = all.fit,
               "all.syntax"   = NULL,
               "all.diff.sub" = FALSE)
  return(list)
}

  setup.out <- setup(data     = data,
                     sep      = sep,
                     header   = header,
                     out      = out,
                     plot     = plot,
                     ar       = ar,
                     paths    = NULL,
                     subgroup = FALSE,
                     agg      = TRUE,
                     deconvolve_hrf = FALSE,
                     control=list(deconvolve_method="bush"))

  agg.internal.out <- agg.internal(setup.out = setup.out)

  wrapup.out <- wrapup(indsem.internal.out = agg.internal.out,
                       setup.out           = setup.out)

  print.gimme.aggSEM(z=setup.out)
}

print.gimme.aggSEM <- function(z){
  writeLines("aggSEM finished running normally")
  writeLines(paste("output is stored in", z$out))
}
################################################################################
################################################################################

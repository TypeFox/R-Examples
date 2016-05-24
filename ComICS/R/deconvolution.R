#'@title DCQ - Digital Cell Quantifier
#'@name dcq
#'@aliases dcq
#'@description DCQ combines genome-wide gene expression data with an immune cell-type
#'  reference data to infer changes in the quantities immune cell subpopulations.
#'
#'@param reference_data a data frame representing immune cell expression profiles.
#'  Each row represents an expression of a gene, and each column represents a
#'  different immune cell type. \code{colnames} contains the name of each immune cell
#'  type and the \code{rownames} includes the genes' symbol. The names of each immune
#'  cell type and the symbol of each gene should be unique. Any gene with
#'  missing expression values must be excluded.
#'
#'@param mix_data a data frame representing RNA-seq or microarray
#'  gene-expression profiles of a given complex tissue. Each row represents an
#'  expression of a gene, and each column represents a different experimental sample.
#'  \code{colnames} contain the name of each sample and \code{rownames} includes the genes' symbol.
#'  The name of each individual sample and the symbol of each gene should be unique.
#'  Any gene with missing expression values should be excluded.
#'
#'@param marker_set data frames of one column, that includes a
#'  preselected list of genes that likely discriminate well
#'  between the immune-cell types given in the reference data.
#'
#'@param alpha_used,lambda_min parameters of the L1 and L2 regularization. It is generally
#' recommended to leave the default value. For more information about this
#' parameter, see the glmnet package.
#'
#'@param number_of_repeats using one repeat will generate only one output model.
#' Using many repeats, DCQ calculates a collection of models, and outputs the
#' average and standard deviation for each predicted relative cell quantity.
#'
#' @param precent_of_data in order to run the analysis over all the cell types use 1.0.
#' For bootstrap purposes, you can use part of the data (e.g, 0.5).
#'
#'@usage dcq(reference_data, mix_data, marker_set, alpha_used=0.05,
#'  lambda_min=0.2, number_of_repeats=3, precent_of_data=1.0)
#'
#'@return a list that contains two matrices
#'  \item{average}{a matrix that contains the average relative quantities for each cell
#'    type in everytest sample.}
#'  \item{stdev}{a matrix that contains the standard deviations over all repeats for
#'    each cell types in each test sample.}
#'
#'@examples
#'data(commons)
#'data(dcqEx)
#'results <- dcq(reference_data=immgen_dat, mix_data=lung_time_series_dat, marker_set=DCQ_mar)
#'
#'@references Altboum Z, Steuerman Y, David E, Barnett-Itzhaki Z, Valadarsky L, Keren-Shaul H, et al.
#'  Digital cell quantification identifies global immune cell dynamics during influenza
#'  infection. Mol Syst Biol. 2014;10: 720. doi:10.1002/msb.134947
#'

#' @export
dcq <- function(reference_data,
               mix_data,
               marker_set,
               alpha_used=0.05,
               lambda_min=0.2,
               number_of_repeats=3,
               precent_of_data=1.0) {

  set.seed(9999)

  #### selecting only specific markers that we want to use for our analysis ####
  print("selecting only specific markers that we want to use for our analysis")
  mix_data_markers_only <- as.matrix(mix_data[row.names(mix_data) %in% marker_set[,1], ])
  reference_data_markers_only <- as.matrix(reference_data[row.names(reference_data) %in% rownames(mix_data_markers_only), ])

  reference_data_markers_only <- reference_data_markers_only[order(match(rownames(reference_data_markers_only),
                                                                      rownames(mix_data_markers_only))), ]

  alphaParam <- alpha_used
  lambda_min <- lambda_min

  holder = c()

  list_of_results_mat = list()
  print(paste("running deconvolution on ", colnames(marker_set), sep=""))
  for (n in c(1:number_of_repeats)) {

    #### shuffling the matrices in the same manner each run ####

    order_in_run <- sample(1:ncol(reference_data_markers_only),
                           precent_of_data*ncol(reference_data_markers_only), replace=FALSE)

    holder <- rbind(holder, order_in_run)

    reference_data_markers_only_sampled <- reference_data_markers_only[, holder[n, ]]
    res_mat <- matrix(NA, nrow=ncol(mix_data_markers_only), ncol= ncol(reference_data_markers_only),
                      dimnames=list(colnames(mix_data_markers_only), colnames(reference_data_markers_only)))


    ##### perform elastic-net
    for (i in c(1:ncol(mix_data_markers_only))) {

      fit2 <- glmnet::glmnet(reference_data_markers_only_sampled,
                            mix_data_markers_only[, i], family = c('gaussian'),
                            alpha = alphaParam, nlambda=100, lambda.min.ratio=lambda_min);
      y <- fit2$beta[,100];
      res_mat[i, match(colnames(reference_data_markers_only_sampled),
                             colnames(res_mat))] <- y
    }

    ### transforming mat to the same original order
    list_of_results_mat <- append(list_of_results_mat, list(res_mat))

  }

  average_mat <- apply(simplify2array(list_of_results_mat), 1:2, mean, na.rm=T)
  stdev_mat <- apply(simplify2array(list_of_results_mat), 1:2, stats::sd, na.rm=T)

  return(list("average"=average_mat, "stdev"=stdev_mat))
}

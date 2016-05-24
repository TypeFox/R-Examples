#' Time-course Gene Set Analyis
#'
#' Wrapper function for performing gene set analysis of longitudinal RNA-seq data
#'
#'@param y a numeric matrix of size \code{n x G} containing the raw RNA-seq counts or
#'preprocessed expressions from \code{n} samples for \code{G} genes.
#'
#'@param x a numeric matrix of size \code{n x p} containing the model covariates from
#'\code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the \code{K} variables
#'to be tested
#'
#'@param genesets either a vector or a list of index or subscript that define
#'
#'@param indiv a vector of length \code{n} containing the information for
#'attributing each sample to one of the studied individuals. Coerced
#'to be a \code{factor}.
#'
#'@param Sigma_xi a matrix of size \code{K x K} containing the covariance matrix
#'of the \code{K} random effects.
#'
#'@param which_weights a character string indicating which method to use to estimate
#'the mean-variance relationship wheights. Possibilities are \code{"loclin"},
#'\code{"voom"} or \code{NULL} (in which case no weighting is performed).
#'Default is \code{"loclin"}.
#'See \code{\link{sp_weights}} and \code{\link{voom_weights}} for details.
#'
#'@param which_test a character string indicating which method to use to approximate
#'the variance component score test, either \code{"permutation"} or \code{"asymptotic"}.
#'Default is \code{"permutation"}.
#'
#'@param n_perm the number of perturbations
#'
#'@param preprocessed a logical flag indicating wether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param doPlot a logical flag indicating wether the mean-variance plot should be drawn.
#' Default is \code{FALSE}.
#'
#'@param bw a character string indicating the smoothing bandwidth selection method to use. See
#'\code{\link[stats]{bandwidth}} for details. Possible values are \code{"ucv"}, \code{"SJ"},
#'\code{"bcv"}, \code{"nrd"} or \code{"nrd0"}
#'
#'@param kernel a character string indicating which kernel should be used.
#'Possibilities are \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"},
#'\code{"triangular"}, \code{"biweight"}, \code{"tricube"}, \code{"cosine"},
#'\code{"optcosine"}. Default is \code{"gaussian"} (NB: \code{"tricube"} kernel
#'corresponds to the loess method).
#'
#'@param exact a logical flag indicating wether the non-parametric weights accounting
#'for the mean-variance relationship should be computed exactly or extrapolated
#'from the interpolation of local regression of the mean against the
#'variance. Default is \code{FALSE}, which uses interporlation (faster).
#'
#'@param padjust_methods multiple testing correction method used if \code{genesets}
#'is a list. Default is "BH", i.e. Benjamini-Hochberg procedure for contolling the FDR.
#'Other possibilities are: \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'\code{"bonferroni"} or \code{"BY"} (for Benjamini-Yekutieli procedure).
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}:
#'   \item \code{preprocessed}:
#'   \item \code{n_perm}:
#'   \item \code{pval}: associated p-value
#' }
#'
#'@seealso \code{\link{sp_weights}} \code{\link{vc_test_perm}} \code{\link{vc_test_asym}} \code{\link{p.adjust}}
#'
#'@references Agniel D, Hejblum BP, Variance component score test for
#'time-course gene set analysis of longitudinal RNA-seq data, \emph{submitted}, 2016.
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts. \emph{Genome
#'Biology}, 15(2), R29.
#'
#'@importFrom stats p.adjust
#'@export
tcgsa_seq <- function(y, x, phi, genesets,
                      indiv = rep(1, nrow(x)), Sigma_xi = diag(ncol(phi)),
                      which_test = c("permutation", "asymptotic"),
                      which_weights = c("loclin", "voom"),
                      n_perm = 1000,
                      preprocessed = FALSE, doPlot = TRUE,
                      bw = "nrd",
                      kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "tricube", "cosine", "optcosine"),
                      exact = FALSE,
                      padjust_methods = c("BH", "BY", "holm", "hochberg", "hommel", "bonferroni")){




  stopifnot(is.matrix(y))

  if(!preprocessed){
    y_lcpm <- apply(y, MARGIN=2,function(v){log2((v+0.5)/(sum(v)+1)*10^6)})
    preprocessed <- TRUE
  }else{
    y_lcpm <- y
  }
  rm(y)

  if(is.data.frame(x)){
    phi <- as.matrix(as.data.frame(lapply(phi, as.numeric)))
  }

  if(is.data.frame(phi)){
    phi <- as.matrix(as.data.frame(lapply(phi, as.numeric)))
  }

  if(det(crossprod(cbind(x, phi)))==0){
    stop("crossprod(x, phi) cannot be inversed. x and phi are likely colinear...")
  }

  if(length(padjust_methods)>1){
    padjust_methods <- padjust_methods[1]
  }
  stopifnot(padjust_methods %in% c("BH", "BY", "holm", "hochberg", "hommel", "bonferroni"))

  if(length(which_weights)>1){
    which_weights <- which_weights[1]
  }
  stopifnot(is.null(which_weights) || which_weights %in% c("loclin", "voom"))

  if(length(which_test)>1){
    which_test <- which_test[1]
  }
  stopifnot(which_test %in% c("asymptotic", "permutation"))

  w <-  switch(which_weights,
               loclin = sp_weights(x = x, y = t(y_lcpm), phi=phi,
                                   preprocessed = preprocessed, doPlot=doPlot,
                                   bw = bw, kernel = kernel,
                                   exact = exact),
               voom = voom_weights(x, y_lcpm, preprocessed = preprocessed, doPlot = doPlot),
               NULL = matrix(1, ncol=ncol(y_lcpm), nrow=nrow(y_lcpm))
  )


  if(which_test == "asymptotic"){
    n_perm <- NA
  }

  if(is.list(genesets)){

    if(class(genesets[[1]])=="character"){
      gene_names_measured <- rownames(y_lcpm)
      prop_meas <- sapply(genesets, function(x){length(intersect(x, gene_names_measured))/length(x)})
      if(sum(prop_meas)!=length(prop_meas)){
        warning("Some genes in the investigated genesets were not measured:\nremoving those genes form the geneset definition...")
        genesets <- lapply(genesets, function(x){x[which(x %in% gene_names_measured)]})
      }
    }

    if(which_test == "asymptotic"){
      rawPvals <- sapply(genesets, FUN = function(gs){
        vc_test_asym(y = y_lcpm[gs, ], x = x, indiv = indiv, phi = phi,
                     w = w[gs, ], Sigma_xi = Sigma_xi)$pval}
      )
    }
    else if(which_test == "permutation"){
      rawPvals <- sapply(genesets, FUN = function(gs){
        vc_test_perm(y = y_lcpm[gs, ], x = x, indiv = indiv, phi = phi,
                     w = w[gs, ], Sigma_xi = Sigma_xi, n_perm=n_perm)$pval
        }
      )
    }

    pvals <- data.frame("rawPval" = rawPvals, "adjPval" = stats::p.adjust(rawPvals, padjust_methods))
    if(!is.null(names(genesets))){
      rownames(pvals) <- names(genesets)
    }

  }else{

    if(class(genesets)=="character"){
      gene_names_measured <- rownames(y_lcpm)
      if((length(intersect(genesets, gene_names_measured))/length(x)) != 1){
        warning("Some genes in the investigated genesets were not measured:\n removing those genes form the geneset definition...")
        genesets <- genesets[which(genesets %in% gene_names_measured)]
      }
    }

    res_test <- switch(which_test,
                       asymptotic = vc_test_asym(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                 w = w[genesets, ], Sigma_xi = Sigma_xi),
                       permutation = vc_test_perm(y = y_lcpm[genesets, ], x = x, indiv = indiv, phi = phi,
                                                  w = w[genesets, ], Sigma_xi = Sigma_xi, n_perm = n_perm)
    )
    pvals <- data.frame("rawPval" = res_test$pval, "adjPval" = NA)
    padjust_methods <- NA

  }

  #print(pvals)

  return(list("which_test" = which_test, "preprocessed" = preprocessed, "n_perm" = n_perm,
              "genesets" = genesets, "pvals" = pvals))

}

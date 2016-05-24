#' Estimate genome-wide prediction accuacy using cross-validation
#'
#' @description \code{x.val} performs cross-validation (CV) to estimate the accuracy of genome-wide prediction (otherwise known as genomic selection) for a specific training population (TP), i.e. a set of individuals for which phenotypic and genotypic data is available. Cross-validation can be conducted via one of two methods within \code{x.val}, see \code{Details} for more information.
#' 
#'              NOTE - \code{x.val}, specifically \code{\link[BGLR]{BGLR}} writes and reads files to disk so it is highly recommended to set your working directory
#' @param G.in \code{Matrix} of genotypic data. First row contains marker names and the first column contains entry (taxa) names. Genotypes should be coded as follows: \itemize{
#'          \item \code{1}: homozygous for minor allele
#'          \item \code{0}: heterozygous
#'          \item \code{-1}: homozygous for major allele
#'          \item \code{NA}: missing data
#'          \item Imputed genotypes can be passed, see \code{impute} below for details
#'          }
#'          TIP - Set header=\code{FALSE} within \code{\link{read.table}} or \code{\link{read.csv}} when importing a tab-delimited file containing data for \code{G.in}.
#' @param y.in \code{Matrix} of phenotypic data. First column contains entry (taxa) names found in \code{G.in}, regardless of whether the entry has a phenotype for any or all traits. Additional columns contain phenotypic data; column names should reflect the trait name(s). TIP - Set header=\code{TRUE} within \code{read.table} or \code{read.csv} when importing a tab-delimited file contianing dat
#' @param impute Options include \code{c("EM", "mean", "pass")}. By default (i.e. \code{"EM"}), after filtering missing genotypic data will be imputed via the EM algorithm implemented in \code{\link{rrBLUP}} (\cite{Endelman, 2011}; \cite{Poland et al., 2012}). If \code{"mean"} missing genotypic data will be imputed via the 'marker mean' method, also implemented in \code{\link{rrBLUP}}. Enter \code{"pass"} if a pre-filtered and imputed genotype matrix is provided to \code{G.in}.
#' @param min.maf Optional \code{numeric} indicating a minimum minor allele frequency (MAF) when filtering \code{G.in}. Markers with an MAF < \code{min.maf} will be removed. Default is \code{0.01} to remove monomorphic markers. Set to \code{0} for no filtering.
#' @param mkr.cutoff Optional \code{numeric} indicating the maximum missing data per marker when filtering \code{G.in}. Markers missing > \code{mkr.cutoff} data will be removed. Default is \code{0.50}. Set to \code{1} for no filtering.
#' @param entry.cutoff Optional \code{numeric} indicating the maximum missing genotypic data per entry alloed when filtering \code{G.in}. Entries missing > \code{entry.cutoff} marker data will be removed. Default is \code{0.50}. Set to \code{1} for no filtering.
#' @param remove.dups Optional \code{logical}. If \code{TRUE} duplicate entries in the genotype matrix, if present, will be removed. This step may be necessary for missing marker imputation (see \code{impute}). Default is \code{TRUE}.
#' @param frac.train Optional \code{numeric} indicating the fraction of the TP that is used to estimate marker effects (i.e. the prediction set) under cross-validation (CV) method 1 (see \code{Details}). The remaining \eqn{(1-frac.trait)} of the TP will then comprise the prediction set.
#' @param nCV.iter Optional \code{integer} indicating the number of times to iterate \emph{CV method 1} described in \code{Details}. Default is \code{100}.
#' @param nFold Optional \code{integer}. If a number is provided, denoting the number of "folds", then CV will be conducted using \emph{CV method 2} (see \code{Details}). Default is \code{NULL}, resulting in the default use of the \emph{CV method 1}.
#' @param nFold.reps Optional \code{integer} indicating the number of times \emph{CV method 2} is repeated. The CV accuracy returned is the average \emph{r} of each rep. Default is \code{1}.
#' @param return.estimates Optional \code{logical}. If \code{TRUE} additional items including the marker effect and beta estimates from the selected prediction model (i.e. highest CV accuracy) will be returned.
#' @param CV.burnIn Optional \code{integer} argument used by \code{\link[BGLR]{BGLR}} when fitting Bayesian models. Default is \code{750}.
#' @param CV.nIter  Optional \code{integer} argument used by \code{\link[BGLR]{BGLR}} (\cite{de los Compos and Rodriguez, 2014}) when fitting Bayesian models. Default is \code{1500}.
#' @param models Optional \code{character vector} of the regression models to be used in CV and to estimate marker effects. Options include \code{rrBLUP, BayesA, BayesB, BayesC, BL, BRR}, one or more may be included at a time. By default all models are tested.
#' @return A list containing: \itemize{
#'            \item \code{CVs} A \code{dataframe} of CV results for each trait/model combination specified
#'            \item If \code{return.estimates} is \code{TRUE} the additional items will be returned: \itemize{
#'                \item \code{models.used} A \code{list} of the models chosen to estimate marker effects for each trait
#'                \item \code{mkr.effects} A \code{vector} of marker effect estimates for each trait generated by the respective prediction model used
#'                \item \code{betas} A \code{list} of beta values for each trait generated by the respective prediction model used
#'            }  
#'          }
#' @details Two CV methods are available within \code{PopVar}: \itemize{
#'            \item \code{CV method 1}: During each iteration a training (i.e. model training) set will be \strong{randomly sampled} from the TP of size \eqn{N*(frac.train)}, where \emph{N} is the size of the TP, and the remainder of the TP is assigned to the validation set. The accuracies of individual models are expressed as average Pearson's correlation coefficient (\emph{r}) between the genome estimated breeding value (GEBV) and observed phenotypic values in the validation set across all \code{nCV.iter} iterations. Due to its amendibility to various TP sizes, \emph{CV method 1} is the default CV method in \code{\link{pop.predict}}.
#'            \item \code{CV method 2}: \code{nFold} \strong{independent} validation sets are sampled from the TP and predicted by the remainder. For example, if \eqn{nFold = 10} the TP will be split into 10 equal sets, each containing \eqn{1/10}-th of the TP, which will be predicted by the remaining \eqn{9/10}-ths of the TP. The accuracies of individual models are expressed as the average (\emph{r}) between the GEBV and observed phenotypic values in the validation set across all \code{nFold} folds. The process can be repeated \code{nFold.reps} times with \code{nFold} new independent sets being sampled each replication, in which case the reported prediction accuracies are averages across all folds and replications.
#'          }
#' @examples
#' \dontrun{
#' ## CV using method 1 with 25 iterations
#' CV.mthd1 <- x.val(G.in = G.in_ex, y.in = y.in_ex, nCV.iter = 25)
#' CV.mthd1$CVs
#' 
#' ## CV using method 2 with 5 folds and 3 replications
#' x.val(G.in = G.in_ex, y.in = y.in_ex, nFold = 5, nFold.reps = 3)
#' }
#' @export

x.val <- function(G.in=NULL, y.in=NULL, min.maf=0.01, mkr.cutoff=0.50, entry.cutoff=0.50, remove.dups=T, impute="EM", frac.train=0.60, nCV.iter=100, nFold=NULL, nFold.reps=1, return.estimates=F, CV.burnIn=750, CV.nIter=1500, models=c("rrBLUP", "BayesA", "BayesB","BayesC", "BL", "BRR")){
  
  ## QC steps
  if(is.null(G.in)) stop("Must provide a genotype (G.in) file.")
  if(is.null(y.in)) stop("Must provide a phenotype (y.in) file.")
  if(!is.null(min.maf) & min.maf >= 1) stop("min.maf must be within the range [0, 1)")
  if(!is.null(entry.cutoff) & entry.cutoff > 1) stop("entry.cutoff must be within the range (0, 1]")
  if(!is.null(mkr.cutoff) & mkr.cutoff > 1) stop("mkr.cutoff must be within the range (0, 1]")
  if(impute == "pass"){min.maf=0; mkr.cutoff=1; entry.cutoff=1}
  
  ### Requird functions found in 'Internal_PopVar_functions_2.20.15.R'  
  
  ###### START HERE ##############
  ## Step 1 - Parse out Geno and Map files
  G.entries <- as.character(G.in[-1, 1])
  entries.removed <- NULL; entries.to.remove <- c() ## This is needed for output, may be replaced by list of entries if filtering is enabled or if duplicate entries found
  G.markers <- as.character(t(G.in[1, -1]))
  mkrs.removed <- NULL; mkrs.to.remove <- c()  ## This is needed for output, may be replaced by list of markers if filtering is enabled
  
  ## Marker, map, and geno QC
  G.mat <- as.matrix(G.in[-1, -1]); class(G.mat) <- "numeric"
  if(impute != "pass" && !all(unique(G.mat[,1]) %in% c(-1, 0, 1, NA))) stop("\nNon-imputed genotypes need to be coded as -1, 0, 1.\nIf imputed genotypes are passed set impute = 'pass'")
  
  ### ### Filter markers first for missing data and MAF
  mkrs.to.remove <- c()
  if(min.maf > 0){maf.list <- apply(G.mat, 2, maf.filt); mkrs.to.remove <- c(mkrs.to.remove, which(maf.list < min.maf))}
  if(mkr.cutoff <1){mkrNA.list <- apply(G.mat, 2, function(M){return(length(which(is.na(M))) / length(M))})
                    mkrs.to.remove <- unique(c(mkrs.to.remove, which(mkrNA.list > mkr.cutoff)))}
  
  if(length(mkrs.to.remove > 0)){
    G.mat <- G.mat[, -mkrs.to.remove]
    mkrs.removed <- G.markers[mkrs.to.remove]
    G.markers <- G.markers[-mkrs.to.remove]
  }
  
  ### Filter for duplicated entries and missing data
  entries.to.remove <- c()
  if(remove.dups) entries.to.remove <- c(entries.to.remove, which(duplicated.array(G.mat)))
  if(entry.cutoff < 1){entryNA.list <- apply(G.mat, 1, function(E){return(length(which(is.na(E))) / length(E))})
  entries.to.remove <- unique(c(entries.to.remove, which(entryNA.list > entry.cutoff)))}
  
  if(length(entries.to.remove > 0)){
    G.mat <- G.mat[-entries.to.remove, ]
    entries.removed <- G.entries[entries.to.remove]
    G.entries <- G.entries[-entries.to.remove]
  }
  
  y <- y.in[match(G.entries, as.character(y.in[,1])),]
  y.entries <- as.character(y[,1])
  traits <- as.character(colnames(y))[-1]; nTraits <- length(traits)
  
  
  ## Imput missing markers with EM... will switch to imputing with the mean if nEntries > nMarkers
  ## Will need to use our own MAF filter so that we can keep track of which markers are removed due to MAF and missing data
  if(impute == "EM") G.imp <- rrBLUP::A.mat(G.mat, min.MAF = 0, max.missing = 1, impute.method = "EM", return.imputed = T)$imputed
  if(impute == "mean") G.imp <- rrBLUP::A.mat(G.mat, min.MAF = 0, max.missing = 1, impute.method = "mean", return.imputed = T)$imputed
  if(impute == "pass") G.imp <- G.mat
  
  ### Start simulation
  for(t in 1:nTraits){
    
    trait <- traits[t]
    
    ### the y.noNA and G.noNA define the TP and are used for CV and marker effect estimation
    y_notNAs <- !is.na(y[,trait])
    y_TP <- as.numeric(y[y_notNAs, trait])
    TP.entries <- y.entries[y_notNAs]
    G_TP <- G.imp[y_notNAs, ] 
    
    cat("\n")
    cat(paste("\nPerforming cross validation for ", trait, sep=""))
    
    if(is.null(nFold)) junk <- capture.output(xval.out <- XValidate.nonInd(y.CV = y_TP, G.CV = G_TP, models.CV = models, frac.train.CV=frac.train, nCV.iter.CV=nCV.iter, burnIn.CV = CV.burnIn, nIter.CV = CV.nIter)$CV.summary)
    if(!is.null(nFold)) junk <- capture.output(xval.out <- XValidate.Ind(y.CV = y_TP, G.CV = G_TP, models.CV = models, nFold.CV = nFold, nFold.CV.reps = nFold.reps, burnIn.CV = CV.burnIn, nIter.CV = CV.nIter)$CV.summary)
    
    if(t == 1){
      CV.results <- list()
      mkr_effs.mat <- matrix(NA, ncol = nTraits, nrow = ncol(G.imp))
      colnames(mkr_effs.mat) <- traits
      beta.list <- c()
      best.models <- c()
    }
    
    CV.results[[t]] <- xval.out
    names(CV.results)[[t]] <- trait
    
    best.model <- as.character(xval.out$Model[which(xval.out$r_avg == max(xval.out$r_avg))])
    best.models[t] <- best.model; names(best.models)[t] <- trait
    
    if(best.model == "rrBLUP"){
      mix.solve.out <- rrBLUP::mixed.solve(y_TP, Z=G_TP, SE=F, return.Hinv=F)
      beta <- as.numeric(mix.solve.out$beta)
      mkr_effects <- mix.solve.out$u
    }else{
      capture.output(bayes.fit <- BGLR::BGLR(y=y_TP, ETA=list(list(X=G_TP, model=best.model)), verbose=F, nIter=CV.nIter, burnIn=CV.burnIn))
      mkr_effects <- as.numeric(bayes.fit$ETA[[1]]$b)
      beta <- bayes.fit$mu  
    }
    
    beta.list[t] <- beta
    mkr_effs.mat[,t] <- mkr_effects
    
  }
  
  if(return.estimates) return(list(CVs=CV.results, models.used=best.models, mkr.effects=mkr_effs.mat, betas=beta.list))
  if(!return.estimates) return(list(CVs=CV.results))
  
} # End of PopVar.CV

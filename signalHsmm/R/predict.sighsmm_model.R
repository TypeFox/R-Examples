#' Predict sighsmm_model object
#' 
#' Predicts the presence of signal peptides using signalHsmm models.
#' @param object \code{sighsmm_model} object.
#' @param newdata unknown sequence of class \code{character} or \code{character}.
#' Alternatively, a \code{list} of sequences in mentioned formats.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @examples
#' #remember to remove it
#' \dontrun{
#' pos_train_ultrahard <- read_uniprot("pos_ultrahard_data.txt", euk = TRUE)
#' model1 <- train_hsmm(pos_train_ultrahard, aa_group = aaaggregation)
#' predict(model1, benchmark_dat[1L:5])
#' }



predict.sighsmm_model <- function(object, newdata, ...){
  if (class(newdata) == "numeric" || class(newdata) == "factor" || 
        class(newdata) == "data.frame" || class(newdata) == "matrix")
    stop("Input data must have class 'SeqFastaAA', 'character' or 'list'.")
  
  if(class(newdata) == "SeqFastaAA" || 
       class(newdata) == "character") {
    #single input
    decisions <- signalHsmm_decision(newdata, aa_group = object[["aa_group"]], 
                                      pipar = object[["pipar"]], 
                                      tpmpar = object[["tpmpar"]], 
                                      od = object[["od"]], 
                                      overall_probs_log = object[["overall_probs_log"]], 
                                      params = object[["params"]])
    decisions <- list(decisions)
    names(decisions) <- attr(newdata, "name")
  } else {
    #list input
    decisions <- lapply(newdata, function(prot)
      signalHsmm_decision(prot, aa_group = object[["aa_group"]], 
                           pipar = object[["pipar"]], 
                           tpmpar = object[["tpmpar"]], 
                           od = object[["od"]], 
                           overall_probs_log = object[["overall_probs_log"]], 
                           params = object[["params"]]))
  }
  class(decisions) <- "hsmm_pred_list"
  decisions
}

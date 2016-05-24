# SIGNAL-HSMM ------------------------------------

#' Predict presence of signal peptide in protein
#'
#' Using the hidden semi-Markov model predict presence of signal peptide in 
#' eukaryotic proteins.
#'
#' @param test_data single protein sequence (\code{character} vector) or list of 
#' sequences. It may be an object of class \code{\link[seqinr]{SeqFastaAA}}.
#' @return An object of class \code{hsmm_pred_list}.
#' @details Function \code{signalHsmm} returns respectively probability of presence of 
#' signal peptide, start of signal peptide and the probable cleavage site localization.
#' If input consists of more than one sequence, result is a data.frame where each column
#' contains above values for different proteins.
#' @note Currently start of signal peptide is naively set as 1 amino acid. The prediction 
#' of a cleavage site is still an experimental feature, use on your own risk.
#' @useDynLib signalHsmm
#' @export
#' @seealso \code{\link{hsmm_pred_list}} \code{\link{hsmm_pred}} 
#' @keywords classif
#' @examples
#' #run signalHsmm on one sequence
#' x1 <- run_signalHsmm(benchmark_dat[[1]])
#' 
#' #run signalHsmm on one sequence, but input is a character vector
#' x2 <- run_signalHsmm(c("M", "A", "G", "K", "E", "V", "I", "F", "I", "M", "A", "L", 
#' "F", "I", "A", "V", "E", "S", "S", "P", "I", "F", "S", "F", "D", 
#' "D", "L", "V", "C", "P", "S", "V", "T", "S", "L", "R", "V", "N", 
#' "V", "E", "K", "N", "E", "C", "S", "T", "K", "K", "D", "C", "G", 
#' "R", "N", "L", "C", "C", "E", "N", "Q", "N", "K", "I", "N", "V", 
#' "C", "V", "G", "G", "I", "M", "P", "L", "P", "K", "P", "N", "L", 
#' "D", "V", "N", "N", "I", "G", "G", "A", "V", "S", "E", "S", "V", 
#' "K", "Q", "K", "R", "E", "T", "A", "E", "S", "L"))
#' 
#' #run signalHsmm on list of sequences
#' x3 <- run_signalHsmm(benchmark_dat[1:3])
#' #see summary of results
#' summary(x3)
#' #print results as data frame
#' pred2df(x3)
#' #summary one result
#' summary(x3[[1]])
#' plot(x3[[1]])

run_signalHsmm <- function(test_data) {
  predict.sighsmm_model(signalHsmm_main_model, test_data)
}

signalHsmm_decision <- function(prot, aa_group, pipar, tpmpar, 
                                od, overall_probs_log, params) {
  if (length(prot) == 1) {
    prot <- strsplit(prot, "")[[1]]
    if ("name" %in% names(attributes(prot)))
      attr(prot, "name") <- "undefined_name"
    if (length(prot) == 1)
      stop("Input sequence is too short.")
  }
  if(!is_protein(prot))
    stop("Atypical aminoacids detected, analysis cannot be performed.")
  
  deg_sample <- as.numeric(degenerate(prot[1L:ifelse(length(prot) > 50, 50, length(prot))], aa_group))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi(deg_sample - 1, pipar, tpmpar, od, params)
  viterbi_path <- viterbi_res[["path"]] + 1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)) + 1, 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  res <- list(sp_probability = rescale(unname(1 - 1/(1 + prob.total))),
              sp_start = 1,
              sp_end = c_site,
              struc = viterbi_path,
              prot = toupper(prot[1L:70]),
              name = attr(prot, "name"),
              str_approx = 0)
  class(res) <- "hsmm_pred"
  
  #structure approximation - if atypical (normally negative signal peptide)
  while(!all(1L:4 %in% res[["struc"]])) {
    res[["struc"]] <- c(res[["struc"]], which.min(1L:4 %in% res[["struc"]]))
    res[["str_approx"]] <- res[["str_approx"]] + 1
  }
  
  res
}


#' GUI for signalHsmm
#'
#' A graphical user interface for predicting presence of signal peptides.
#' @return null.
#' @export
#' @seealso \code{\link{run_signalHsmm}}
#' @note
#' Any ad-blocking software may be cause of malfunctions.

gui_signalHsmm <- function() {
  runApp(system.file("signal_gui", package = "signalHsmm"))
}

signalHsmm_main_model <- structure(list(aa_group = structure(list(`1` = c("R", "N", "D", "Q", "E", "H", "K"), 
                                                                  `2` = c("G", "P", "S", "T", "Y"), 
                                                                  `3` = c("I", "L", "M", "F", "W", "V"), 
                                                                  `4` = c("A", "C")), .Names = c("1", "2", "3", "4")), 
                                        pipar = c(1, 0, 0, 0), 
                                        tpmpar = structure(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0), .Dim = c(4L, 4L)), 
                                        od = structure(c(0.282182230869001, 0.0128884115403133, 0.20406633561323, 
                                                         0.32847743971351, 0.288424124513619, 0.119143735417207, 0.413098472781786, 
                                                         0.310808047593652, 0.336251621271077, 0.702603607273805, 0.233673756207252, 
                                                         0.260085051824138, 0.0931420233463035, 0.165364245768675, 0.149161435397733, 
                                                         0.1006294608687), .Dim = c(4L, 4L)), 
                                        overall_probs_log = structure(c(-1.11328712008703, -1.1685797676654, -1.34674657905895, -2.29631021261161), 
                                                                      .Names = c("1",  "2", "3", "4")), 
                                        params = structure(c(0.0946771378708551, 0.18804537521815, 
                                                             0.191972076788831, 0.0824607329842932, 0.0977312390924956, 0.087260034904014, 
                                                             0.0536649214659686, 0.043630017452007, 0.0344677137870855, 0.0253054101221641, 
                                                             0.0218150087260035, 0.018760907504363, 0.0117801047120419, 0.0104712041884817, 
                                                             0.0109075043630017, 0.00610820244328098, 0.00567190226876091, 
                                                             0.00392670157068063, 0.00305410122164049, 0.00261780104712042, 
                                                             0.00043630017452007, 0.00087260034904014, 0.00174520069808028, 
                                                             0.00261780104712042, 0, 0, 0, 0, 0, 0, 0, 0, 0.000864677907479464, 
                                                             0.000864677907479464, 0.00302637267617812, 0.0012970168612192, 
                                                             0.00345871162991786, 0.0142671854734112, 0.0380458279290964, 
                                                             0.0626891482922611, 0.0912235192390834, 0.126675313445741, 0.155642023346304, 
                                                             0.137916126242975, 0.131863380890618, 0.0860354517942067, 0.0514483354950281, 
                                                             0.041936878512754, 0.0229139645482058, 0.0116731517509728, 0.00778210116731518, 
                                                             0.00475572849113705, 0.00172935581495893, 0.000864677907479464, 
                                                             0.00172935581495893, 0.000432338953739732, 0.000864677907479464, 
                                                             0, 0, 0, 0, 0, 0, 0, 0, 0.142237640936687, 0.146140503035559, 
                                                             0.293148308759757, 0.108846487424111, 0.167389418907199, 0.0633130962705984, 
                                                             0.0450997398091934, 0.0177797051170859, 0.0108412836079792, 0.00216825672159584, 
                                                             0.00303555941023417, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                             0, 0, 0, 0, 0, 0, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 
                                                             0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125), .Dim = c(32L, 4L), 
                                                           .Dimnames = list(NULL, c("n", "h", "c", "")))), 
                                   .Names = c("aa_group", "pipar", "tpmpar", "od", "overall_probs_log", "params"), class = "sighsmm_model")
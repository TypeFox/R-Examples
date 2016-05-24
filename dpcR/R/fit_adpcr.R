fit_adpcr <- function(pcr_data, cyc = 1, fluo = NULL, model = l5, norm = FALSE, 
                      iter_tr = 50){
  if (class(pcr_data) == "list") {
    lapply(pcr_data, function (i) fit_single_adpcr(i, cyc, fluo, model, norm, iter_tr))
  } else {
    fit_single_adpcr(pcr_data, cyc, fluo, model, norm, iter_tr)
  }
}

#not for users
fit_single_adpcr <- function(pcr_data, cyc, fluo, model, norm, iter_tr) {
  if (is.null(fluo)) {
    all_fluos <- (1L:ncol(pcr_data))[-cyc]
  } else {
    all_fluos <- fluo
  }
  anal_fluo <- all_fluos[vapply(all_fluos, function(x) 
    valid_amp(pcr_data[[x]]), TRUE)]
  
  #only validated columns are fitted
  all_fits <- modlist(pcr_data, cyc, anal_fluo, model = model, norm = norm, 
                      verbose=FALSE, remove = "KOD")
  
  #temporary solution, should change nls.lm.control instead
  good_fits_ind <- vapply(all_fits, function(x) 
    x[["convInfo"]][["finIter"]], 0) < iter_tr
  
  good_fits <- all_fits[(1L:length(all_fits))[good_fits_ind]]
  class(good_fits) <- c("modlist", "pcrfit")
  names(good_fits) <- vapply(good_fits, function(x) x[["names"]], "a")
  good_fits
}
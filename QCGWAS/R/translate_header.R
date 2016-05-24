translate_header <-
function(header, standard = c("MARKER","CHR","POSITION","EFFECT_ALL","OTHER_ALL","STRAND","EFFECT","STDERR","PVALUE","EFF_ALL_FREQ","HWE_PVAL","CALLRATE","N_TOTAL","IMPUTED", "USED_FOR_IMP", "IMP_QUALITY"), alternative) {
  if(any(duplicated(alternative[ ,2]))) stop("duplicated elements in alternative, column 2")
  capitalized <- toupper(header)
  unknowns	<- !logical(length = length(header))
  missings	<- logical(length = length(standard))
  for(forI in 1:length(standard)) {
    column_no <- identify_column(std_name = standard[forI], alt_names = alternative, header = capitalized)
    if(length(column_no) == 0L) {
      missings[forI] <- TRUE
    } else {
      header[column_no] <- standard[forI]
      unknowns[column_no] <- FALSE
    }
  }	
  return(list(header_N = length(header), header_h = header,
              missing_N = sum(missings), missing_h = if(sum(missings) == 0L) NULL else standard[missings],
              unknown_N = sum(unknowns), unknown_h = if(sum(unknowns) == 0L) NULL else header[unknowns]	) )
}

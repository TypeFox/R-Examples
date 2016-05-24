#' Mixed effect model formula parser
#'  
#' Parses a mixed effect model in the lme4 structure of  ~fixed +(re1|g1) +...+(ren|gn)
#' 
#' @param f formula for mixed effect mode in the form used in lme4; ~fixed +(re1|g1) +...+(ren|gn)
#' @return A list with elements fix.model and re.model. fix.model contains the formula for the fixed effects;
#' re.model contains elements sub, the grouping formula and model the design formula for the
#' random effect. Each formula is of type character and must be wrapped with as.formula in use with model.matrix
#' @author Devin Johnson <devin.johnson@@noaa.gov>
#' @import lme4
proc.form <- function(f){
  fix.model = paste(deparse(lme4::nobars(f),width.cutoff = 500L),collapse="")
  if(fix.model=="NULL") fix.model="~ 1"
  re.lst = lme4::findbars(f)
	if(length(re.lst)==0){
    re.model <- NULL
	} else{
		re.model <- lapply(re.lst, reSplit)
		names(re.model) <- sapply(re.lst, function(x){deparse(x)})
	}
	return(list(fix.model=fix.model, re.model=re.model))
}
reSplit = function(x){
#  s=strsplit(deparse(x), " | ")[[1]][-2]
  s=strsplit(deparse(x), "|",fixed=TRUE)[[1]]
  model=paste("~", s[1])
  sub=paste("~", s[2], "- 1")
  return(list(model=model, sub=sub))
}
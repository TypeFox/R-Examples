freqCI <- function(x, level = .95){
  if(missing(x)) stop('"x" must be specified')
  this_call <- match.call()

  # x to res (working object) with checks
  if(class(x) == "table"){
    res <- x
    if(length(dim(res)) != 1L) stop('if "x" is a table, it must not have more than one dimension.')
    if(!all((res %% 1) == 0)) warning('"x" is a table and contains non-integer values. Check if this is intended.', call. = FALSE, immediate. = TRUE)
  } else if(is.null(dim(x))){
    if(is.numeric(x) && length(x) < 10L) warning('"x" looks like it is already a table.\nUse as.table(x) in this case.', call. = FALSE, immediate. = TRUE)
    if(is.character(x)) warning('"x" was supplied as a "character" vector. A "factor"-type variable was assumed.', call. = FALSE, immediate. = TRUE)
    res <- table(x)
  } else {
    stop('"x" must either be a vector of individual observations or a "table" object with frequencies.')
  }

  # error checks level
  if(!is.null(dim(level)) || !is.numeric(level)) stop('"level" must be a number or a vector of numbers.')
  if(any(level <= 0.0) || any(level >= 1.0)) stop('CI level(s) must be in (0, 1).')

  level <- sort(level)
  
  n <- sum(res)
  rel_freq <- as.numeric(res / n)
  const <- lapply(level, function(CI_lev){
    abs(qnorm((1-CI_lev)/2)) * sqrt(rel_freq*(1-rel_freq)/n)
  })

  res <- structure(
    list("call"      = this_call,
         "x"         = x,
         "level"     = level,
         "freq"      = unclass(res),
         "n"         = n,
         "rel_freq"  = rel_freq,
         "cat_names" = names(res),
         "CIs_low"   = matrix(unlist(lapply(rev(const), function(ci_cst){ rel_freq - ci_cst })), nrow = length(res)),
         "CIs_high"  = matrix(unlist(lapply(const, function(ci_cst){ rel_freq + ci_cst })), nrow = length(res))),
    "class" = "freqCI"
  )

  return(res)

}

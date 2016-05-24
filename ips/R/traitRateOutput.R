traitRateOutput <- function(file, what = "likelihoods"){
  
  what <- match.arg(what, c("rates", "likelihoods"))
  sep <- ifelse(what == "rates", "", "\n")
  res <- scan(file, what = "c", sep = sep, 
              quiet = TRUE)
  if ( what == "rates" ){
    res <- res[grep("charModelParam", res)]
    res <- as.numeric(gsub("^.+=", "", res))
    names(res) <- c("r01", "r10")
  }
  if ( what == "likelihoods" ){
    id <- grep("LogLikelihood Model 0", res)
    res <- as.numeric(gsub("^.+= ", "", res[-1:0 + id]))
    names(res) <- c("logL", "logL0")
  }
  res
}
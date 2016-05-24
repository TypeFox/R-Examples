## check specification of SE and SP

checkBinPrior <-
function(x, name){
  if (class(x) == "numeric"){
    out <- checkSeSp(list(dist = "fixed", par = x))

  } else if (class(x) == "formula"){
    ## 'x' should be of length 2 ('~' + dist)
    if (length(x) != 2)
      stop("Formula specification of ", name, " is incorrect.\n",
           "See ?truePrev for more details.")

    call <- as.character(x)[[2]]
    dist2list(call, type = "prob")

  } else if (class(x) == "list"){
    check <- checkSeSp(x)

  } else {
    stop(name, " should be specified as a list or a formula.\n",
         "See ?truePrev for more details.")
  }
}
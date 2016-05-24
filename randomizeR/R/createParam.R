#' @include randPar.R
#' @include randSeq.R
NULL

#' Representing any randomization procedure 
#' 
#' Represents any randomization procedure for a two-armed clinical trial.
#' 
#' @details Dending on the input of the user, \code{createParam} creates an object
#' representing a randomization procedures for a two-armed clinical trial
#' (see also \code{\link{randPar}}).
#' 
#' @family randomization procedures
#'
#' @inheritParams overview
#'
#' @return
#' S4object of the corresponding randomization procedure class.
#' 
#' @export
createParam <- function(method, N, mti, bc, rb, p, ini, add) {
  Cls <-  names(getClass("randPar")@subclasses)
  RPs <- sub("Par", "", Cls)
  ## query for the method, if not specified!
  repeat{
    if(missing("method") || !(method %in% toupper(RPs))){
      method <- readline( 
        cat(paste("Which randomization method do you want to use?\n(Possible values are: ",
                  paste(toupper(RPs), sep = "", collapse = ", "), ") \n>", sep = "")))
    } else {
      break
    }
  }
  check <- TRUE
  dec <- numeric(0)
  # choosing the right randomization method
  rpPar <- paste(tolower(method), "Par", sep = "")
  # defined slotnames for the chosen randomization method
  slotns <- slotNames(rpPar)
  # queries for the slotnames parameters
  robEval <- function(x, y) {
    tryCatch(eval(parse(text = (paste(x, " <- ", y, sep = "")))),
             error = function(e) {"error"})    
    tryCatch(eval(parse(text = (paste(x, " <<- ", y, sep = "")))),
             error = function(e) {"error"})
    tryCatch(eval(parse(text = (paste("is.numeric(", x, ")", sep = "")))),
             error = function(e) {FALSE})
  }
  slotns <- slotns[!(slotns %in% c("K", "ratio", "groups"))]
  for(i in 1:length(slotns)){
    eval(parse(text = paste("dec <- missing(", slotns[i], ")", sep = "\"")))
    if(dec){
      repeat{
        # checking if the value was ok
        param <- readline(cat(paste("Value for the paramter ",
                                    slotns[i], ": \n>", sep = "")))
        test <- robEval(slotns[i], param)
        if(!test)
          next
                
        c1 <- paste(slotns[i], "=", slotns[i], sep = "")
        c2 <- paste("check <- paramErrors(", c1, ")", sep = "")
        eval(parse(text = c2))
        
        if(!check) break
      }  
    }   
  }

  constr <- numeric(0)
  ## Output:
  output1 <- paste(tolower(method), "Par(", sep = "")
  output2 <- paste(slotns, "=", slotns, sep = "", collapse = ",")
  output3 <- paste("constr <- ", output1, output2,")", sep = "")
  eval(parse(text = output3))

  constr  
}



#' Query to create a randomization sequence of a particular randomization procedure 
#' 
#' This function is a query to create an corresponding randomization sequence
#' for a two-armed clinical trial. If
#' \code{file} is defined, the generated sequence is automatically saved to the
#' corresponding path.
#' 
#' @inheritParams overview
#'
#' @return an object \code{Param}, which is available 
#' 
#' @export
createSeq <- function(file) {
  constr <- createParam()
  Seq <- genSeq(constr)
  if(!missing(file))
    saveRand(Seq, file)
  return(Seq)
}  



# Function for errors requesting 
# 
# This function is a query to make sure that the parameters are all in the right range.
# 
# @inheritParams overview
#
# @return returns a \code{TRUE} if everything is fine, otherwise a \code{FALSE}
paramErrors <- function(method, N, mti, bc, rb, p, ini, add) {
  Cls <-  names(getClass("randPar")@subclasses)
  RPs <- sub("Par", "", Cls)
  out <- FALSE
  # error request for the method
  if(!missing("method"))
    if(!(method %in% toupper(RPs)))
      out <- TRUE
  # error request for N
  if(!missing("N"))   
    if(!(length(N) == 1 && is.numeric(N) && N > 0 && (N %% 2 == 0)))
      out <- TRUE
  # error request for mti
  if(!missing("mti"))
    if(!(length(mti) == 1 && is.numeric(mti) && mti == ceiling(mti)))
      out <- TRUE
  # error request for bc
  if(!missing("bc"))
    if(!(all(is.numeric(bc)) && all(bc > 0) && all(bc %% 2 == 0)))
      out <- TRUE
  # error request for rb
  if(!missing("rb"))
    if(!(all(is.numeric(rb)) && all(rb > 0) && all(rb%%2 == 0) ))
      out <- TRUE  
  # error request for p
  if(!missing("p"))
    if(!(length(p) == 1 && is.numeric(p) && p >= 0.5 && p <= 1))
      out <- TRUE
  # error request for ini
  if(!missing("ini"))
    if(!(length(ini) == 1 && is.numeric(ini) && ini == round(ini)))
      out <- TRUE
  # error request for add
  if(!missing("add"))
    if(!(length(add) == 1 && is.numeric(add) && add == round(add)))
      out <- TRUE
  return(out)
}

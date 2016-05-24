writeModel <- function(model, con = "model.txt", digits = 5)
{
  if (is.R()){
    model.text <- c("model", replaceScientificNotationR(body(model), digits = digits))
    # "[\+\-]?\d*\.?[Ee]?[\+\-]?\d*"
  } else {
    ## In S-PLUS the source code of a function can be obtained with
    ## as.character(function_name).  This omits the "function_name <- function()" piece
    model.text <- paste("model", as.character(model))
  }
  model.text <- gsub("%_%", "", model.text)
  if (!is.R()){
    ## In S-PLUS, scientific notation is different than it is in WinBUGS.
    ## Change the format of any numbers in scientific notation.
    model.text <- replaceScientificNotationS(model.text)

    ## remove the "invisible()" line.  
    model.text <- gsub("invisible[ ]*\\([ ]*\\)", "", model.text)
  }
  writeLines(model.text, con = con)
}


replaceScientificNotationR <- function(bmodel, digits = 5){
    env <- new.env()
    assign("rSNRidCounter", 0, envir=env)
    replaceID <- function(bmodel, env, digits = 5){
        for(i in seq_along(bmodel)){
            if(length(bmodel[[i]]) == 1){
                if(as.character(bmodel[[i]]) %in% c(":", "[", "[[")) return(bmodel)
                if((typeof(bmodel[[i]]) %in% c("double", "integer")) && ((abs(bmodel[[i]]) < 1e-3) || (abs(bmodel[[i]]) > 1e+4))){
                    counter <- get("rSNRidCounter", envir=env) + 1
                    assign("rSNRidCounter", counter, envir=env)
                    id <- paste("rSNRid", counter, sep="")
                    assign(id, formatC(bmodel[[i]], digits=digits, format="E"), envir=env)
                    bmodel[[i]] <- id
                }
            } else {
                bmodel[[i]] <- replaceID(bmodel[[i]], env, digits = digits)
            }
        }
        bmodel
    }
    bmodel <- deparse(replaceID(bmodel, env, digits = digits), control = NULL)
    for(i in ls(env)){
        bmodel <- gsub(paste('"', i, '"', sep=''), get(i, envir=env), bmodel, fixed=TRUE)
    }
    bmodel
}



replaceScientificNotationS <- function(text){
## Change the format of any numbers in "text" that are in S-PLUS 
## scientific notation to WinBUGS scientific notation

  ## First, handle the positive exponents
  ## Find the first instance
  ## Note that the number may or may not have a decimal point.
  sciNoteLoc <- regexpr("[0-9]*\\.{0,1}[0-9]*e\\+0[0-9]{2}", text)
    
  ## For every instance, replace the number
  while(sciNoteLoc > -1){
    sciNoteEnd <- sciNoteLoc + attr(sciNoteLoc, "match.length")-1
    sciNote <- substring(text, sciNoteLoc, sciNoteEnd)
    text <- gsub(sciNote, toSingleS4(sciNote), text)
    sciNoteLoc <- regexpr("[0-9]*\\.{0,1}[0-9]*e\\+0[0-9]{2}", text)
  }

  ## Then, handle the negative exponents
  ## Find the first instance
  sciNoteLoc <- regexpr("[0-9]*\\.{0,1}[0-9]*e\\-0[0-9]{2}", text)

  ## For every instance, replace the number
  while(sciNoteLoc > -1){
    sciNoteEnd <- sciNoteLoc + attr(sciNoteLoc, "match.length")-1
    sciNote <- substring(text, sciNoteLoc, sciNoteEnd)
    text <- gsub(sciNote, toSingleS4(sciNote), text)
    sciNoteLoc <- regexpr("[0-9]*\\.{0,1}[0-9]*e\\-0[0-9]{2}", text)
  }

  text
}

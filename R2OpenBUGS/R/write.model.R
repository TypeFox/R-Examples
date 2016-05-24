write.model <- function(model, con = "model.bug", digits = 5)
{
  model.text <- c("model", replaceScientificNotationR(body(model), digits = digits))
  # "[\+\-]?\d*\.?[Ee]?[\+\-]?\d*"
  model.text <- gsub("%_%", "", model.text)
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

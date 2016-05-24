#' @include model-bbinchoice.R
#' @include model-blogit.R
#' @include model-bprobit.R
#' @include model-ologit.R
#' @include model-oprobit.R
#' @include model-mlogit.R

#library(jsonlite)


createJSONzeligchoice <- function(){

  z5blogit <- zblogit$new()
  z5blogit$toJSON()

  z5bprobit <- zbprobit$new()
  z5bprobit$toJSON()

  z5mlogit <- zmlogit$new()
  z5mlogit$toJSON()

  z5ologit <- zologit$new()
  z5ologit$toJSON()

  z5oprobit <- zoprobit$new()
  z5oprobit$toJSON()

  zeligchoicemodels <- list(zelig5choicemodels = list("blogit" = z5blogit$ljson,
                                                    "bprobit" = z5bprobit$ljson,
                                                    "mlogit" = z5mlogit$ljson,
                                                    "ologit" = z5ologit$ljson,
                                                    "oprobit" = z5oprobit$ljson))

  # cat(jsonlite::toJSON(zeligchoicemodels, pretty = TRUE),
  #     file = file.path("inst/JSON", "zelig5choicemodels.json"))

  cat(toJSON(zeligchoicemodels, pretty = TRUE), file = file.path("zelig5choicemodels.json"))
  file.rename(from = file.path("zelig5choicemodels.json"),
            to = file.path("inst", "JSON", "zelig5choicemodels.json"))
  file.remove(file.path("zelig5choicemodels.json"))

  return(TRUE)
}
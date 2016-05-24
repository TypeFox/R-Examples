
.experiment <- function(data, method, problem, params, outs, name, historic){
  
  tags <- .metaTags(title = name, context = name, alias = "experiment")
  
  result <- list( "data" = data,
                  "method" = method,
                  "problem" = problem,
                  "parameters" = params,
                  "outputs" = outs,
                  "configuration" = list(),
                  "tags" = tags,
                  "historic" = historic)
  
  class(result) <- c("experiment", "reportable")
  result
}


is.experiment <- function(x) {
  is(x, "experiment")
}

#' @export
toString.experiment <- function (x, ...) {
  d <- x[["data"]]
  
  # Print experiment name
  result <- paste0("#Experiment name: ", x$tags$title,"\n\n")
  
  # Print method list
  result <- paste0(result,
                  sprintf("#%s: %s\n", 
                          x$method, paste(levels(x$data[[x$method]]),
                                          collapse = ', ') ), "\n")
  # Print problem list
  result <- paste0(result,
                  sprintf("#%s: %s\n",
                          x$problem, paste(levels(x$data[[x$problem]]),
                                     collapse = ', ') ), "\n")
  
  result <- paste0(result,"#parameters:\n")
  
  # Print the parameters list if any
  params <- c()
  if (length(x$parameters) != 0) 
    for (p in x$parameters)
      params <- c(params, paste0(p, ' [', paste0(levels(x$data[[p]]), collapse = ","), ']'))
    
  if (length(x$configuration) != 0) 
    params <- c(params, x$configuration)
  
    result <- paste0(result, .nestedList2String(params, numbered=FALSE))
  
  result <- paste0(result, sprintf("\n#outputs: %s\n", paste(x[["outputs"]], collapse = ', ') ))
}

#' @export
print.experiment <- function(x, ...){
  cat( toString(x) )
}

#' @export
summary.experiment <- function (object, ...) {
  h <- object$historic
  
  cat(.nestedList2String(h))
  
  cat("\n")
  print(object)
}

names.experiment <- function(e){
  e$name
}
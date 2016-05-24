#' Run the Timeseries shiny app.
#'
#' @export
#' @examples
#' \dontrun{s_timeseries()}
s_timeseries<-function(){
  shiny::runApp(system.file("timeseries", package="learnstats"),quiet=TRUE)
}

#' Run the proportion confidence interval shiny app.
#'
#' @export
#' @examples
#' \dontrun{s_propci()}
s_propci<-function(){
  shiny::runApp(system.file("propconfint", package="learnstats"),quiet=TRUE)
}

#' Run the Two Normal Distributions shiny app.
#' @export
#' @examples
#' \dontrun{s_twonorm()}
s_twonorm<-function(){
    shiny::runApp(system.file("twonorm", package="learnstats"))
  }



#' Run the T distribution shiny app.
#' @export
#' @examples
#' \dontrun{s_tdist()}
s_tdist<-function(){
  shiny::runApp(system.file("tdist", package="learnstats"),quiet=TRUE)
}


#' Run the F distribution shiny app.
#' @export
#' @examples
#' \dontrun{s_fdist()}
s_fdist<-function(){
  shiny::runApp(system.file("fdist", package="learnstats"),quiet=TRUE)
}



#' Run the Normal Approximation to the Binomial shiny app.
#'  @export
#' @examples
#' \dontrun{s_normbinom()}
s_normbinom<-function(){
  shiny::runApp(system.file("normapprox2bin", package="learnstats"),quiet=TRUE)
}


#' Run the Intro Binomial shiny app.
#'  @export
#' @examples
#' \dontrun{s_introbinom()}
s_introbinom<-function(){
  shiny::runApp(system.file("introbinomial", package="learnstats"),quiet=TRUE)
}




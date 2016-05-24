# This is file ../spam/R/profile.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

".onLoad" <- function (lib, pkg) {
#    if (R.version$minor<paste(14))    require(methods) 

}

# Framework introduced with much input from Roger Bivand for 0.13-2 and higher.

spam.Version <- function() {
  release <- utils::packageDescription("spam",field="Version")
  date <- utils::packageDescription("spam",field="Date")
  list(status="",
              major=sub("-","",substr(release,1,4)),
              minor=substr(sub("-","",substr(release,5,7)),1,1),
              year=substr(date,1,4),
              month=substr(sub("200.-","",date),1,2),
              day=sub("200.-..-","",date),
              version.string= paste("Spam version ",
                utils::packageDescription("spam",field="Version")," (",
                utils::packageDescription("spam",field="Date"),")",sep='')
              )
}
                
spam.version <- spam.Version()
class(spam.version) <- "simple.list"
  
".Spam" <- list(eps=.Machine$double.eps,   # smaller than this is considered as zero
                drop=FALSE,                # drop passed to subset functions
                printsize=100,             # the max size which we print as regular matrices
                imagesize=10000,           # the max size which we display as regular matrices
                 cex=1200,                  # scaling factor for scatter displays

                structurebased=!FALSE,      # calculating on nonzero entries only...
                
                inefficiencywarning=1e6,  # tell when something inefficient is done
                
               trivalues=FALSE,           # with upper./lower/.tri return values (TRUE) or only structure?
                listmethod='PE',           # method to be used when using spam.list

                NAOK=FALSE,      #  passing of !is.finite to fortran
                safemodevalidity=TRUE,  # verify while S4 construction
                dopivoting=TRUE,           # what type of back/forwardsolve?
                cholsymmetrycheck=TRUE,     # Should symmetry be tested in the cholesky factorization
                cholpivotcheck=TRUE,        # Should the pivot be tested?
                cholupdatesingular="warning",     # ("error", "warning","NULL")
                cholincreasefactor=c(1.25,1.25),
                nearestdistincreasefactor=1.3,
                nearestdistnnz=c(500^2,500)
                )
#noquote(unlist(format(.Spam[-1])) )

"inefficiencywarning" <- function(msg,size) {
  maxsize <- if (is.logical(.Spam$inefficiencywarning)) {
    ifelse(.Spam$inefficiencywarning,1,Inf) } else { 
    .Spam$inefficiencywarning
  }
  if (size>maxsize) warning(msg, call. = FALSE)
}
    
".onAttach" <- function (lib, pkg) {
   packageStartupMessage( spam.version$version.string," is loaded.",
#       "\nType demo( spam) for some demos,",
#       " help( Spam) for an overview\nof this package.",
#       "\nHelp for individual functions is obtained by ",
#       "adding the\nsuffix '.spam' to the function name, e.g. 'help(chol.spam)'.")
       "\nType 'help( Spam)' or 'demo( spam)' for a short introduction ",
       "\nand overview of this package.",
       "\nHelp for individual functions is also obtained by ",
       "adding the\nsuffix '.spam' to the function name, e.g. 'help( chol.spam)'.")
   unlockBinding(".Spam", asNamespace("spam"))
 }

"spam.getOption" <- function (x)
  spam.options(x)[[1]]


"spam.options" <- function (...) {
    if (nargs() == 0) return(.Spam)
    current <- .Spam
    temp <- list(...)
    if (length(temp) == 1L && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.Spam[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
#    changed <- current[n]  #rf
    current[n] <- temp
    if (sys.parent() == 0) env <- asNamespace("spam") else env <- parent.frame()
    assign(".Spam", current, envir = env)
    invisible(current)
}

powerboost <- function(flag="on") {
    if (sys.parent() == 0)
            env <- asNamespace("spam")
    else env <- parent.frame()

    current <- spam.options()
    current[c("NAOK","safemodevalidity","cholsymmetrycheck","cholpivotcheck","eps")] <-
    if (tolower(flag) %in% c("true","on","an","ein")) {
            list(!FALSE,FALSE,FALSE,FALSE,1e-8)
    } else { list(!TRUE,TRUE,TRUE,TRUE,.Machine$double.eps)
    }
    assign(".Spam", current, envir = env)
    invisible( current)
}




timeDateOptions <-
function(...) {
  current <- if(exists(".splusTimeDateOptions", envir=.splusTimeDateEnv)) {
      .splusTimeDateEnv$.splusTimeDateOptions
  } else {
      .defaultSplusTimeDateOptions
  }
  if(!nargs()) {
    return(current)
  }
  a <- list(...)
  if(length(a) == 1 && is.null(names(a)) && is.list(a[[1]])) {
    a <- a[[1]]
    if(is.null(names(a))) {
            stop("list argument has no valid names")
          }
  }
  r <- vector("list", length(a))
  n <- names(a)
  show <- FALSE
  changed <- current
  for(i in seq_along(a)) {
    ni <- n[[i]]
    if(!is.null(ni) && nchar(ni)) {
      changed[[ni]] <- a[[i]]
    } else {
      ni = a[[i]]
      if(!is.character(ni)) stop("invalid argument")
      show <- TRUE
    }
    v <- current[[ni]]
    if(!is.null(v)) {
      r[[i]] <- if(is.null(a[[i]])) FALSE else v
    }
    names(r)[[i]] <- ni
  }
  assign(".splusTimeDateOptions", changed, envir=.splusTimeDateEnv)
  if(show) r else invisible(r)
}

.defaultSplusTimeDateOptions <- 
  list(ts.eps = 1e-05,
       sequence.tol = 1e-06,
       time.month.name = c("January", "February", "March", "April", "May",
         "June", "July", "August", "September", "October", "November",
         "December" ),
       time.month.abb = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
         "Aug", "Sep", "Oct", "Nov", "Dec"),
       time.day.name = c("Sunday", "Monday", "Tuesday", "Wednesday",
         "Thursday", "Friday", "Saturday" ),
       time.day.abb = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"),
       time.am.pm = c("AM", "PM"),
       time.century = 1930,
       time.zone = "GMT",
       time.out.format = "%02m/%02d/%04Y %02H:%02M:%02S.%03N",
       time.out.format.notime = "%02m/%02d/%04Y",
       time.in.format = "%m[/][.]%d[/][,]%y [%H[:%M[:%S[.%N]]][%p][[(]%3Z[)]]]",
       tspan.out.format = "%dd %Hh %Mm %Ss %NMS",
       tspan.in.format = "[%yy[ear[s]][,]] [%dd[ay[s]][,]] [%Hh[our[s]][,]] [%Mm[in[ute][s]][,]] [%Ss[ec[ond][s]][,]] [%NM[s][S]]"
       )

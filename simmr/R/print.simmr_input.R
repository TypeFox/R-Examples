print.simmr_input <-
function(x,...) {
  cat('This is a valid simmr input object with ')
  cat(paste(x$n_obs,'observations, '))
  cat(paste(x$n_tracers,'tracers, and '))
  cat(paste(x$n_sources,'sources.\n'))
  if(x$n_groups>1) cat(paste('There are',x$n_groups,'groups.\n'))
  cat('The source names are: ')
  cat(x$source_names,sep=', ')
  cat('.\n\n')
}

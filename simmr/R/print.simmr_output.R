print.simmr_output <-
function(x,...) {
  cat('This is a valid simmr output object with ')
  cat(paste(x$input$n_obs,'observations, '))
  cat(paste(x$input$n_tracers,'tracers, and '))
  cat(paste(x$input$n_sources,'sources.\n'))
  if(x$input$n_groups>1) cat(paste('There are',x$input$n_groups,'groups.\n'))
  cat('The source names are: ')
  cat(x$input$source_names,sep=', ')
  cat('.\n\n')
  cat('The input data has been run via simmr_mcmc and has produced ')
  cat(nrow(x$output[[1]][[1]]),'iterations over',length(x$output[[1]]),'MCMC chains.')
  cat('\n\n')
}

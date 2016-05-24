tv.home <- function(recheck = FALSE) {
  if(is.null(getOption('tv_home')) | recheck) {
    if(.Platform$OS.type == "unix") 
     if('Turbowin' %in% list.dirs(path=paste(Sys.getenv('HOME'),'/.wine/drive_c', sep=''), full.names=FALSE, recursive = FALSE))
          tv_home <- file.path(Sys.getenv('HOME'),'.wine/drive_c/Turbowin') else
          tv_home <- NA

  if(.Platform$OS.type == "windows") {
    if(file.access('O:/Turbowin/Popup/tvscale.dbf')==0) tv_home <- 'O:/Turbowin' else 
    	if(file.access('C:/Turbowin/Popup/tvscale.dbf')==0) tv_home <- 'C:/Turbowin' else 
    	  if(file.access('C:/Programs/Turbowin/Popup/tvscale.dbf') ==0) tv_home <- 'C:/Programs/Turbowin' else
    	    if(file.access('C:/Programme/Turbowin/Popup/tvscale.dbf') ==0) tv_home <- 'C:/Programme/Turbowin' else
    	      if(file.access('D:/Programme/Turbowin/Popup/tvscale.dbf') ==0) tv_home <- 'D:/Programme/Turbowin' else
              tv_home <- NA
  }
  if(is.na(tv_home)) {
   message('\nNo Turbowin installation path found. \n')
   if(interactive()) {
   ANSWER <- readline("Should I use \n 1) the vegdata package path (recommended),  or \n 2) a temporary folder? ")
   tv_home <- switch(substr(ANSWER, 1, 1),
         "1" = file.path(path.package('vegdata'), 'tvdata'),
         tempdir()
   ) } else tv_home <- tempdir()
   options(tv_home = tv_home)
   if(!file.exists(file.path(tv_home, 'tvdata', 'Popup', 'tvscale.dbf')))
     for(d in c('Popup', 'Data', 'Species')) {
     dir.create(file.path(tv_home, d), showWarnings = FALSE)
     if(d == 'Data') {
       wd <- getwd()
       setwd(file.path(path.package('vegdata'), 'tvdata', 'Data')  )
       dbs <- list.dirs(, recursive=TRUE, full.names=FALSE)
       for(l in 2:length(dbs)) {
         dir.create(file.path(tv_home, 'Data', dbs[l]), showWarnings = FALSE)
       file.copy(from =  list.files(dbs[l], recursive=TRUE, full.names=TRUE, include.dirs=TRUE), to = file.path(tv_home, 'Data', dbs[l]))
       }
       setwd(wd)
     } else
     file.copy(from =  list.files(file.path(path.package('vegdata'), 'tvdata', d), recursive=TRUE, full.names=TRUE, include.dirs=TRUE), to = file.path(tv_home, d))
   }
  }
  message('############################################################',
      '\nTurboveg root directory is set to "', tv_home, '"',
      '\nIf you want to change this use: options(tv_home=\"<path_to_your_Turbowin_root>\")',
      '\n############################################################')
  options(tv_home=tv_home)
  } else {
  tv_home <- getOption('tv_home')
#  message('Turboveg root directory has already been set to "', tv_home,'".\n', sep='')
  }
  invisible(tv_home)
}




# For calculating scores in spectrogram cross correlation
# Modified: 6 Sept 2015

corMatch <-
function(
  survey,                 # Complete survey which is to be analyzed for calls. Wave object or vector.
  templates,              # Template list, made with makeTemplates
  parallel=FALSE,         # If TRUE, mclapply is used for correlation calculations, for parallel processing (Linux or Mac OS X only). If FALSE lapply is used.
  show.prog=FALSE,        # If TRUE, progress is displayed during correlation calculations 
  cor.method='pearson',   # Method used by cor function (see ?cor)
  warn=TRUE,              # Set to FALSE to surpress warnings
  time.source='filename', # 'filename' or 'fileinfo' as the mtime source
  rec.tz=NA,              # Time zone setting for recorders 
  write.wav=FALSE,        # Set to TRUE to allow creation of file of survey in working directory
  ...                     # Additional arguments to the spectro function
) {

  # Check arguments
  if(missing(survey)) stop('Required argument survey is missing.')
  if(missing(templates)) stop('Required argument templates is missing.')

  # Packages
  if(parallel) {
    lapplyfun <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores=parallel::detectCores())
  } else lapplyfun <- lapply

  # Start tracking time (after loading packages)
  t.start <- Sys.time()

  # Work with survey outside template loop
  # Creates a wav file for survey if it isn't already a file
  survey <- getClip(survey, name=deparse(substitute(survey)), write.wav=write.wav)

  if(time.source == 'fileinfo') {
     file.time <- file.info(survey)$mtime
     if(is.na(rec.tz)) rec.tz <- format(file.time, format='%Z')
     file.time <- as.POSIXct(format(file.time, tz=rec.tz), tz=rec.tz)
  } else if(time.source == 'filename') {
     survey.short <- strsplit(survey, '/')[[1]][length(strsplit(survey, '/')[[1]])]
     date.time.info <- regmatches(survey.short, regexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}[ _][0-9]{6}[ _][A-Z0-9]{1,7}', survey.short))
     date.time.info <- gsub("_", " ", date.time.info)
     if(length(date.time.info) == 1 && nchar(date.time.info) %in% 19:23)
       file.time <- as.POSIXct(substr(date.time.info, start=1, stop=17), tz=substr(date.time.info, start=19, stop=length(date.time.info)), format='%Y-%m-%d %H%M%S') 
     else {
       warning('time.source was set to \"filename\" but file name does not have date and time info, so using \"fileinfo\" instead')
       file.time <- file.info(survey)$mtime
     }
  } else stop('time.source argument, ', time.source, ' not recognized.')

  survey.path <- survey
  survey <- readClip(survey)

  # score.L is a list for storing results
  score.L <- list()
  survey.data <- list()

  # Loop through templates
  for(i in names(templates@templates)) {
    cat('\nStarting ', i,'. . .')

    # Working with a single template
    template <- templates@templates[[i]]

    if(i == names(templates@templates)[1] || any(template@wl != wl, template@ovlp != ovlp, template@wn != wn)) {
      cat('\n\tFourier transform on survey . . .')
      wl <- template@wl
      ovlp <- template@ovlp
      wn <- template@wn
      # Perform Fourier transform on survey
      survey.spec <- spectro(wave=survey, wl=wl, ovlp=ovlp, wn=wn, ...)
      # NTS arbitrary adjustment to eliminate -Inf
      survey.spec$amp[is.infinite(survey.spec$amp)] <- min(survey.spec$amp[!is.infinite(survey.spec$amp)]) - 10
      frq.bins <- survey.spec$freq
      t.bins <- survey.spec$time
      t.survey <- length(survey@left)/survey@samp.rate
      t.step <- t.bins[2] - t.bins[1]
      frq.step <- frq.bins[2] - frq.bins[1]
      cat('\n\tContinuing. . .\n')
    }

    # Switch the order of columns in pt.on and pt.off to use them directly for indexing
    pts <- template@pts[, c(2:1, 3)]

    # Adjust pts if step sizes differ
    if(!isTRUE(all.equal(template@t.step, t.step, tolerance=t.step/1E4))) {
      pts[, 't'] <- round(pts[, 't']*template@t.step/t.step)
      if(warn) warning('For ', i,' time step doesn\'t match survey time step: ', t.step, ' != ', template@t.step)
    }
    if(!isTRUE(all.equal(template@frq.step, frq.step, tolerance=frq.step/1E6))) {
      pts[, 'frq'] <- round(pts[, 'frq']*template@frq.step/frq.step)
      if(warn) warning(i, ' frequency step does\'t match survey frequency step, ', frq.step, ' != ', template@frq.step)
    }

    # Determine the frequency limits from the template points
    frq.lim <- frq.bins[range(pts[, 'frq'])] 

    # Get number of time windows/bins in frequency domain data
    n.t.survey <- length(survey.spec$time)
  
    # Pare down amplitude matrix based on filter frequencies 
    which.frq.bins <- which(survey.spec$freq >= frq.lim[1] & survey.spec$freq <= frq.lim[2])
    amp.survey <- survey.spec$amp[which.frq.bins, ]

    # Shift frq indices in pts. The t indices already start at 1.
    pts[, 'frq'] <- pts[, 'frq'] - min(which.frq.bins) + 1
    n.t.template <- max(pts[, 't'])
    n.frq.template <- max(pts[, 'frq'])

    # Translate pts matrix of indices into a vector index so indexing is faster within the lapplyfun call
    pts.v <- (pts[, 't'] - 1)*n.frq.template + pts[, 'frq']
    amp.template <- pts[, 'amp']
    amp.survey.v <- c(amp.survey)  

    # Create progress bar object if requested
    if(show.prog & !parallel) pb <- txtProgressBar(max=(n.t.survey-n.t.template)*n.frq.template + 1, char='.', width=0, style=3)

    # Perform analysis for each time value (bin) of survey 
    # Starting time value (bin) of correlation window, set up as a list to use mclapply

    c.win.start <- as.list(1:(n.t.survey-n.t.template)*n.frq.template) # Starting position of subset of each survey amp matrix  
    score.survey <- unlist(
      lapplyfun(X=c.win.start, FUN=function(x) 
        {
        if(!parallel && show.prog) setTxtProgressBar(pb, x)
        # Unpack columns of survey amplitude matrix for correlation analysis
        cor(amp.template, amp.survey.v[x + pts.v], method=cor.method, use='complete.obs')  
        }
      )
    )

    # Collect score results and time (center of time bins) in data frame
    score.L[[i]] <- data.frame(
      date.time=file.time + survey.spec$time[1:(n.t.survey-n.t.template)+n.t.template/2] - t.survey, 
      time=survey.spec$time[1:(n.t.survey-n.t.template)+n.t.template/2], 
      score=score.survey
    )
    survey.data[[i]] <- list(amp=survey.spec$amp, t.bins=t.bins, frq.bins=frq.bins)
    cat('\n\tDone.\n')
  }

  # Calculate total run time
  t.run <- signif(as.numeric(difftime(Sys.time(), t.start, units='secs')), 4)
  time.info <- c(t.exe=as.character(t.run), RTfactor=paste(signif(t.survey/as.numeric(t.run, units='secs'), 4), 'x', sep=''))

  # Return results 
  scores.obj <- new('templateScores', survey.name=survey.path, survey=survey, survey.data=survey.data, templates=templates@templates, scores=score.L, time=time.info)
  return(scores.obj)
}

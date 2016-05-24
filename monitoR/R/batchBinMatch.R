# For running multiple surveys and templates from specified directories
# Modified: 2013 JUNE 14

batchBinMatch <-
function(
  dir.template,         # Directory with template files
  dir.survey='.',       # Directory with surveys 
  ext.template='bt',    # Extension for template files
  ext.survey='wav',     # Extension for survey files
  templates,            # Or provides template list directly
  parallel=FALSE,       # If TRUE, mclapply is used for correlation calculations, for parallel processing (Linux or Mac OS X only). If FALSE lapply is used. (corMatch & findPeaks)
  show.prog=FALSE,      # If TRUE, progress is displayed during correlation calculations  (binMatch)
  warn=TRUE,            # Set to FALSE to surpress warnings (binMatch)
  time.source='filename', # 'filename' or 'fileinfo' as the mtime source (binMatch)
  fd.rat=1,             # Factor to multiply template duration by for determining frame width (findPeaks)
  ...                   # Additional arguments to the spectro function
) {

  # Read in templates
  if(missing(templates)) {
    if(missing(dir.template)) stop('Both arguments dir.template and templates are missing--one is needed')
    cat('Reading in templates. .')
    templates <- readBinTemplates(dir=dir.template, ext=ext.template, parallel=parallel)
    cat(' . done.\n')
  } 

  # Create a vector of survey file names
  survey.files <- list.files(path=dir.survey, pattern=paste0('*\\.', ext.survey, '$')) 
  survey.files.full <- list.files(path=dir.survey, pattern=paste0('*\\.', ext.survey, '$'), full.names=TRUE) 

  detects <- NULL
  for(i in 1:length(survey.files)) {
    cat('\nCalling binMatch for', survey.files[i], '. . .\n')
    scores <- binMatch(survey=survey.files.full[i], templates=templates, parallel=parallel, show.prog=show.prog, warn=warn, time.source=time.source)
    cat('\nCalling findPeaks for output from', survey.files[i], '. . .\n')
    pks <- findPeaks(score.obj=scores, parallel=parallel)
    detects <- rbind(detects, getDetections(pks, id=survey.files[i]))
  }

  return(detects)
}

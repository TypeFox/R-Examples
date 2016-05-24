# Get the current working directory of running script from location of call
current_directory <-
  (function(frames = sys.frames()) {
    # http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
    frame_files <- Filter(Negate(is.null), lapply(frames, function(x) x$ofile))
    if (length(frame_files) == 0)
      stop("Ramd cannot identify what directory you are currently executing ",
           "in. Please specify it using base::setwd. (Note that certain Ramd",
           " features will not work from the console.")
    dirname(frame_files[[length(frame_files)]])
  })


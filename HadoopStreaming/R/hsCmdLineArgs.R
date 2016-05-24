hsCmdLineArgs <- function(spec=c(),openConnections=TRUE,args=commandArgs(TRUE)) {
  basespec = c(
    'mapper',     'm', 0, "logical","Runs the mapper.",FALSE,
    'reducer',    'r', 0, "logical","Runs the reducer, unless already running mapper.",FALSE,
    'mapcols',    'a', 0, "logical","Prints column headers for mapper output.",FALSE,
    'reducecols', 'b', 0, "logical","Prints column headers for reducer output.",FALSE,
    'infile'   ,  'i', 1, "character","Specifies an input file, otherwise use stdin.",NA,
    'outfile',    'o', 1, "character","Specifies an output file, otherwise use stdout.",NA,
    'skip',       's',1,"numeric","Number of lines of input to skip at the beginning.",0,
    'chunksize',  'c',1,"numeric","Number of lines to read at once, a la scan.",-1,
    'memlim',     'z',1,"numeric","Max number of bytes allowed for use in carry.",-1,
    'numlines',   'n',1,"numeric","Max number of lines to read from input, per mapper or reducer job.",0,
    'sepr',       'e',1,"character","Separator character, as used by scan.",'\t',
    'insep',      'f',1,"character","Separator character for input, defaults to sepr.",NA,
    'outsep',     'g',1,"character","Separator character output, defaults to sepr.",NA,
    'debug',      'd',0,"logical","Turn on debugging output.",FALSE,
    'help',       'h',0,"logical","Get a help message.",FALSE
    )

  specmat = matrix(c(spec,basespec),ncol=6,byrow=TRUE)

  opt = getopt(specmat[,1:5],opt=args)
  ## Set Default parameter values
  for (p in seq(along=specmat[,1])) {
    s = specmat[p,1]
    if (is.null(opt[[specmat[p,1]]]) ) {
      opt[[specmat[p,1]]] = specmat[p,6]
      storage.mode( opt[[specmat[p,1]]] ) <- specmat[p,4]
    }
  }

  ## Set separator character
  if (is.na(opt$insep)) opt$insep=opt$sep
  if (is.na(opt$outsep)) opt$outsep=opt$sep

  ## Print help, if necessary
  if ( opt$help || (!opt$mapper && !opt$reducer && !opt$mapcols && !opt$reducecols)  ) {
    ##Get the script name (only works when invoked with Rscript).
    ## self = commandArgs()[1];
    cat(getopt(specmat,usage=TRUE))
    opt$set = FALSE
##    return(opt)
  }

  if (openConnections) {
    if (is.na(opt$infile) && (opt$numlines==0)) {
      opt$incon = file(description="stdin",open="r")
    } else if  (is.na(opt$infile) && (opt$numlines>0)) {
      opt$incon = pipe( paste("head -n",opt$numlines), "r" )
    } else if  (opt$numlines>0) {
      opt$incon = pipe( paste("head -n",opt$numlines,opt$infile), "r" )
    } else if (opt$numlines==0) {
      opt$incon = file(description=opt$infile,open="r")
    } else {
      stop("I can't figure out what's going on with this input stuff.")
    }

    if (is.na(opt$outfile) ) {
      opt$outcon = ""
    } else {
      opt$outcon = file(description=opt$outfile,open="w")
    }
  }

  opt$set = TRUE
  return(opt)
}

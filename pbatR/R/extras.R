catwd <- function() {
  cat( paste("Working directory changed to '", getwd(), "'.", sep="", "\n") )
  return(invisible())
}

pbat.work <- function(pedOrPhe=NULL) {
  prevwd <- getwd()

  if( is.null(pedOrPhe) || !is.sym(pedOrPhe) ) {
    if( !is.null(pedOrPhe) ) warning( "pedOrPhe not symbolic, ignored." )
    wd <- getwd()
  }else{
    wd <- str.getpath( get.sym(pedOrPhe) )
  }
  setwd( wd )
  if( !file.exists("./pbatRwork") )
    dir.create( "./pbatRwork" )
  wd <- paste( wd, "/pbatRwork", sep="" )
  cat( wd, "\n" )

  setwd(wd)
  catwd()

  return( prevwd )
}

pbat.unwork <- function(cur=NULL) {
  if( is.null(cur) ) {
    setwd("..")
  }else{
    setwd(cur)
  }
  catwd()
}

pbat.status <- function(n=1,workFirst=FALSE) {
  ofile <- "pbatstatus.txt"
  if( workFirst ) ofile <- "pbatRwork/pbatstatus.txt"

  if( !file.exists(ofile) ) {
    ofile <- "pbatRwork/pbatstatus.txt"
    if( workFirst ) ofile <- "pbatstatus.txt"

    if( !file.exists(ofile) )
      return("'pbatstatus.txt' does not exist, no status information available.")
  }

  f <- file( ofile )
  lines <- readLines(f)
  close(f)

  if( n==0 ) return(lines)  ## return everything

  return( tail(lines,n=n) )
}

#extrasDebug <- function() {
#  library( pbatR )
#  is.sym <- getFromNamespace( "is.sym", "pbatR" )
#  get.sym <- getFromNamespace( "get.sym", "pbatR" )
#  str.getpath <- getFromNamespace( "str.getpath", "pbatR" )
#
#  catwd()
#  pbat.work()
#  pbat.unwork()
#
#  ped <- read.ped( "~/pbatGUI/data/junk" )
#  cur <- pbat.work(ped)
#  pbat.unwork(cur)
#}

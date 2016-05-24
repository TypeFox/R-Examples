##' Create an Emu segment list from a file
##' 
##' Create an Emu segment list from a file saved by the Emu query tools.
##' 
##' Reads segment lists created by programs external to R/Splus and stored in
##' text files on disk.
##' 
##' @param file The name of the file to read
##' @return An Emu segment list.
##' @author Steve Cassidy
##' @seealso \code{\link{query}}
##' @keywords IO
##' @examples
##' 
##' ## create a segment list file and write it out
##' # seglist.txt <- "database:demo"\
##' # query:Phonetic=vowel\
##' # type:segment\
##' #\
##' # @@:	3059.65	3343.65	msdjc001\
##' # e:	5958.55	6244.55	msdjc002\
##' # @@u	8984.75	9288.75	msdjc003\
##' # A	11880.8	12184.8	msdjc004\
##' # E	17188.3	17366.4	msdjc005\
##' # ei	20315.2	20655.2	msdjc006"
##' 
##' \dontrun{cat(seglist.txt, file="seglist.txt")}
##' 
##' # now read it back as a segment list
##' \dontrun{segs <- read.emusegs("seglist.txt")}
##' \dontrun{segs}
##' ## and clean up
##' \dontrun{unlink("seglist.txt")}
##' 
##' 
##' @export read.emusegs
read.emusegs <- function(file)
{
  ## scan the lines of the file into a vector
  
  ## R 1.4 introduced comment.char="#" arg to scan, grrr
  if( is.R() && as.numeric(version$minor) > 3.0 ) {
    ## in R, we need to avoid skipping the # as a comment line
    lines <- scan(file, what = "", sep="\n", comment.char="")
  } else {
    lines <- scan(file, what = "", sep="\n")
  }
  
  ## first three lines are header followed by a hash line
  inheader <- 1
  i <- 1
  labels <- start <- end <- utts <- NULL
  while( i < length(lines) && inheader) {
    if( lines[i] == "#" ) {
      inheader <- 0
    } else {
      foo <- splitstring( lines[i], ":" )
      if( foo[1] == "database" ) database <- paste(foo[-1], sep=":")
      if( foo[1] == "query" ) {
        query <- paste(foo[-1], sep=":")
      }
      if( foo[1] == "type" ) type <- paste(foo[-1], sep=":")
      i <- i + 1
    }
  }
  
  if (inheader) {
    stop( "End of header (#) not found in segment file" )
  }
  
  ## now slurp the body of the segment list
  mat <- matscan( file, 4, what="", sk=i )
  
  segs <- make.seglist(mat[,1], mat[,2], mat[,3], mat[,4], 
                       query, type, database )
  
  segs
}

if( version$major >= 5 ) {
  setOldClass(c("emusegs", "data.frame"))
}









##' Make an Emu segment list from the various components
##' 
##' This is the appropriate way to make an Emu segment list and ensure that it
##' has all of the required components.
##' 
##' An Emu segment list is the result of a query to a speech database (see
##' \code{\link{query}}) and has one row per matching segment or event from
##' the query. Each row lists the label, start and end times (in milliseconds)
##' and utterance name for the segment.  This information is used by
##' \code{\link{get_trackdata}} and other functions to extract data corresponding
##' to these segments.
##' 
##' In order to ensure the proper format for segment lists and to ensure
##' against future changes to the format, \code{make.seglist} should be used
##' whenever you wish to create a segment list.  Another function,
##' \code{\link{modify.seglist}} can be used to change some part of an existing
##' segment list. The functions \code{\link{label.emusegs}},
##' \code{\link{start.emusegs}}, \code{\link{end.emusegs}} and
##' \code{\link{utt.emusegs}} can be used to access the different columns of
##' the segment list.
##' 
##' @param labels A character vector of labels for each segment
##' @param start A vector of start times
##' @param end A vector of end times
##' @param utts A character vector of utterance names
##' @param query A query string
##' @param type \code{segment} or \code{event}
##' @param database The database name associated with the segment list
##' @return An Emu segment list.
##' @author Steve Cassidy
##' @seealso \code{\link{modify.seglist}}, \code{\link{label.emusegs}}
##' @keywords misc
##' @examples
##' 
##' 
##'    l <- c("A", "B", "C")
##'    s <- 1:3
##'    e <- 2:4
##'    u <- c("u1", "u1", "u1")
##'    segs <- make.seglist(l, s, e, u, "Fake Query", "segment", "fake")
##'    segs
##'    ## summary gives an overview of the data in the segment list
##'    summary(segs)
##'    
##' 
##'    # The following should be TRUE
##'    label(segs) == l
##'    dur(segs) == s
##'    end(segs) == e
##'    utt(segs) == u
##'    emusegs.database(segs) == "fake"
##'    emusegs.type(segs) == "segment"
##'    emusegs.query(segs) == "Fake Query"
##' 
##'    # segment durations should all be 1
##'    dur(segs) == c(1,1,1)
##' 
##' 
##' @export make.seglist
make.seglist <- function(labels, start, end, utts, query, type, database)
{
  seglist <- data.frame(labels=I(as.character(labels)),
                        start=as.numeric(start), 
                        end=as.numeric(end), 
                        utts=I(as.character(utts)))
  
  if( version$major >= 5 ) {
    oldClass(seglist) <- "emusegs"
  } else {
    class(seglist) <- c("emusegs", "data.frame")
  }
  
  attr(seglist, "query") <- query
  attr(seglist, "type") <- type
  attr(seglist, "database") <- database
  
  seglist
}









##' is seglist
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export is.seglist
is.seglist <- function(object) {
  return( inherits(object, "emusegs") )
}

## modify a segment list by changing one or more of the fields








##' Modify one of the components of an Emu segment list
##' 
##' This function can be used to modify one of the parts of an Emu segment list
##' while leaving the other parts unchanged.
##' 
##' An Emu segment list has a number of components and is stored as an R object
##' of class \code{emusegs}.  This function can be used to modify a segment
##' list while retaining all of the proper structures.
##' 
##' Any new vectors passed to the function must have the same length as the
##' segment list itself for this call to succeed.
##' 
##' All arguments are optional and default to not modifying the segment list if
##' not supplied.
##' 
##' The original segment list is not modified, instead, a modified copy is
##' returned.
##' 
##' @param segs A segment list to modify, a modified copy is returned
##' @param labels A new label vector
##' @param start A new start time vector
##' @param end A new end time vector
##' @param utts A new vector of utterance labels
##' @param query A new query string to associate with the segment list
##' @param type A new type string
##' @param database A new database name
##' @return An Emu segment list.
##' @author Steve Cassidy
##' @seealso \code{\link{query}}
##' @keywords misc
##' @examples
##' 
##' data(vowlax)
##' segs = vowlax
##' # extend the start times by 10ms
##' newsegs <- modify.seglist( segs, start=start(segs)+10 )
##' 
##' # change the associated database name
##' # this will affect where emu.track looks to find data
##' newsegs <-  modify.seglist( segs, database="notdemo" )
##' 
##' 
##' @export modify.seglist
"modify.seglist" <- function( segs,
                              labels=label.emusegs(segs),
                              start=start.emusegs(segs),
                              end=end.emusegs(segs),
                              utts=utt.emusegs(segs),
                              query=emusegs.query(segs),
                              type=emusegs.type(segs),
                              database=emusegs.database(segs))
{
  make.seglist( labels, start, end, utts,
                query, type, database )
}









##' emusegs database
##' 
##' Returns the database attribute from a segmentlist
##' 
##' 
##' @keywords internal
##' @export emusegs.database
"emusegs.database" <- function(sl) 
{ 
  if(is.seglist(sl))
    attr(sl, "database")
  else 
    stop( "not an emu segment list" )
}









##' segment list type
##' 
##' Gives SEGEMENT or EVENT
##' 
##' 
##' @keywords internal
##' @export emusegs.type
"emusegs.type" <- function(sl) 
{ 
  if(is.seglist(sl))
    attr(sl, "type")
  else 
    stop( "not an emu segment list" )
}









##' emusegs query
##' 
##' sends a emu query to EMU
##' 
##' 
##' @keywords internal
##' @export emusegs.query
"emusegs.query" <- function(sl) 
{ 
  if(is.seglist(sl))
    attr(sl, "query")
  else 
    stop( "not an emu segment list" )
}









##' print emusegs
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export
"print.emusegs" <-  function(x, ...)
{
  cat(attributes(x)$type, " list from database: ", attributes(x)$database, "\n")
  cat("query was: ", attributes(x)$query, "\n" )
  if( version$major >= 5 ) {
    oldClass(x) <- "data.frame"
  } else {
    class(x) <- "data.frame"
  }
  print.data.frame(x, ...)
}


##' summary emusegs
##' 
##' summarizes data in emu segment lists
##' 
##' 
##' @param object the segmentlist
##' @param \dots nothing special
##' @keywords internal
##' @method summary emusegs
##' @export
summary.emusegs <- function(object, ...)
{
  cat(attributes(object)$type, " list from database: ", attributes(object)$database, "\n")
  cat("query was: ", attributes(object)$query, "\n" )
  cat(" with", length(object$start), "segments\n\n")
  cat("Segment distribution:\n")
  print(table(object$label))
  invisible()
}









##' Get labels / utterances from segment list
##' 
##' label: extracts the labels from the segment list.  utt: extracts the
##' utterances from the segment list.
##' 
##' 
##' @aliases label.emusegs label utt.emusegs utt
##' @param segs segment list
##' @return label / utterance vector
##' @author Jonathan Harrington
##' @seealso \code{\link{segmentlist} \link{start} \link{end}}
##' @keywords methods
##' @examples
##' 
##'    data(dip)
##'    #dip is a segment list - first ten segments only
##'    dip[1:10,]
##'    
##' 
##'    #extract labels from the segment list
##'    dips.labs = label(dip)
##'    dips.labs 
##'    
##' 
##' @export label
"label" <- function(segs) {
  UseMethod("label")
}




##' @export
"label.emusegs" <- function(segs)
{
  as.character(segs$label)
}









##' as.matrix.emusegs
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export
"as.matrix.emusegs" <- function(x, ...)
{
  cbind( as.character(x$label), x$start, x$end, as.character(x$utt) )
}









##' Write an Emu segment list to a file
##' 
##' Writes an Emu segment list to a file
##' 
##' 
##' @param seglist An Emu segment list
##' @param file The name of a file to write the segment list into.
##' @return None.
##' @section Side Effects: The segment list is written to a file in the
##' standard format, suitable for input to \code{gettrack} or other Emu utility
##' programs.
##' @seealso \code{\link{query}}
##' @keywords misc
##' @examples
##' 
##'    data(dip)
##'    #dip a segment list - first 10 segments only
##'    dip[1:10,]
##'    \dontrun{write.emusegs(dip, "write.emusegs.example.txt")}
##'    
##'    #The file write.emusegs.example.txt would have been written to R_HOME
##'    \dontrun{unlink("write.emusegs.example.txt")}
##' 
##' @export write.emusegs
"write.emusegs" <- function(seglist, file)
{
  if(inherits(seglist,"emuRsegs")){
    warning("You are using the write function of the legacy class emusegs for an emuRsegs object. The persisted object cannot be read back as emuRsegs object. It is recommended to use standard R function save() instead to persist an emuRsegs object.")
  }
  cat(paste("database:", attributes(seglist)$database, "\n", sep=""), file=file)
  cat(paste("query:", attributes(seglist)$query, "\n", sep=""), file=file, append=TRUE)
  cat(paste("type:", attributes(seglist)$type, "\n", sep=""), file=file, append=TRUE)
  cat("#\n", file=file, append=TRUE)
  write(t(as.matrix(seglist)), file, ncolumns = 4, append=TRUE)
}


##' @export
"start.emusegs" <- function(x, ...)
{
  as.numeric(x$start)
}


##' @export
"end.emusegs" <- function(x, ...)
{
  as.numeric(x$end)
}


##' @export
"utt" <- function(x) {
  UseMethod("utt")
}

##' @export
"utt.emusegs" <- function(x)
{
  as.character(x$utts)
}


##' duration
##' 
##' calculates durations
##' 
##' @param x ???
##' @export
"dur" <- function(x) {
  UseMethod("dur")
}









##' Duration of segments
##' 
##' duration of segments is calculated for each segment in the segment list
##' 
##' 
##' @param x a segment list
##' @return a vector of durations
##' @author Jonathan Harrington
##' @keywords internal
##' @export
"dur.emusegs" <- function (x) 
{
  if(all(end(x)==0))
    d <- end(x)
  else
    d <-  end(x) - start(x)
  d
}

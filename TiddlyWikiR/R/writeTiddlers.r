##writeTiddlers.r
##2013-06-13 dmontaner@cipf.es
##TiddlyWikiR library

## Some auxiliary functions:

################################################################################

## @name read.tiddler.titles
## @author David Montaner \email{dmontaner@@cipf.es}
## 
## @keywords read tiddler names
## @seealso \code{\link{writeTiddlers}}
## 
## @title Reads the titles from the tiddlers in an existing TiddlyWiki file. 
## 
## @description Utility function not to overwrite tiddlers or create several of them with the same name.
## 
## @param file an input TiddlyWiki file.
## @param verbose verbose.

## @return a character vector with the titles.

read.tiddler.titles <- function (file, verbose = TRUE) {
  
  lineas <- readLines (file)
  
  ##Read and split the template
  if (verbose) cat ("Reading file:", file, fill = TRUE)
  templ <- readLines (con = file)
  ##
  pivot.line <- grep ('<div id="storeArea">', templ)  ##unique tag
  total.lines <- length (templ)
  
  s.line <- grep ('<!--POST-SHADOWAREA-->', templ)
  e.line <- grep ('<!--POST-STOREAREA-->', templ)
  
  t.lines <- templ[s.line:e.line]

  t.lines <- grep ("<div title=", templ[s.line:e.line], value = TRUE)
  t.lines <- gsub ('<div title=\"', '', t.lines)
  
  t.lines <- strsplit (t.lines, split = '\"')
  t.lines <- sapply (t.lines, function (x) x[1])

  t.lines <- as.character (t.lines) ##just in case it is an empty list
  
  return (t.lines)
}

################################################################################

paste.tags <- function (x) {
  ##To Do: a trimming would be nice before the grep
  conespacio <- grep (" ", x)
  x[conespacio] <- paste ("[[", x[conespacio], "]]", sep = "")
  x <- paste (x, collapse = " ")
  return (x)
}

########################################

wikify.tiddler.content <- function (tid) {
  out <- unlist (lapply (tid@content, wikify))
}


wikify.tiddler.left <- function (tid) {
  out <- paste ('<div title="', tid@title,
                '" creator="',  tid@creator,
                '" modifier="', tid@modifier,
                '" created="',  tid@created,
                ifelse (length (tid@modified) == 0, "", paste ('" modified="', tid@modified, sep = "")),
                ifelse (length (tid@tags)     == 0, "", paste ('" tags="', paste.tags (tid@tags), sep = "")),
                '" changecount="', tid@changecount,
                '">', sep = "")
  return (out)
}

wikify.tiddler <- function (tid) {
  out <- c (wikify.tiddler.left (tid),
            "<pre>",
            wikify.tiddler.content (tid),
            "</pre>",   ##this does the
            "</div>")   ##wikify.tiddler.right
  return (out)
}

########################################

wikify.tiddler.list <- function (tiddlerList) {
  out <- lapply (tiddlerList, wikify.tiddler)
  out <- unlist (out)
  return (out)
}

########################################

setMethod ("wikify", "tiddler", function (object) { ## may be it should not be here
  out <- wikify.tiddler (object)
  return (out)
})

################################################################################

##' @name writeTiddlers
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases insert
##' @keywords insert tiddler
##' @seealso \code{\link{wikify}}
##' 
##' @title Inserts one or more tiddlers into a TiddlyWiki file.
##' 
##' @description Inserts a tiddler into the desired TiddlyWiki file.
##'
##' @details By default the TiddlyWiki template "file" will be overwritten.
##' The options "infile" and "outfile" may be used not to overwrite the template.
##' 
##' @param tid a tiddler or a list of tiddlers to be inserted.
##' @param file TiddlyWiki template file.
##' @param infile the template file where the tiddler is to be inserted.
##' @param outfile the output file. Intended not to overwrite the template.
##' @param verbose verbose.
##' @param check check that there is not already a tiddler with the same title.
##' @export

writeTiddlers <- function (tid = NULL,
                           file, infile = file, outfile = file,
                           check = TRUE, verbose = TRUE) {
  
  if (is.null (tid)) { ## NULL tiddler is intended just for copying the template file.
    if (infile != outfile) {
      cp <- file.copy (from = infile, to = outfile, overwrite = TRUE)  ##is it dangerous the overwrite???
      if (verbose & cp) {
        cat ("Writing file:", outfile, fill = TRUE)
      }
    }
  } else { ##check the object format
    if (is.list (tid)) {
      if (any (sapply (tid, class) != "tiddler")) {
        stop ("all elements in the list should be of class 'tiddler'.")
      }
    } else {
      if (class (tid) == "tiddler") {
        tid <- list (tid)
      } else {
        stop ("'tid' must be of class 'tiddler' or a list of them.")
      }
    }
    
    if (check) {## check if tiddler names are already in the TiddlyWiki template
      titulos <- sapply (tid, tdTitle)
      titulos.comunes <- intersect (titulos, read.tiddler.titles (infile, verbose = FALSE))
      ##
      if (length (titulos.comunes) > 0) {
        stop ("Tiddler with titles: \n",
              paste (paste ("   ", titulos.comunes), concatenate = "\n"),
              "already exist in file ", infile)
      }
    }
    
    ## #########################################################################
    
    ##read, split and insert the template
    if (verbose) {
      cat ("Reading file:", infile, fill = TRUE)
    }
    
    templ <- readLines (con = infile)
    pivot.line <- min (grep ('<div id="storeArea">', templ))  ##is a unique tag unless is also written as content in the tiddlers
    total.lines <- length (templ)
    ##
    templ <- c (templ[1:(pivot.line-1)],
                '<div id="storeArea">',
                wikify.tiddler.list (tid),
                sub ('<div id="storeArea">', '', templ[pivot.line]),  ##just in case there is something additional in the storeArea line
                templ[(pivot.line + 1):total.lines])
    
    ##SAVE
    if (verbose) {
      cat ("Writing file:", outfile, fill = TRUE)
    }
    writeLines (text = templ, con = outfile)
  }
}
################################################################################

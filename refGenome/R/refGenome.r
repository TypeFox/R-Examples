
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   refGenome                                                  ##
##  Created   :   19.03.2012                                                 ##
##  Author    :   W. Kaisers                                                 ##
##  Content   :   Funktionality for importing and managing genomic reference ##
##                                                                  data     ##
##                (Ucsc, Ensembl, Genbank)                                   ##
##                for usage in R                                             ##
##  Version   :   1.1.1                                                      ##
##                                                                           ##
##  Changelog :                                                              ##
##  05.06.12  :   addIsoforms and addEnsemble delete qualifier rows before   ##
##                                                      (re-) inserting      ##
##  06.06.12  :   Added correction in get_ens_attribute_df which removed a   ##
##                                                      memory leak          ##
##                (Also the token_list module is now valgrind tested)        ##
##  08.06.12  :   Added class refdb         (encapsulates database access)   ##
##  09.06.12  :   Added class refDataLocations                               ##
##                (encapuslates directory and file management)               ##
##  13.06.12  :   Implementation updates for refGenome,ensembl,ucsc and      ##
##                                                  genbank finished         ##
##  26.11.12  :   Added function 'extractSeqids' and 'tableSeqids'           ##
##  08.05.13  :   refGenome_1.0.0 on CRAN                                    ##
##  09.05.13  :   Added class and function ensemblExons                      ##
##  06.06.13  :   Added strand and frame in 'getGenePositions'               ##
##  17.07.13  :   Changed signature for 'extractByGeneName' 1.0.4            ##
##                (so generic can be used in 'spliceSites')                  ##
##  01.08.13  :   C-routines valgrind tested                                 ##
##  04.08.13  :   getGenePositions changed (doBy): >116 sec to 3.9 sec       ##
##                                                  runtime (1.0.6)          ##
##  04.08.13  :   Added getSpliceTable, unifyJuncs (1.0.7)                   ##
##  05.08.13  :   Added getUnifiedJuncs, updated vignette (1.0.8)            ##
##  06.08.13  :   New getSplicSite and unifyJuncs in C (1.0.10),             ##
##                                              valgrind tested              ##
##  07.08.13  :   refGenome_1.1.0 on CRAN                                    ##
##  08.08.13  :   Corrected generic for extractByGeneName (refGenome_1.1.0)  ##
##  04.12.13  :   refGenome_1.2.0 on CRAN                                    ##
##                                                                           ##
##  18.09.14  :   refGenome_1.3.0 on CRAN                                    ##
##                                                                           ##
##  07.07.14  :   Added R_init_refGenome                                     ##
##  08.07.14  :   Added overlapJuncs function (1.2.5)                        ##
##  10.07.14  :   Added tests which are executed in R CMD check              ##
##                                                                           ##
##  30.06.15  :   New C++ based implementation of GTF import (1.4.2)         ##
##  20.07.15  :   refGenome_1.5.6 on CRAN                                    ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

.onUnload<-function(libpath) { library.dynam.unload("refGenome", libpath) }

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
## refGenome      Functionality for processing Sequence and Annotation data  ##
##                in existing File (and dir) structures.                     ##
##                                                                           ##
##                refdir/(input data)                                        ##
##                refdir/seqs                                                ##
##                refdir/dbname.db3                                          ##
##                                                                           ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
## Split gtf Attribute column data                                           ##
## and return list with two data.frames                                      ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## A: The Brent Lab (Washington University, St.Louis)                        ##
## http:#mblab.wustl.edu/GTF2.html                                           ##
## [attributes] All four features have the same two mandatory attributes at  ##
## the end of the record:                                                    ##
##                                                                           ##
## gene_id value;       A globally unique identifier for the genomic source  ##
##                      of the transcript                                    ##
## transcript_id value; A globally unique identifier for the predicted       ##
##                      transcript.                                          ##
##                                                                           ##
## These attributes are designed for handling multiple transcripts from the  ##
## same genomic region.                                                      ##
## Any other attributes or comments must appear after these two and will be  ##
## ignored.                                                                  ##
##                                                                           ##
## Attributes must end in a semicolon which must then be separated from the  ##
## start of any subsequent                                                   ##
## attribute by exactly one space character (NOT a tab character). Textual   ##
## attributes !should! be surrounded by doublequotes.                        ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## B: Wellcome Trust Sanger Institute                                        ##
## http:#www.sanger.ac.uk/resources/software/gff/spec.html                   ##
##                                                                           ##
## Free text values *must* be quoted with double quotes.                     ##
## Note: all non-printing characters in such free text value strings (e.g.   ##
## newlines, tabs, control characters, etc) must be explicitly represented   ##
## by their C (UNIX) style backslash-escaped representation                  ##
## (e.g. newlines as '\n', tabs as '\t').                                    ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Static functions                                                          ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
ucPrimAssembly <- function() {return("^chr[0-9XYM]{1,2}$")}
ensPrimAssembly <- function(){return("^([0-9]{1,2})$|^[XY]|MT$")}
strandlevels <- c("+", "-", "*")


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
## Class declarations                                                        ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Genome classes
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

# con can also be NULL
setClass("refGenome",
         representation(
             ev="environment",
             basedir="character",
             "VIRTUAL"),
         prototype=prototype(
             ev=new.env(),
             basedir="."
         )
)

setClass("ucscGenome", contains = "refGenome")
setClass("ensemblGenome", contains = "refGenome")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Exon clases
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setClass("refExons", representation("VIRTUAL"), contains = "refGenome")
setClass("ensemblExons", contains = "refExons")
setClass("ucscExons", contains = "refExons")
setClass("unifiedExons", contains="refExons")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Junction classes
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

# Environment contains 'gtf' table which contains junction data.
setClass("refJunctions", representation("VIRTUAL"), contains = "refGenome")
setClass("ensemblJunctions", contains = "refJunctions")
setClass("ucscJunctions", contains = "refJunctions")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
## Declarations of Generics                                                  ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setGeneric("basedir", function(object) standardGeneric("basedir"))
setGeneric("basedir<-", function(object, value) standardGeneric("basedir<-"))

setGeneric("getGtf", function(object) standardGeneric("getGtf"))
setGeneric("setGtf", function(object, value) standardGeneric("setGtf"))

setGeneric("getGeneTable", function(object) standardGeneric("getGeneTable"))

setGeneric("read.gtf", function(object, filename="transcripts.gtf",
                                    sep = "\t", useBasedir = TRUE,
                                    comment.char = "#",
                                    progress=100000L,
                                    ...) standardGeneric("read.gtf"))


setGeneric("saveGenome", function(object, filename, useBasedir = TRUE, ...)
                        standardGeneric("saveGenome"))

setGeneric("writeDB", function(object, filename, useBasedir = TRUE, ...)
                        standardGeneric("writeDB"))

setGeneric("addIsoforms", function(object, filename = "ucsc_knownisoforms.csv",
            sep = "\t", useBasedir = TRUE, ...) standardGeneric("addIsoforms"))

setGeneric("addEnsembl", function(object, filename = "ucsc_knownToEnsembl.csv",
            sep = "\t", useBasedir = TRUE, ...) standardGeneric("addEnsembl"))

setGeneric("addXref", function(object, filename = "kgXref.csv",
                sep = "\t", useBasedir = TRUE, ...) standardGeneric("addXref"))

setGeneric("getXref", function(object) standardGeneric("getXref"))

setGeneric("tableFeatures", function(object) standardGeneric("tableFeatures"))

setGeneric("extractFeature", function(object, feature)
                                            standardGeneric("extractFeature"))

setGeneric("extractTranscript", function(object, transcripts)
                                        standardGeneric("extractTranscript"))

setGeneric("extractByGeneName", function(object, geneNames, src, ...)
                                        standardGeneric("extractByGeneName"))

setGeneric("extractByGeneId", function(object, geneids, ...)
                                        standardGeneric("extractByGeneId"))

setGeneric("getGenePositions", function(object, by, force = FALSE, ...)
                                        standardGeneric("getGenePositions"))

setGeneric("tableTranscript.id",
                        function(object) standardGeneric("tableTranscript.id"))

setGeneric("tableTranscript.name", function(object)
                                        standardGeneric("tableTranscript.name"))

setGeneric("tableSeqids",
                        function(object, regex) standardGeneric("tableSeqids"))

setGeneric("extractSeqids",
                    function(object, regex) standardGeneric("extractSeqids"))

setGeneric("extractPaGenes", function(object) standardGeneric("extractPaGenes"))

setGeneric("addExonNumber", function(object) standardGeneric("addExonNumber"))

setGeneric("refExons", function(object) standardGeneric("refExons"))

setGeneric("getSpliceTable", function(object) standardGeneric("getSpliceTable"))

setGeneric("unifyJuncs", function(object) standardGeneric("unifyJuncs"))

setGeneric("addIsCoding", function(object, ens) standardGeneric("addIsCoding"))




## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Section: refGenome - class                                               ##
##  Virtual base class for ucscGenome and ensemblGenome                      ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("initialize", "refGenome", function(.Object)
{
  # Has to be done explicitly here
  # because otherwise obscure copies appear
  .Object@ev <- new.env()
  assign("rgClass", class(.Object), envir = .Object@ev)
  return(.Object)
})

loadGenome <- function(src, basedir)
{
    if(is.character(src))
    {
        if(missing(basedir))
        basedir <- dirname(src)
    else
    {
        if(dirname(src) == ".")
            src <- file.path(basedir, src)
    }
    # Has to be done this way
    # because class needs to be read from ev
    ev <- new.env()
    load(src, envir=ev)
    rg <- new(ev$rgClass)
    basedir(rg) <- basedir
    rg@ev <- ev

    return(invisible(rg))
  }
  if(is(src, "url"))
  {
    ev <- new.env()
    load(src, envir = ev)
    rg <- new(ev$rgClass)
    basedir(rg) <- basedir
    rg@ev <- ev
  }
  stop("Invalid 'src' argument!")
}

setMethod("show", "refGenome", function(object)
{
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    cat("An empty object of class '", class(object), "'.\n", sep="")
  else
  {
    n<-min(nrow(object@ev$gtf), 6L)
    bm<-Sys.localeconv()[7]
    cat("Object of class '", class(object), "' with ",
                format(nrow(object@ev$gtf), big.mark = bm), " rows and ",
                ncol(object@ev$gtf), " columns.\n", sep="")

    print(object@ev$gtf[1:n, ])
  }
})

setMethod("basedir", "refGenome", function(object) {return(object@basedir)})

setReplaceMethod("basedir", "refGenome", function(object, value)
{
  if(!file.exists(value))
    cat("[basedir.refGenome] Directory '", value, "' does not exist!\n", sep="")

  object@basedir<-value
  return(object)
})


setMethod("getGtf", "refGenome", function(object)
{
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    return(NULL)

  return(object@ev$gtf)
})


setMethod("setGtf", "refGenome", function(object, value)
{
  assign("gtf", value, envir=object@ev)
  return(invisible())
})

setMethod("read.gtf", "refGenome",
          function(object, filename="transcripts.gtf",
                   sep = "\t", useBasedir=TRUE,
                   comment.char = "#",
                   progress=100000L,
                   ...)
{
    if(useBasedir && (length(object@basedir)>0))
        filename <- file.path(object@basedir, filename)

    if(!file.exists(filename))
        stop("File '", filename, "' does not exist!\n", sep="")

    bm <- Sys.localeconv()[7]

    cat("[read.gtf.refGenome] Reading file '",
        basename(filename), "'.\n", sep="")

    param <- c(path.expand(filename), comment.char[1], sep)
    tbl <- .Call("read_gtf", param, progress, PACKAGE="refGenome")

    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Ensembl 76 has introduced 'gene' feature type:
    ## Second item in attributes is "gene_name"
    ## (unlike 'transcript_id').
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    cat("[read.gtf.refGenome] Extracting genes table.\n")
    genl <- tbl$feature=="gene"
    ngen <- sum(genl)

    if(ngen > 0)
    {
        cat("[read.gtf.refGenome] Found",
            format(ngen, big.mark = bm), "gene records.\n")
        genes <- tbl[genl, ]
        genes$feature <- NULL
        tbl <- tbl[!genl, ]
    }


    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Copy result tables into object environment
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    assign("gtf", tbl, envir = object@ev)
    if(ngen > 0)
        assign("genes", genes, envir = object@ev)


    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Print out report and terminate.
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    cat("[read.gtf.refGenome] Finished.\n")
    return(invisible())
})

setMethod("saveGenome", "refGenome",
                            function(object, filename, useBasedir=TRUE, ...)
{
    if(!is.character(filename))
        stop("'filename' must be character!")

    if(useBasedir && length(object@basedir)>0)
    {
        save(file=file.path(object@basedir, filename),
             list=ls(envir=object@ev), envir=object@ev, ...)
    }else{
        save(file=filename, list=ls(envir=object@ev), envir=object@ev, ...)
    }

    return(invisible())
})



setMethod("writeDB", "refGenome", function(object, filename, useBasedir=TRUE, ...)
{
  bm <- Sys.localeconv()[7]
  if(missing(filename))
    stop("filename argument is not optional ")

  if(useBasedir && (length(object@basedir)>0))
    filename <- file.path(object@basedir, filename)

  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, filename)

  # Write basic data
  base <- data.frame(item="class", value=class(object))
  dbWriteTable(con, "base", base, append=FALSE, overwrite=TRUE)

  # No open database connection
  if(exists("gtf", where=object@ev, inherits=FALSE))
  {
    dbWriteTable(con, "gtf", object@ev$gtf, append=FALSE, overwrite=TRUE)
    cat("[writeDB.refGenome]", format(nrow(object@ev$gtf), big.mark=bm),
                                            "rows written to table 'gtf'.\n")
  }

  if(exists("xref", where=object@ev, inherits=FALSE))
  {
    if(nrow(object@ev$xref)>0)
    {
      dbWriteTable(con, "xref", object@ev$xref, append=FALSE, overwrite=TRUE)
      cat("[writeDB.refGenome]", format(nrow(object@ev$xref), big.mark=bm),
                                        "rows written to table 'xref'.\n")
    }
  }
  dbDisconnect(con)
  return(invisible())
})

loadGenomeDb <- function(filename)
{
  if(!file.exists(filename))
    stop("File '", basename(filename), "' does not exist")
  dbdrv <- dbDriver("SQLite")
  dbcon <- dbConnect(dbdrv, filename)
  bm <- Sys.localeconv()[7]

  base <- dbReadTable(dbcon, "base")
  classitem <- match("class", base$item)

  ref <- new(base$value[classitem])
  basedir(ref) <- dirname(filename)

  copy_table <- function(con, tablename)
  {
    if(dbExistsTable(con, tablename))
    {
      assign(tablename, dbReadTable(con, tablename), envir=ref@ev)
      cat("[loadGenomeDb] ", format(nrow(ref@ev$gtf), big.mark=bm),
                                " rows copied to '", tablename, "'.\n", sep="")
    }
  }
  copy_table(dbcon, "gtf")
  copy_table(dbcon, "xref")
  copy_table(dbcon, "genes")
  copy_table(dbcon, "ujs")
  return(ref)
}

setMethod("extractByGeneId", c("refGenome", "character"), function(object, geneids, ...)
{
    if(!exists("gtf", where = object@ev, inherits = FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")

    if(is.na(match("gene_name", names(object@ev$gtf))))
        stop("'gene_name' column does not exist in 'gtf' table!")

    #if(!is.character(geneNames))
    #    stop("geneNames must be character!")

    mtc <- match(geneids, object@ev$gtf$gene_id)
    if(any(is.na(mtc)))
    {
        message("Missing matches for ", sum(is.na(mtc)),
                " gene_ID(s):\n", sep = "")

        print(geneids[is.na(mtc)])
        if(all(is.na(mtc)))
            return(invisible(NULL))
    }

    mtc <- mtc[!is.na(mtc)]
    # Retrieving gene_id's for geneids
    dtb <- data.frame(gene_name = object@ev$gtf$gene_name[mtc])

    # Returning all rows that match with found gene_name's
    gtf <- merge(object@ev$gtf, dtb)

    # Re-calibrate factor levels
    fc <- which(unlist(lapply(gtf,class)) == "factor")
    gtf[ , fc] <- data.frame(lapply(gtf[ , fc], factor))

    gtf <- gtf[order(gtf$gene_name,gtf$start), ]

    # Assemble result object
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf, envir = res@ev)
    
    if(exists("genes", envir=object@ev))
    {
        mtc <- match(geneids, object@ev$genes$gene_id)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})

setMethod("extractByGeneId", c("refJunctions", "character"), function(object, geneids, ...)
{
    if(!exists("gtf", where = object@ev, inherits = FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")
    
    if(is.na(match("gene_name", names(object@ev$gtf))))
        stop("'gene_name' column does not exist in 'gtf' table!")
    
    #if(!is.character(geneNames))
    #    stop("geneNames must be character!")
    
    mtc <- match(geneids, object@ev$gtf$gene_id)
    if(any(is.na(mtc)))
    {
        message("Missing matches for ", sum(is.na(mtc)),
                " gene_ID(s):\n", sep = "")
        
        print(geneids[is.na(mtc)])
        if(all(is.na(mtc)))
            return(invisible(NULL))
    }
    
    mtc <- mtc[!is.na(mtc)]
    # Retrieving gene_id's for geneids
    dtb <- data.frame(gene_name = object@ev$gtf$gene_name[mtc])
    
    # Returning all rows that match with found gene_name's
    gtf <- merge(object@ev$gtf, dtb)
    
    # Re-calibrate factor levels
    fc <- which(unlist(lapply(gtf,class)) == "factor")
    gtf[ , fc] <- data.frame(lapply(gtf[ , fc], factor))
    
    gtf <- gtf[order(gtf$gene_name,gtf$lend), ]
    
    # Assemble result object
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf, envir = res@ev)
    
    if(exists("genes", envir=object@ev))
    {
        mtc <- match(geneids, object@ev$genes$gene_id)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})


setMethod("extractByGeneName", c("refGenome", "character"), function(object, geneNames, ...)
{
    if(!exists("gtf", where = object@ev, inherits = FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")

    if(is.na(match("gene_name", names(object@ev$gtf))))
        stop("'gene_name' column does not exist in 'gtf' table!")

    mtc <- match(geneNames,object@ev$gtf$gene_name)
    if(any(is.na(mtc)))
    {
        message("Missing matches for ", sum(is.na(mtc)),
                                            " gene_name(s):\n", sep = "")

        print(geneNames[is.na(mtc)])
        if(all(is.na(mtc)))
            return(invisible(NULL))
    }

    mtc <- mtc[!is.na(mtc)]
    # Retrieving gene_id's for geneNames
    dtb <- data.frame(gene_name = object@ev$gtf$gene_name[mtc])

    # Returning all rows that match with found gene_name's
    gtf <- merge(object@ev$gtf, dtb)

    # Re-calibrate factor levels
    fc <- which(unlist(lapply(gtf,class)) == "factor")
    gtf[ , fc] <- data.frame(lapply(gtf[ , fc], factor))

    gtf <- gtf[order(gtf$gene_name, gtf$start), ]

    # Assemble result object
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf, envir = res@ev)
    
    if(exists("genes", envir=object@ev))
    {
        mtc <- match(geneNames, object@ev$genes$gene_name)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})

setMethod("extractByGeneName", c("refJunctions", "character"), function(object, geneNames, ...)
{
    if(!exists("gtf", where = object@ev, inherits = FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")
    
    if(is.na(match("gene_name", names(object@ev$gtf))))
        stop("'gene_name' column does not exist in 'gtf' table!")
    
    mtc <- match(geneNames,object@ev$gtf$gene_name)
    if(any(is.na(mtc)))
    {
        message("Missing matches for ", sum(is.na(mtc)),
                " gene_name(s):\n", sep = "")
        
        print(geneNames[is.na(mtc)])
        if(all(is.na(mtc)))
            return(invisible(NULL))
    }
    
    mtc <- mtc[!is.na(mtc)]
    # Retrieving gene_id's for geneNames
    dtb <- data.frame(gene_name = object@ev$gtf$gene_name[mtc])
    
    # Returning all rows that match with found gene_name's
    gtf <- merge(object@ev$gtf, dtb)
    
    # Re-calibrate factor levels
    fc <- which(unlist(lapply(gtf,class)) == "factor")
    gtf[ , fc] <- data.frame(lapply(gtf[ , fc], factor))
    
    gtf <- gtf[order(gtf$gene_name, gtf$lend), ]
    
    # Assemble result object
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf, envir = res@ev)
    
    if(exists("genes", envir=object@ev))
    {
        mtc <- match(geneNames, object@ev$genes$gene_name)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})



## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  ucscGenome                                                               ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

ucscGenome<-function(basedir)
{
  obj<-new("ucscGenome")
  if(!missing(basedir))
  {
    if(!is.character(basedir))
      stop("'basedir' must be character!")
    if(!file.exists(basedir))
      stop("'basedidr' does not exist!")
    obj@basedir<-basedir
  }
  return(obj)
}



setMethod("addIsoforms", "ucscGenome",
                function(object, filename="ucsc_knownisoforms.csv",
                                    sep="\t", useBasedir=TRUE, ...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir, filename)

  if(!exists("gtf", where=object@ev, inherits=FALSE))
    stop("No gtf table present! Use 'read.gtf' or 'setGtf'!\n")

  if(!file.exists(filename))
    stop("File '", filename, "' does not exist!\n", sep="")

  message("Reading file '", basename(filename), "'.", sep="")

  dt <- read.table(filename,  sep = sep,  quote = "",
                                stringsAsFactors = FALSE,  comment.char = "#")
  if(ncol(dt)!=2)
    stop("Wrong number of columns: ", ncol(dt), "\n")

  names(dt)<-c("clusterId", "transcript_id")

  mtc<-match(object@ev$gtf$transcript_id, dt$transcript_id)
  if(any(is.na(mtc)))
    message("[addIsoforms.ucscGenome] Warning: gtf$transcript_id misses matches in isoforms table!")

  object@ev$gtf$clusterId<-dt$clusterId[mtc]
  message("[addIsoforms.ucscGenome] Finished.")
  return(invisible())
})



setMethod("addEnsembl", "ucscGenome", function(object,
            filename="ucsc_knownToEnsembl.csv", sep="\t", useBasedir=TRUE, ...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir, filename)

  if(!exists("gtf", where=object@ev, inherits=FALSE))
    stop("No gtf table present! Use 'read.gtf' or 'setGtf'!\n")

  if(!file.exists(filename))
    stop("File '", filename, "' does not exist!\n", sep="")

  message("[addEnsembl.ucscGenome] Reading file '", basename(filename), "'.", sep="")

  dt <- read.table(filename, sep=sep, quote="", stringsAsFactors=FALSE, comment.char="#")

  if(ncol(dt)!=2)
    stop("[addEnsembl.ucscGenome] Wrong number of columns: ", ncol(dt), "\n")

  names(dt) <- c("transcript_id", "ensembl")
  mtc <- match(object@ev$gtf$transcript_id, dt$transcript_id)

  object@ev$gtf$ensembl <- factor(dt$ensembl[mtc])
  message("[addEnsembl.ucscGenome] Finished.")
  return(invisible())
})



setMethod("addXref", "ucscGenome",
                    function(object, filename = "kgXref.csv",
                                sep = "\t", useBasedir = TRUE, ...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename <- file.path(object@basedir, filename)

  if(!exists("gtf", where=object@ev, inherits=FALSE))
    stop("No gtf table present! Use 'read.gtf' or 'setGtf'!")

  if(!file.exists(filename))
    stop("File '", filename, "' does not exist!", sep="")

  message("[addXref.ucscGenome] Reading file '", basename(filename), "'.", sep="")

  dt <- read.table(filename, sep=sep, quote="", comment.char="#")
  if(ncol(dt) != 10)
    stop("Wrong number of columns: ", ncol(dt))

  names(dt) <- c("kgID", "mRNA", "spID", "spDisplayID", "gene_name",
               "refseq", "protAcc", "description", "rfamAcc", "tRnaName")

  mtc <- match(object@ev$gtf$gene_id, dt$kgID)
  object@ev$gtf$gene_name <- factor(dt$gene_name[mtc])
  assign("xref", dt, envir=object@ev)

  message("[addXref.ucscGenome] Finished.")
  return(invisible())
})


setMethod("getXref", "ucscGenome", function(object)
{
  if(!exists("xref", where = object@ev, inherits = FALSE))
    stop("No xref table present! Use 'addXref'!\n")

  return(invisible(object@ev$xref))
})



setMethod("getGenePositions", "ucscGenome", function(object, by, force=FALSE, ...)
{
  # UCSC has can have many gene_id's for one gene_name
  # Differing by 'gene_name' makes sense for UCSC
  if(missing(by))
    by <- "gene_name"

  if(!is.character(by))
    stop("[getGenePositions.ucscGenome] by must be character!")
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    return(NULL)
  if(!is.data.frame(object@ev$gtf))
    stop("[getGenePositions.ucscGenome] gtf-table must be data.frame!")

  if(by=="gene_id")
  {
    if(!is.element("gene_id", names(object@ev$gtf)))
      stop("gtf-table must contain 'gene_id' column!")

    # Min start position (table has same order as genes!)
    mig <- summaryBy(start~gene_id, data=object@ev$gtf, FUN=min)
    # Max end   position (table has same order as genes!)
    mxg <- summaryBy(end~gene_id, data=object@ev$gtf, FUN=max)

    # Get (sorted) gene_id's
    # (as.numeric(genes)==1:n! asc sorted!)
    genes <- factor(levels(object@ev$gtf$gene_id))
    n <- length(genes)
    # Point back into source table
    mtc <- match(genes, object@ev$gtf$gene_id)

    # Assemble result
    if(is.na(match("gene_name", names(object@ev$gtf))))
    {
      res <- data.frame(id=1:n, gene_id=genes,
                      seqid=object@ev$gtf$seqid[mtc], start=mig[, 2], end=mxg[, 2],
                      strand=object@ev$gtf$strand[mtc])
    } else {
      res <- data.frame(id=1:n, gene_id=genes, gene_name=object@ev$gtf$gene_name[mtc],
                      seqid=object@ev$gtf$seqid[mtc], start=mig[, 2], end=mxg[, 2],
                      strand=object@ev$gtf$strand[mtc])
    }
  }
  else if(by=="gene_name")
  {
    if(!is.element("gene_name", names(object@ev$gtf)))
      stop("gtf-table must contain 'gene_name' column. Use 'addXref'!")

    # Min start position (table has same order as genes!)
    mig <- summaryBy(start~gene_name, data=object@ev$gtf, FUN=min)
    # Max end   position (table has same order as genes!)
    mxg <- summaryBy(end~gene_name, data=object@ev$gtf, FUN=max)

    # Get (sorted) gene_id's
    # (as.numeric(genes)==1:n! asc sorted!)
    genes <- factor(levels(object@ev$gtf$gene_name))
    n <- length(genes)
    # Point back into source table
    mtc <- match(genes, object@ev$gtf$gene_name)

    res <- data.frame(id=1:n, gene_id=object@ev$gtf$gene_id[mtc], gene_name=genes,
              seqid=object@ev$gtf$seqid[mtc], start=mig[, 2], end=mxg[, 2],
              strand=object@ev$gtf$strand[mtc])
  }
  else
    stop("'by' must be 'gene_id' or 'gene_name'!")

  message("[getGenePositions.ucscGenome] Adding 'start_codon' and 'stop_codon' positions.")
  strt <- extractFeature(object, "start_codon")@ev$gtf
  mtc <- match(res$gene_id, strt$gene_id)
  stap <- strt$start[mtc]
  stam <- strt$end[mtc]
  res$start_codon <- ifelse(res$strand=='+', stap, stam)
  stpp <- extractFeature(object, "stop_codon")@ev$gtf
  mtc <- match(res$gene_id, stpp$gene_id)
  sttp <- stpp$start[mtc]
  sttm <- stpp$end[mtc]
  res$stop_codon <- ifelse(res$strand=='+', sttp, sttm)

  res <- res[order(res$seqid, res$start), ]
  assign("genes", res, envir=object@ev)
  return(invisible(res))
})

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  ensemblGenome                                                            ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

ensemblGenome <- function(basedir)
{
  obj <- new("ensemblGenome")
  if(!missing(basedir))
  {
    if(!is.character(basedir))
      stop("[ensemblGenome] basedir must be character!")
    if(!file.exists(basedir))
      stop("[ensemblGenome] basedidr does not exist!")
    obj@basedir <- basedir
  }
  return(obj)
}

setMethod("getGeneTable", "ensemblGenome", function(object)
{
    if(!exists("genes", where=object@ev, inherits=FALSE))
        return(NULL)

    return(object@ev$genes)
})

setMethod("tableFeatures", "ensemblGenome", function(object)
{return(table(object@ev$gtf$feature))})

setMethod("tableFeatures", "ucscGenome", function(object)
{return(table(object@ev$gtf$feature))})


setMethod("extractFeature", c("refGenome", "character"), function(object, feature="exon")
{
    if(!exists("gtf", where=object@ev, inherits=FALSE))
        return(NULL)
    
    if(!is.data.frame(object@ev$gtf))
        stop("gtf-table must be data.frame!")
    
    if(length(feature) > 1)
    {
        warning("[extractFeature.ensemblGenome] Only first feature is used!")
        feature <- feature[1]
    }
    
    # Returning all rows that has "CDS" feature
    ev <- new.env()
    gtf <- object@ev$gtf[object@ev$gtf$feature==feature, ]
    
    # Re-calibrate factor levels
    fc <- which(unlist(lapply(gtf, class))=="factor")
    gtf[, fc] <- data.frame(lapply(gtf[, fc], factor))
    
    # Assemble result object
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf[order(gtf$seqid, gtf$start), ], envir=res@ev)
    
    # Also copy an extract of gene-table into target
    if(exists("genes", envir=object@ev))
    {
        geneids <- unique(object@ev$gtf$gene_id)
        mtc <- match(geneids, object@ev$genes$gene_id)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})


setMethod("extractTranscript", c("refGenome", "character"), function(object, transcripts)
{
    if(!exists("gtf", where=object@ev, inherits=FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")
    
    if(is.na(match("transcript_id", names(object@ev$gtf))))
        stop("'transcript_id' column does not exist in 'gtf' table!")
    
    mtc <- match(transcripts, object@ev$gtf$transcript_id)
    if(any(is.na(mtc)))
    {
        cat("[extractTranscript.ensemblGenome] Missing matches for ", sum(is.na(mtc)), " transcripts:\n", sep="")
        print(transcripts[is.na(mtc)])
    }
    # reorder (transcrpt_id)
    dtb <- data.frame(transcript_id=transcripts)
    
    # Extract and re-calibrate factor levels
    gtf <- merge(object@ev$gtf, dtb)
    fc <- which(unlist(lapply(gtf, class))=="factor")
    gtf[, fc] <- data.frame(lapply(gtf[, fc], factor))
    
    res <- new(class(object))
    basedir(res) <- basedir(object)
    assign("gtf", gtf, envir=res@ev)
    
    # Also copy an extract of gene-table into target
    if(exists("genes", envir=object@ev))
    {
        geneids <- unique(object@ev$gtf$gene_id)
        mtc <- match(geneids, object@ev$genes$gene_id)
        assign("genes", object@ev$genes[mtc, ], envir=res@ev)
    }
    
    return(res)
})


setMethod("tableTranscript.id", "ensemblGenome", function(object)
{
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    stop("gtf table does not exist! Use 'read.gtf'.")

  if(is.na(match("transcript_id", names(object@ev$gtf))))
    stop("gtf does not contain column 'transcript_id'!")

  return(table(object@ev$gtf$transcript_id))
})
setMethod("tableTranscript.id", "ucscGenome", function(object)
{
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    stop("gtf table does not exist! Use 'read.gtf'.")

  if(is.na(match("transcript_id", names(object@ev$gtf))))
    stop("gtf does not contain column 'transcript_id'.")

  return(table(object@ev$gtf$transcript_id))
})


setMethod("tableTranscript.name", "ensemblGenome", function(object)
{
    if(!exists("gtf", where=object@ev, inherits=FALSE))
        stop("gtf table does not exist! Use 'read.gtf'.")
    
    if(is.na(match("transcript_name", names(object@ev$gtf))))
        stop("gtf table does not contain column 'transcript_name'!")
    
    return(table(object@ev$gtf$transcript_name))
})



setMethod("getGenePositions", c("ensemblGenome","character"), function(object, by, force = FALSE, ...)
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(!is.logical(force))
        stop("[getGenePositions.ensemblGenome] force must be logical!")
    
    # Copy of table will be in ev -> positions
    # need only once be calculated.
    if(exists("genes", where=object@ev, inherits=FALSE) & !force)
        return(object@ev$genes)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Differing by 'gene_id' makes sense for Ensembl
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(missing(by))
        by <- "gene_id"
    if(!is.character(by))
        stop("'by' must be character!")
    
    
    if(!exists("gtf", where=object@ev, inherits=FALSE))
        return(NULL)
    if(!is.data.frame(object@ev$gtf))
        stop("gtf-table must be data.frame!")
    
    if(is.na(match("gene_name", names(object@ev$gtf))))
        stop("No 'gene_name' data found!")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Different gene-identifications
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(by=="gene_id")
    {
        # Get (sorted) gene_id's
        # (as.numeric(genes)==1:n! asc sorted!)
        genes <- factor(levels(object@ev$gtf$gene_id))
        n <- length(genes)
        # Point back into source table
        mtc <- match(genes, object@ev$gtf$gene_id)
        
        # Min start position (table has same order as genes!)
        mig <- summaryBy(start~gene_id, data=object@ev$gtf, FUN=min)
        # Max end   position (table has same order as genes!)
        mxg <- summaryBy(end~gene_id, data=object@ev$gtf, FUN=max)
    }
    else if(by=="gene_name")
    {
        # Get (sorted) gene_id's
        # (as.numeric(genes)==1:n! asc sorted!)
        genes <- factor(levels(object@ev$gtf$gene_name))
        n <- length(genes)
        # Point back into source table
        mtc <- match(genes, object@ev$gtf$gene_name)
        
        # Min start position (table has same order as genes!)
        mig <- summaryBy(start~gene_name, data=object@ev$gtf, FUN=min)
        # Max end   position (table has same order as genes!)
        mxg <- summaryBy(end~gene_name, data=object@ev$gtf, FUN=max)
    }
    else
        stop("[getGenePositions.ensemblGenome] by must be 'gene_id' or 'gene_name'!")
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Assemble result
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(is.na(match("gene_biotype", names(object@ev$gtf))))
    {
        res <- data.frame(id = 1:n,  gene_id = object@ev$gtf$gene_id[mtc],
                          gene_name = object@ev$gtf$gene_name[mtc],
                          seqid = object@ev$gtf$seqid[mtc],
                          start = mig[, 2],
                          end = mxg[, 2],
                          strand = object@ev$gtf$strand[mtc])
        
    }else{
        res <- data.frame(id = 1:n, gene_id = object@ev$gtf$gene_id[mtc],
                          gene_name = object@ev$gtf$gene_name[mtc],
                          seqid = object@ev$gtf$seqid[mtc],
                          start = mig[, 2],
                          end=mxg[, 2],
                          strand = object@ev$gtf$strand[mtc],
                          gene_biotype = object@ev$gtf$gene_biotype[mtc])
    }
    
    message("[getGenePositions.ensemblGenome] Adding 'start_codon' and 'stop_codon' positions.")
    strt <- extractFeature(object, "start_codon")@ev$gtf
    mtc <- match(res$gene_id, strt$gene_id)
    stap <- strt$start[mtc]
    stam <- strt$end[mtc]
    res$start_codon <- ifelse(res$strand=='+', stap, stam)
    stpp <- extractFeature(object, "stop_codon")@ev$gtf
    mtc <- match(res$gene_id, stpp$gene_id)
    sttp <- stpp$start[mtc]
    sttm <- stpp$end[mtc]
    res$stop_codon <- ifelse(res$strand=='+', sttp, sttm)
    
    res <- res[order(res$seqid, res$start), ]
    assign("genes", res, envir=object@ev)
    return(invisible(res))
})

setMethod("getGenePositions", c("ensemblGenome","missing"), function(object, by, force = FALSE, ...)
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Differing by 'gene_id' makes sense for Ensembl
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    return(getGenePositions(object=object, by="gene_id", force=force, ...))
})



## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##  tableSeqids and extract seqids, targeted  on regex:
##  Primary intention is to retreave the primary assembly and remove haplotypes
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##



setMethod("tableSeqids", "refGenome", function(object, regex)
{
  if(missing(regex))
    return(table(object@ev$gtf$seqid))
  else
  {
    if(!is.character(regex))
      stop("[tableSeqids.refGenome] regex must be character!")
    tbl <- table(object@ev$gtf$seqid)
    seqNames <- names(tbl)

    return(data.frame(seqid = seqNames,  nAnnotations = as.numeric(tbl),
                    lgl = grepl(regex, seqNames)))
  }
})


setMethod("extractSeqids", "refGenome", function(object, regex)
{
  if(!is.character(regex))
    stop("[extractSeqids.refGenome] regex must be character!")
  lgl <- grepl(regex, object@ev$gtf$seqid)
  ans <- new(class(object))
  basedir(ans) <- basedir(object)

  # Extract and remove unused factor levels
  gtf <- object@ev$gtf[lgl, ]
  fc <- which(unlist(lapply(gtf, class))=="factor")
  gtf[, fc] <- data.frame(lapply(gtf[, fc], factor))
  ans@ev$gtf <- gtf

  if(exists("xref", where=object@ev, inherits=FALSE))
  {
    # Should only be present in uscsGenome
    tbl <- data.frame(kgID=unique(object@ev$gtf$gene_id))
    ans@ev$xref <- merge(object@ev$xref, tbl)
  }
  if(exists("genes", where=object@ev, inherits=FALSE))
  {
    # Should only be present in ensemblGenome
    tbl <- object@ev$genes
    lgl <- grepl(regex, tbl$seqid)
    tbl <- tbl[lgl, ]

    # Extract and remove unused factor levels
    fc <- which(unlist(lapply(tbl, class))=="factor")
    tbl[, fc] <- data.frame(lapply(tbl[, fc], factor))
    ans@ev$genes <- tbl
  }
  return(ans)
})

setMethod("extractPaGenes", "ensemblGenome", function(object)
{
  if(is.na(match("gene_name", names(object@ev$gtf))))
    stop("gtf table does not contain column 'gene_name'!")
  if(is.na(match("exon_number", names(object@ev$gtf))))
    stop("gtf table does not contain column 'exon_number'!")

  enpa <- extractSeqids(object, ensPrimAssembly())
  return(getGenePositions(enpa))
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## ref-Exons
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

ensemblExons <- function(basedir)
{
    obj <- new("ensemblExons")
    if(!missing(basedir))
    {
      if(!is.character(basedir))
        stop("[ensemblExons] basedir must be character!")
      if(!file.exists(basedir))
        stop("[ensemblExons] basedidr does not exist!")
      obj@basedir <- basedir
    }
    return(obj)
}

ucscExons <- function(basedir)
{
  obj <- new("ucscExons")
  if(!missing(basedir))
  {
    if(!is.character(basedir))
      stop("[ucscExons] basedir must be character!")
    if(!file.exists(basedir))
      stop("[ucscExons] basedidr does not exist!")
     obj@basedir <- basedir
  }
  return(obj)
}


setMethod("addExonNumber", "refGenome", function(object)
{
    dfr <- data.frame(tr = as.integer(object@ev$gtf$transcript_id),
                    sq = as.integer(object@ev$gtf$seqid),
                    st = as.integer(object@ev$gtf$start),
                    en = as.integer(object@ev$gtf$end),
                    id = object@ev$gtf$id
                    )

    dfr <- dfr[order(dfr$tr, dfr$sq, dfr$st, dfr$en), ]

    res <- .Call("get_exon_number", dfr$tr, dfr$sq, dfr$st,
                                            dfr$en, PACKAGE = "refGenome")

    mtc <- match(dfr$id,object@ev$gtf$id)
    object@ev$gtf$exon_number <- res[mtc]
    return(invisible())
})


setMethod("refExons", "refGenome", function(object)
{
  if(!is.element("exon_number", names(object@ev$gtf)))
     addExonNumber(object)

  cat("[refExons.refGenome] Extracting tables.\n")
  cds <- extractFeature(object, "CDS")@ev$gtf
  exons <- extractFeature(object, "exon")@ev$gtf
  stacod <- extractFeature(object, "start_codon")@ev$gtf
  stocod <- extractFeature(object, "stop_codon")@ev$gtf

  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ## CDS
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  cat("[refExons.refGenome] Adding 'CDS'.\n")
  cds <- cds[, c("start", "end", "transcript_id", "exon_number")]
  names(cds)[1:2] <- c("cds_start", "cds_end")
  ex_cds <- merge(exons, cds, by=c("transcript_id", "exon_number"), all=T)
  ex_cds$cds_start <- ex_cds$cds_start-ex_cds$start
  ex_cds$cds_end <- ex_cds$end-ex_cds$cds_end

  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ## start_codon
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  cat("[refExons.refGenome] Adding 'start_codon'.\n")
  start_codons <- stacod[, c("start", "transcript_id", "exon_number")]
  names(start_codons)[1] <- "start_codon"
  ex_cds_st <- merge(ex_cds, start_codons, by=c("transcript_id", "exon_number"), all=T)
  ex_cds_st$start_codon <- ex_cds_st$start_codon-ex_cds_st$start

  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  ## stop_codon
  ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  cat("[refExons.refGenome] Adding 'stop_codon'.\n")
  stop_codons <- stocod[, c("start", "transcript_id", "exon_number")]
  names(stop_codons)[1] <- "stop_codon"
  ex_cds_ststp <- merge(ex_cds_st, stop_codons, by=c("transcript_id", "exon_number"), all=T)
  ex_cds_ststp$stop_codon <- ex_cds_ststp$stop_codon-ex_cds_ststp$start

  nCol <- ncol(exons)
  mtc <- match(names(exons), names(ex_cds_ststp)[1:nCol])
  ex_cds_ststp <- ex_cds_ststp[, c(mtc, nCol+1:4)]
  ex_cds_ststp <- ex_cds_ststp[order(ex_cds_ststp$seqid, ex_cds_ststp$start), ]
  ex_cds_ststp$exon_number <- as.numeric(ex_cds_ststp$exon_number)

  # remove feature column (contains only "exon" entries)
  ex_cds_ststp$feature <- NULL

  if(is(object, "ensemblGenome"))
    res <- ensemblExons(basedir(object))
  else
    res <- ucscExons(basedir(object))

  res@ev$gtf <- ex_cds_ststp
  #assign("gtf", ex_cds_ststp, envir=res@ev)

  # Should only be present in ensemblGenomes
  if(exists("genes", where=object@ev, inherits=FALSE))
      res@ev$genes <- object@ev$genes

  cat("[refExons.refGenome] Finished.\n")
  return(res)
})


ensemblJunctions <- function(basedir)
{
  obj <- new("ensemblJunctions")
  if(!missing(basedir))
  {
    if(!is.character(basedir))
      stop("[ensemblJunctions] basedir must be character!")
    if(!file.exists(basedir))
      stop("[ensemblJunctions] basedidr does not exist!")
    obj@basedir <- basedir
  }
  return(obj)
}

ucscJunctions <- function(basedir)
{
  obj <- new("ucscJunctions")
  if(!missing(basedir))
  {
    if(!is.character(basedir))
      stop("[ucscJunctions] basedir must be character!")
    if(!file.exists(basedir))
      stop("[ucscJunctions] basedidr does not exist!")
    obj@basedir <- basedir
  }
  return(obj)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Unified ranges are intended to be used as preparation of exonic coordinates
# for counting of alignments (reads) inside these regions.
# In order to prevent ambiguities, the function removes overlaps:
# A) Ranges which are contained in a precedent range are ignored
# B) For ranges which overlap wich a precedent range, the start position
#       is shifted to one position after the end of the precedent range.
#       The size of the shift is contained in the 'ubs' value.
#
# For an alternative approach, see:
# http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("unifyRanges", function(object) standardGeneric("unifyRanges"))
setMethod("unifyRanges", "refExons", function(object){

    ex <- object@ev$gtf
    uex <- .Call("unify_genomic_ranges",
                    ex$id, ex$seqid, ex$start, ex$end, PACKAGE="refGenome")


    # Unifying ranges may produce ranges of size 0 (or < 0):
    egb <- uex$end > uex$begin
    nrem <- sum(!egb)
    if(nrem > 0)
    {
        cat("[unifyRanges] Removing", nrem, " ranges due to size zero (or < 0).\n" )
        uex <- uex[egb, ]
    }

    # Create returned object
    res <- new("unifiedExons")
    res@basedir <- object@basedir
    uex$seqid <- factor(uex$seqid)
    levels(uex$seqid) <- levels(object@ev$gtf$seqid)


    # Create annotation table
    res@ev$gtf <- data.frame(id=uex$id, seqid=uex$seqid,
                             begin=uex$begin, end=uex$end)

    # Additional columns
    mtc <- match(res@ev$gtf$id, ex$id)
    res@ev$gtf$strand <- ex$strand[mtc]
    res@ev$gtf$gene_id <- ex$gene_id[mtc]
    res@ev$gtf$transcript_id <- ex$transcript_id[mtc]
    res@ev$gtf$gene_name <- ex$gene_name[mtc]
    res@ev$gtf$exon_number <- ex$exon_number[mtc]
    res@ev$gtf$ubs <- uex$ubs
    # (Eventually add some more columns...)

    return(res)
})


setGeneric("plotUbs", function(object, ...) standardGeneric("plotUbs"))
setMethod("plotUbs", "unifiedExons", function(object,
        xlim=c(0,2000), ylim=c(0,300),
        lwd=2, col="darkslateblue", bty="n", las=1,
        main="ubs (unification begin shift)", ...)
{
    ubs <- object@ev$gtf$ubs
    ubs <- ubs[ubs!=0]
    tbl <- table(ubs)
    xval <- as.numeric(names(tbl))
    yval <- as.numeric(tbl)
    plot(xval, yval, type="l", xlim=xlim, ylim=ylim, las=las, bty=bty,
         xlab="value", ylab="frequency", lwd=lwd, col=col, main=main, ...)
    return(invisible())
})


setMethod("getSpliceTable", "refGenome", function(object)
{
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Prepare input table
    # Filter for exons
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(is(object, "refExons"))
        gtf <- object@ev$gtf
    else
        gtf <- object@ev$gtf[object@ev$gtf$feature=="exon", ]

    # Shape and sort input values
    gtf <- gtf[!is.na(gtf$transcript_id), ]
    dfr <- data.frame(tr=as.integer(factor(gtf$transcript_id)),
                    sq=as.integer(factor(gtf$seqid)),
                    st=as.integer(gtf$start),
                    en=as.integer(gtf$end),
                    id=gtf$id
    )
    dfr <- dfr[order(dfr$tr, dfr$sq, dfr$st), ]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Junction assembly in C
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    res  <-  .Call("get_splice_juncs",  dfr$tr,  dfr$id,  dfr$st,
                                                dfr$en, PACKAGE="refGenome")


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Assembly of result table
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    n <- nrow(res)
    mtc <- match(res$lexid, object@ev$gtf$id)

    junc <- data.frame(id=1:n, seqid=factor(object@ev$gtf$seqid[mtc]),
                   lstart=res$lstart, lend=res$lend,
                   rstart=res$rstart, rend=res$rend,
                   gene_id=factor(object@ev$gtf$gene_id[mtc]),
                   gene_name=factor(object@ev$gtf$gene_name[mtc]),
                   strand=factor(object@ev$gtf$strand[mtc]),
                   transcript_id=factor(object@ev$gtf$transcript_id[mtc]),
                   lexid=res$lexid, rexid=res$rexid)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Eventually add column (used in unifyJuncs for determination of
    # NMD = nonsense mediated decay)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(!is.na(match("transcript_biotype", names(object@ev$gtf))))
        junc$transcript_biotype <- object@ev$gtf$transcript_biotype[mtc]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Construct instance of returned class
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(is(object, "ensemblGenome")|| is(object, "ensemblExons"))
        res <- ensemblJunctions(basedir(object))
    else
        res <- ucscJunctions(basedir(object))

    res@ev$gtf <- junc

    # Should only be present in ensemblGenomes
    if(exists("genes", where=object@ev, inherits=FALSE))
        res@ev$genes <- object@ev$genes


    cat("[getSpliceTable.refGenome] Finished.\n")
    return(res)
})


setMethod("unifyJuncs", "refJunctions", function(object){

    # Create shortcut reference
    gtf <- object@ev$gtf

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Check whether feature is not of type NMD (nonsense mediated decay).
    # This will later be used to determine whether a junction
    # is present in any non-NMD Transcript at all.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    mtc <- match("transcript_biotype", names(gtf))
    if(!is.na(mtc))
    {
        nnmd <- as.integer(gtf$transcript_biotype != "nonsense_mediated_decay")
    }else{
        nnmd <- rep(0L, nrow(gtf))
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Prepare input table
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    dfr <- data.frame(seqid = as.integer(gtf$seqid),
                lstart = gtf$lstart,
                lend = gtf$lend,
                rstart = gtf$rstart,
                rend = gtf$rend,
                id = gtf$id,
                gene_id = as.integer(gtf$gene_id),
                strand = as.integer(gtf$strand),
                nnmd = nnmd)

    dfr <- dfr[order(dfr$seqid, dfr$lend, dfr$rstart), ]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Junction unifications in C
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    ujc <- .Call("unify_splice_juncs", dfr$seqid, dfr$lstart, dfr$lend,
                dfr$rstart, dfr$rend,
                dfr$id, dfr$gene_id, dfr$strand, dfr$nnmd, PACKAGE="refGenome")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Assembly of result table
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Create factors from numeric values and remove unused levels
    seql <- 1:length(levels(object@ev$gtf$seqid))
    genl <- 1:length(levels(object@ev$gtf$gene_id))
    strl <- 1:length(levels(object@ev$gtf$strand))

    ujs <- data.frame(id=ujc$id,
            seqid = factor(ujc$seqid,  levels = seql,
                                labels = levels(object@ev$gtf$seqid)),
            lstart = ujc$lstart,
            lend = ujc$lend,
            rstart = ujc$rstart,
            rend = ujc$rend,
            nSites = ujc$nSites,
            gene_id = factor(ujc$gene_id, levels = genl,
                                        labels = levels(object@ev$gtf$gene_id)),
            strand = factor(ujc$strand,  levels=strl,
                                        labels = levels(object@ev$gtf$strand)),
            fexid = ujc$fexid,
            cnNmd=ujc$cnnmd)

    # Add gene name
    mtc <- match(ujs$gene_id, gtf$gene_id)
    ujs$gene_name = factor(gtf$gene_name[mtc])

    # Write a copy of the table back into source object
    assign("ujs", ujs, envir=object@ev)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Construct instance of returned class
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

    res <- new(class(object))
    basedir(res) <- basedir(object)
    res@ev$gtf <- ujs

    return(res)
})


setMethod("getGenePositions", "refJunctions",
                                    function(object, by, force = FALSE, ...)
{
    # Works the same way as version for ensemblGenome
    # Only start -> lstart, end -> rend changed.

  # + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  if(!is.logical(force))
    stop("[getGenePositions.refJunctions] force must be logical!")

  # Copy of table will be in ev -> positions
  # need only once be calculated.
  if(exists("genes",where=object@ev,inherits=FALSE) & !force)
    return(object@ev$genes)

  # + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  # Differing by 'gene_id' makes sense for Ensembl
  if(missing(by))
    by <- "gene_id"
  if(!is.character(by))
    stop("[getGenePositions.refJunctions] by must be character!")

  # + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  if(!exists("gtf", where=object@ev, inherits=FALSE))
    return(NULL)
  if(!is.data.frame(object@ev$gtf))
    stop("[getGenePositions.refJunctions] gtf-table must be data.frame!")
  if(is.na(match("gene_name", names(object@ev$gtf))))
    stop("[getGenePositions.refJunctions] No 'gene_name' data found!")

  # + + + + + + + + + + + + + + + + + + + + + + + + + + + +
  if(by=="gene_id")
  {
    if(!is.element("gene_id", names(object@ev$gtf)))
      stop("gtf-table must contain 'gene_id' column!")
    # Get (sorted) gene_id's
    # (as.numeric(genes)==1:n! asc sorted!)
    genes <- factor(levels(object@ev$gtf$gene_id))
    n <- length(genes)
    # Point back into source table
    mtc <- match(genes, object@ev$gtf$gene_id)

    # Min start position (table has same order as genes!)
    mig <- summaryBy(lstart~gene_id, data=object@ev$gtf, FUN=min)
    # Max end   position (table has same order as genes!)
    mxg <- summaryBy(rend~gene_id, data=object@ev$gtf, FUN=max)

    # Assemble result
    res <- data.frame(id=1:n, gene_id=genes,
                    seqid=object@ev$gtf$seqid[mtc], start=mig[, 2], end=mxg[, 2],
                    strand=object@ev$gtf$strand[mtc])
  }
  else if(by=="gene_name")
  {
    if(!is.element("gene_id", names(object@ev$gtf)))
      stop("gtf-table must contain 'gene_id' column!")
    # Get (sorted) gene_id's
    # (as.numeric(genes)==1:n! asc sorted!)
    genes <- factor(levels(object@ev$gtf$gene_name))
    n <- length(genes)
    # Point back into source table
    mtc <- match(genes, object@ev$gtf$gene_name)

    # Min start position (table has same order as genes!)
    mig <- summaryBy(lstart~gene_name, data=object@ev$gtf, FUN=min)
    # Max end   position (table has same order as genes!)
    mxg <- summaryBy(rend~gene_name, data=object@ev$gtf, FUN=max)

    # Assemble result
    res <- data.frame(id=1:n, gene_name=genes,
                    seqid=object@ev$gtf$seqid[mtc], start=mig[, 2], end=mxg[, 2],
                    strand=object@ev$gtf$strand[mtc])
  }
  else
    stop("'by' must be 'gene_id' or 'gene_name'!")

  res <- res[order(res$seqid, res$start), ]

  assign("genes", res, envir=object@ev)
  return(invisible(res))
})


setMethod("addIsCoding", "ensemblJunctions", function(object, ens)
{
  if(!is(ens, "ensemblGenome"))
    stop("[addIsCoding.ensemblJunctions] ens must be 'ensemblGenome'!")

  cds <- extractFeature(ens, "CDS")@ev$gtf
  if(nrow(cds)==0)
    stop("[addIsCoding.ensemblJunctions] ens contains no 'CDS' entries!")

  lc <- cds[, c("end", "transcript_id", "feature")]
  names(lc)[1] <- "lend"
  rc <- cds[, c("start", "transcript_id", "feature")]
  names(rc)[1] <- "rstart"
  jc <- object@ev$gtf[order(object@ev$gtf$transcript_id, object@ev$gtf$lend),
            c("id", "lend", "rstart", "transcript_id")]

  message("[addIsCoding.ensemblJunctions] Adding left  coding.")
  lcd <- merge(jc, lc, all.x=TRUE)
  mtc <- match(object@ev$gtf$id, lcd$id)
  object@ev$gtf$licd <- ifelse(is.na(lcd$feature[mtc]), FALSE, TRUE)

  message("[addIsCoding.ensemblJunctions] Adding right coding.")
  rcd <- merge(jc, rc, all.x=TRUE)
  mtc <- match(object@ev$gtf$id, rcd$id)
  object@ev$gtf$ricd <- ifelse(is.na(rcd$feature[mtc]), FALSE, TRUE)
  message("[addIsCoding.ensemblJunctions] Finished.")
  return(invisible())
})



## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
## overlap :    Takes tables of query and reference ranges,                  ##
##              overlaps both and reports overlaps (+ misses)                ##
##              for every query range ()                                     ##
##                                                                           ##
##            Results encoding                                               ##
##            no      OVERLAP_NO_OVER     1     query misses reference       ##
##            r       OVERLAP_R_OVER      2     query overhangs on           ##
##                                                                right side ##
##                                                                           ##
##            b       OVERLAP_B_OVER      3     query overhangs on           ##
##                                              both sides (right and left)  ##
##                                                                           ##
##            n       OVERLAP_N_OVER      4     query overhangs on neither   ##
##                                                                  side     ##
##            l       OVERLAP_L_OVER      5     query overhangs on left side ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

overlap <- function(qry, ref)
{
    res <- .Call("overlap_ranges", as.integer(qry$id), as.integer(qry$start),
                as.integer(qry$end), as.integer(ref$id), as.integer(ref$start),
                as.integer(ref$end), PACKAGE = "refGenome")
  return(res)
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Unexported function
## Does overlapping for all present seqid's separately
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

overlap_seqs <- function(qry, ref, verbose = FALSE)
{
    ## + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## expects qry and ref as data.frames with columns
    ## id, start, end, seqid
    ## + + + + + + + + + + + + + + + + + + + + + + + + + ##

    qRefNames <- levels(qry$seqid)
    nQryRefs <- length(qRefNames)

    l <- list()
    k <- 1
    bm <- Sys.localeconv()[7]
    for(i in 1:nQryRefs)
    {
        qrs <- qry[qry$seqid == qRefNames[i], ]
        qrs <- qrs[order(qrs$start), ]
        rfs <- ref[ref$seqid == qRefNames[i], ]
        rfs <- rfs[order(rfs$start), ]

        nq<-nrow(qrs)
        nr<-nrow(rfs)
        if(verbose)
        {
            message("[overlap_genome] i: ", format(i, width = 2), "\tseqid: ",
                    format(qRefNames[i], width = 4), "\tquery set: ",
                    format(nq, width = 7, big.mark = bm),
                    "\tref set:", format(nr, width = 6, big.mark = bm))

            if(nq == 0 | nr == 0)
                message("\t\tskipped!")
            else
            {
                message("")
                l[[k]] <- overlap(qry = qrs, ref = rfs)
                k <- k + 1
            }
        }
        else # verbose=FALSE
        {
            if(nq > 0 & nr > 0)
            {
                l[[k]] <- overlap(qry = qrs, ref = rfs)
                k <- k + 1
            }
        }
    }
    return(do.call(rbind, l))
}


overlapJuncs <- function(qry, junc)
{
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Section 0: Prerequisites
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ##  qry     :
    ##  class   :   data.frame
    ##  columns :   id, seqid, lstart, lend, rstart, rend
    ##
    ##  junc    :
    ##  class   :   ensemblJunctions (or ucscJunctions?)
    ##
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

    # Static declaration
    bm <- Sys.localeconv()[7]

    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Section 1: Preparation
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

    message("[overlapJuncs] Input preparation.")
    column_names <- c("id","lstart","lend","rstart","rend","seqid")

    ## + + + + + + + + + + + + + + + ##
    ## Prepare given query data
    ## + + + + + + + + + + + + + + + ##
    mtc <- match(column_names, names(qry))
    if(any(is.na(mtc)))
        stop("'qry' table must contain names: id, seqid,
                                            lstart, lend, rstart, rend")

    # Extract qry in expected column order
    qry <- qry[ , mtc]

    # Numeric colums in qry are all integer?
    if(!all(unlist(lapply(qry[ ,1:5], is.integer))))
        stop("All numeric columns in 'qry' must be integer.")


    ## + + + + + + + + + + + + + + + ##
    ## Prepare given annotation
    ## + + + + + + + + + + + + + + + ##
    if(!is(junc, "refJunctions"))
        stop("'junc' must be of class 'refJunctions'")

    # Extract table of junction coordinates
    ref <- junc@ev$gtf[ , column_names]


    message("[overlapJuncs] Query size          :",
                format(nrow(qry), big.mark = bm, width = 10),
                "\tRef  size : ",
                format(nrow(ref), big.mark = bm, width = 10)
    )
    message("[overlapJuncs] - - - - - - - - - - - - - - - - - -",
                            " - - - - - - - - - - - -")


    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Section 3:   Do query junction table separately for
    ##              each sequence
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

    qRefNames <- levels(qry$seqid)
    nQryRefs  <- length(qRefNames)

    if(all(is.na(match(qRefNames, levels(ref$seqid)))))
        stop("No matches between query and reference seqid's (wrong genome?)")

    l <- list()
    k <- 1

    for(i in 1:nQryRefs)
    {
        qrs <- qry[qry$seqid == qRefNames[i], 1:5]
        qrs <- qrs[order(qrs$lstart, qrs$rend), ]
        rfs <- ref[ref$seqid == qRefNames[i], 1:5]
        rfs <- rfs[order(rfs$lstart, rfs$rend), ]

        nq <- nrow(qrs)
        nr <- nrow(rfs)

        message("[overlapJuncs] i: ", format(i, width = 2), "\tseqid: ",
                format(qRefNames[i], width = 4), "\tquery set: ",
                format(nq, width = 7, big.mark = bm),
                "\tref set:", format(nr, width = 6, big.mark = bm))

        rfs$rMaxRend <-.Call("get_cum_max", rfs$rend, PACKAGE = "refGenome")
        if(nq > 0 & nr > 0)
        {
            l[[k]] <- .Call("gap_overlap", qrs$id, qrs$lstart, qrs$lend,
                        qrs$rstart, qrs$rend,
                        rfs$id, rfs$lstart, rfs$lend,
                        rfs$rstart, rfs$rend,
                        rfs$rMaxRend, PACKAGE = "refGenome")
            k <- k + 1
        }else{
            message("\t\tskipped!")
        }
    }

    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Section 4: Assemble result.
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    res <- do.call(rbind, l)

    message("[overlapJuncs] Result assemly.")

    # Columns
    mtc <- match(res$refid,junc@ev$gtf$id)
    res$strand <- junc@ev$gtf$strand[mtc]
    res$gene_id <- junc@ev$gtf$gene_id[mtc]
    res$transcript_id <- junc@ev$gtf$transcript_id[mtc]


    # Check for presence of gene_names and eventually add column:
    if(is.na(match("gene_name",names(junc@ev$gtf))))
    {
        message("[gap_overlap] 'gene_name' column is not present.")
        if(!is(junc,"ensemblJunctions"))
            message("[gap_overlap] Import 'kgXref' table.")
    }else{
        res$gene_name <- junc@ev$gtf$gene_name[mtc]
    }

    # Check for presence of cnNmd (from unifyJuncs)
    if(!is.na(match("cnNmd", names(junc@ev$gtf))))
    {
        res$nSites <- junc@ev$gtf$nSites[mtc]
        res$cnNmd <- junc@ev$gtf$cnNmd[mtc]
    }

    # Transform strand values in order to apply 'strandlevels' as factor levels
    nstrand <- as.numeric(res$strand)
    nstrand[is.na(nstrand)] <- 0
    nstrand <- 3 - nstrand
    res$strand <- factor(nstrand, levels=1:3, labels = strandlevels)


    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Finish.
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

    message("[overlapJuncs] Finished.")
    return(invisible(res))
}

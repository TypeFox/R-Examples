"odfWeave" <- 
function(file, dest, workDir=odfTmpDir(), control=odfWeaveControl())
{ 

   # check that we can write to the target path
   if(file.exists(dest))
   {
      if(file.access(dest, mode = 2) < 0) stop(paste("cannot write to", dest))
   }
   
   # configure
   currentLoc <- getwd()
   zipCmd <- control$zipCmd
   zipCmd <- gsub("$$file$$", shellQuote(basename(file)), zipCmd, fixed=TRUE)
   
   dest <- canonicalFilePath(dest)
   
   verbose <- control$verbose

   # check for an unzipping utility
   if(all(zipCmd == c("zip -r $$file$$ .", "unzip -o $$file$$")))
   {
      errorText <- paste(
         "unzip could not be found.",
         "If installed, check your path.",
         "If not installed, either go to",
         "www.info-zip.org and install or",
         "use another utility (like jar)")
      if(.Platform$OS.type == "windows")
      {
         zipTest <- class(try(system("unzip", intern=TRUE, invisible=TRUE),
            silent=TRUE))
         if(class(zipTest) == "try-error") stop(errorText)
      } else {
         zipTest <- system("unzip", intern = TRUE)[1]
         if(is.na(zipTest) || length(grep("UnZip", zipTest)) == 0)
            stop(errorText)
      }
   }

   # Create temp dir to work in.  Remember if we created it so
   # we will know whether to delete it or not.
   if(!file.exists(workDir))
   {
      announce(verbose, "  Creating ", workDir, "\n")
      dir.create(workDir, showWarnings = TRUE, recursive = FALSE)
      if(!file.exists(workDir)) stop("Error creating working directory")
      created.workDir <- TRUE
   } else {
      warning("working in existing directory")
      created.workDir <- FALSE
   }

   # Make absolute version of workDir so it will be valid after
   # we cd to it
   workDir <- canonicalDirPath(workDir)

   workingCopy <- basename(file)
   
   # copy the user's ODT file to the tmp dir
   if(!file.exists(file)) stop(paste(file, "does not exist"))
   announce(verbose, "  Copying ", file, "\n")
   if(!file.copy(file, paste(workDir, "/", workingCopy, sep = ""), overwrite = TRUE)) stop("Error copying odt file")
   
   announce(verbose, "  Setting wd to ", workDir, "\n")
   setwd(workDir)
   on.exit(setwd(currentLoc))

   # unzip the copied ODT file into tmp dir
   announce(verbose, "  Unzipping ODF file using", zipCmd[2], "\n")
   if(.Platform$OS.type == "windows")
   {
      if(system(zipCmd[2], invisible = TRUE) != 0) stop("Error unzipping file")
   } else {
      if(system(zipCmd[2]) != 0) stop("Error unzipping odt file")
   }

   # remove copy of ODT file now that it is unzipped
   announce(verbose, "\n  Removing ", workingCopy, "\n")
   file.remove(workingCopy)
   if (file.exists(workingCopy)) stop("Error removing original file")

   # Make sure we can convert content.xml to the current locale,
   # otherwise the output from the 'Sweave' function will be invalid XML.
   text <- iconv(readLines('content.xml', warn=FALSE), 'UTF-8', '')
   if (any(is.na(text)))
   {
       announce(verbose, "\n  Unable to convert", workingCopy, "to the current locale\n")
       announce(verbose, "  You may need to process this file in a UTF-8 locale\n")
       stop(sprintf("unable to convert %s to the current locale",
                    workingCopy))
   }
   rm(text)

   # create Pictures directory if it was not created by unzipping the ODT file
   if(!file.exists(paste(workDir, "/Pictures", sep = "")))
   {
      announce(verbose, "  Creating a Pictures directory\n")
      picDir <- dir.create("Pictures", showWarnings = TRUE, recursive = FALSE)
      if(!picDir)  stop("Error creating Pictures directory")
   }
   assign(
      "picPath",
      paste(workDir, "/Pictures", sep = ""),
      pos = .odfEnv)

   # create a character vector in .odfEnv that will collect the
   # names of the image files inserted into the document via
   # odfInsertPlot.  This is needed so that we can add those
   # files to META-INF/manifest.xml later in this function.
   assign(
      "picVector",
      character(),
      pos = .odfEnv)
 
   # Create the "Style Name Environment" which will be used to create
   # unique style names during the Sweave phase.  Names of the existing
   # style definitions will be put into this environment during the
   # preprocessing stage.
   # Note that this assignment may overwrite an environment created on
   # an earlier call to odfWeave.  We need a new, clean one at this point.
   styleNameEnv <- new.env(hash=TRUE, parent=emptyenv())
   assign('styleNameEnv', styleNameEnv, pos=.odfEnv)

   # Parse content.xml
   content.xml.doc <- parseXML("content.xml")
   top <- getTopNode(content.xml.doc)

   # Initialize the "Style Name Environment"
   initStyleNames(top, styleNameEnv)

   # Initialize seqInfo
   seqInfo <- getSeqInfo(top)
   assign('seqInfo', seqInfo, pos=.odfEnv)

   announce(verbose, "\n  Pre-processing the contents\n")
   # pre-process content.xml in preparation for sweaving
   rnwFileName <- "content.Rnw"
   preproc(top, rnwFileName)

   # Create the "New Style Environment" which will be used to register
   # styles (in the form of XMLNode objects) that will later be added
   # to the document during post processing.
   # Note that this assignment may overwrite an environment created on
   # an earlier call to odfWeave.  We need a new, clean one at this point.
   assign('newStyleEnv', new.env(hash=TRUE, parent=emptyenv()), pos=.odfEnv)

   # Sweave results to new xml file
   announce(verbose, "  Sweaving ", rnwFileName, "\n\n")

   Sweave(file=rnwFileName, output="content_1.xml",
      quiet=!control$verbose, driver=RweaveOdf(), control=control,
      encoding="UTF-8")

   # reset the figure captions
   
   .odfEnv$fig.caption <- NULL

   if (!control$debug)
   {
      # remove the original content.xml
      announce(verbose, "\n  Removing content.xml\n")
      file.remove("content.xml")
      if (file.exists("content.xml")) stop("Error removing content.xml file")
   } else {
      announce(verbose, "  Renaming original content.xml to content_orig.xml\n")
      file.rename("content.xml", "content_orig.xml")
      if (file.exists("content.xml")) stop("Error renaming content.xml file")
   }

   announce(verbose, "\n  Post-processing the contents\n")
   # post-process the output from Sweave
   content_1.xml.doc <- parseXML("content_1.xml")
   top <- getTopNode(content_1.xml.doc)
   postproc(top, "content.xml")

   if (!control$debug)
   {
      # remove the input to Sweave
      announce(verbose, "  Removing", rnwFileName, "\n")
      file.remove(rnwFileName)
      if (file.exists(rnwFileName)) stop("Error removing xml file")
   } else {
      announce(verbose, "  Not removing", rnwFileName, "\n")
   }

   # process styles.xml
   styles.xml.doc <- parseXML("styles.xml")
   stylestop <- getTopNode(styles.xml.doc)
   procstyles(stylestop, "styles_2.xml")

   if (!control$debug)
   {
      # remove original styles.xml file
      announce(verbose, "  Removing styles.xml\n")
      file.remove("styles.xml")
      if (file.exists("styles.xml")) stop("Error removing styles.xml file")
   } else {
      announce(verbose, "  Renaming original styles.xml to styles_orig.xml\n")
      file.rename("styles.xml", "styles_orig.xml")
      if (file.exists("styles.xml")) stop("Error renaming styles.xml file")
   }

   # rename post-processed file to styles.xml ready for zipping
   announce(verbose, "  Renaming styles_2.xml to styles.xml\n")
   file.rename("styles_2.xml", "styles.xml")
   if (!file.exists("styles.xml")) stop("Error renaming styles xml file")

   # verbose listing of picture files
   for (picFile in .odfEnv$picVector)
   {
      announce(verbose, paste("  Collected picture file ", picFile, "\n", sep=""))
   }

   # process META-INF/manifest.xml
   manifest.xml.doc <- parseXML("META-INF/manifest.xml")
   manifesttop <- getTopNode(manifest.xml.doc)
   procmanifest(manifesttop, "META-INF/manifest_2.xml")

   if (!control$debug)
   {
      # remove original manifest.xml file
      announce(verbose, "  Removing manifest.xml\n")
      file.remove("META-INF/manifest.xml")
      if (file.exists("META-INF/manifest.xml")) stop("Error removing manifest.xml file")
   } else {
      announce(verbose, "  Renaming original manifest.xml to manifest_orig.xml\n")
      file.rename("META-INF/manifest.xml", "META-INF/manifest_orig.xml")
      if (file.exists("META-INF/manifest.xml")) stop("Error renaming manifest.xml file")
   }

   # rename post-processed file to manifest.xml ready for zipping
   announce(verbose, "  Renaming manifest_2.xml to manifest.xml\n")
   file.rename("META-INF/manifest_2.xml", "META-INF/manifest.xml")
   if (!file.exists("META-INF/manifest.xml")) stop("Error renaming manifest xml file")

   # do final cleanup if debugging isn't enabled
   if (!control$debug)
   {
      announce(verbose, "  Removing extra files\n")
      if(file.exists("content_1.xml")) try(file.remove("content_1.xml"), silent = TRUE)
   } else {
      announce(verbose, "  Not removing extra files\n")
   }

   # zip up the new ODT file
   announce(verbose, "\n\  Packaging file using", zipCmd[1], "\n")
   if(.Platform$OS.type == "windows") {
      if(system(zipCmd[1], invisible=TRUE) != 0)  stop("Error zipping file")
   } else {
      if(system(zipCmd[1]) != 0) stop("Error zipping file")
   }

   # copy new ODT file to destination
   if(!file.exists(workingCopy))  stop(paste(workingCopy, "does not exist"))
   announce(verbose, "  Copying ", workingCopy, "\n")
   if(!file.copy(workingCopy, dest, overwrite=TRUE))  stop("Error copying odt file")

   # set the current working directory to the original location
   announce(verbose, "  Resetting wd\n")
   setwd(currentLoc)

   assign(
      "picPath",
      NA,
      pos = .odfEnv)
   assign(
      "picVector",
      NA,
      pos = .odfEnv)

   # Delete tmp dir that we were working in if we created it
   if(control$cleanup && created.workDir)
   {
      announce(verbose, "  Removing ", workDir, "\n")
      unlink(workDir, recursive=TRUE)
      # make repeated attempts to remove the directory
      if (file.exists(workDir)) unlink(shQuote(workDir), recursive=TRUE)
      if (file.exists(workDir)) unlink(shellQuote(workDir), recursive=TRUE)
      if (file.exists(workDir)) stop("Error removing work dir")
   } else {
      announce(verbose, " Not removing ", workDir, "\n")
   }
   announce(verbose, "\n  Done\n")   
   invisible(NULL)
}

# This works like normalizePath except that it handles the case
# where the path doesn't exist.  The value of dirname(path) must
# exist, however.
"canonicalFilePath" <- function(path)
{
   d <- dirname(path)
   if (! file.exists(d))
      stop('directory does not exist: ', d)
   if (! file_test("-d", d))
      stop(d, ' is not a directory')
   file.path(normalizePath(d), basename(path))
}

# This is a wrapper for normalizePath, but issues specific
# errors rather than a warning if the path doesn't exist.
"canonicalDirPath" <- function(path)
{
   if (! file.exists(path))
      stop('directory does not exist: ', path)
   if (! file_test("-d", path))
      stop(path, ' is not a directory')
   normalizePath(path)
}

"announce" <- function (verbose = TRUE, ...) 
{
   if (verbose) cat(...)
   flush.console()
   invisible()
} 

"attR<-" <- function (x, name, value) {
   #apply an attribute lists in a list
   #value always gets unlisted before assignment
   if( (length(value) != length(x)) || (length(value) == 1))
      lapply(x, function(x) {attr(x, name) <- value[[1]]; x})
   else
      lapply(seq(x), function(y) { attr(x[[y]], name) <- value[[y]]; x[[y]]})
}

"attR" <- function (x, name) {
   #return the values of attribute "name" items in list
   lapply(x, function(x) attr(x, name))
}

"getByAttR" <- function (x, name, value) {
 #return items in list x with attribute called "name"
 x[unlist(lapply(x, function(y) attr(y, name) == value))]
}

"debugWrite" <- function(filename, content){
   #used at breakpoints during debugging
   outFile <- file(filename, "wb")
   writeBin(content, outFile)
   close(outFile)
}

"shellQuote" <- function(x) {
   if (.Platform$OS.type == "windows"){
      return(paste('"', x, '"', sep=""))
   } else {
      return(gsub(" ", "\\\\ ", x))
   }
}

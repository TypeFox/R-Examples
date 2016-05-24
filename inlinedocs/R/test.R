test.file <- function
### Check an R code file with inlinedocs to see if the
### extract.docs.file parser accurately extracts all the code inside!
### The code file should contain a variable .result which is the
### documentation list that you should get when you apply
### extract.docs.file to the file. We check for identity of elements
### of elements of the list, so the order of elements should not
### matter, and thus this should be a good robust unit test.
(f,
### File name of R code file with inlinedocs to parse and check.
 verbose=TRUE
### Show output?
 ){
  ##seealso<< \code{\link{save.test.result}}
  e <- new.env()
  suppressWarnings(sys.source(f,e))
  ## these are the items to check for, in no particular order
  .result <- e$.result
  parsers <- e$.parsers
  result <- extract.docs.file(f,parsers)
  for(FUN in names(.result)){
    if(verbose)cat(FUN,"")
    fun <- result[[FUN]]
    .fun <- .result[[FUN]]
    ## first check to make sure all the stored items are there
    for(N in names(.fun)){
      .res <- .fun[[N]]
      res <- fun[[N]]
      if(is.null(res) || is.na(res) || is.na(.res) || .res!=res){
        cat("\n-----\n",res,"\n-----\nin ",FUN,
            "$",N,", expected:\n-----\n",.res,"\n-----\n")
        stop("docs mismatch in ",f)
      }
    }
    ## now check and see if there are no additional items!
    additional <- !names(fun)%in%names(.fun)
    show <- fun[additional] ##ignore NULL extracted items
    show <- show[!sapply(show,is.null)]
    if(length(show)){
      cat("\n")
      print(show)
      stop("extracted some unexpected docs!")
    }
  }
  ## make sure there are no unexpected outer lists
  not.expected <- names(result)[!names(result)%in%names(.result)]
  if(length(not.expected)){
    print(not.expected)
    stop("extracted some unexpected documentation objects!")
  }
  ## finally make a package using this file and see if it passes
  ## without warnings TDH 27 May 2011 added !interactive() since
  ## recursive calling R CMD check seems to systematically
  ## fail... ERROR: startup.Rs not found. This file is usually copied
  ## to the check directory and read as a .Rprofile, as done in
  ## tools:::.runPackageTests ... is this a bug in R? Anyway for now
  ## let's just not run the R CMD check.
  if(!is.null(e$.dontcheck) || !interactive())return()
  make.package.and.check(f,parsers,verbose)
  if(verbose)cat("\n")
}

make.package.and.check <- function
### Assemble some R code into a package and process it using R CMD
### check, stopping with an error if the check resulted in any errors
### or warnings.
(f, ##<< R code file name from which we will make a package
 parsers=default.parsers,
### Parsers to use to make the package documentation.
 verbose=TRUE
### print the check command line?
 ){
  pkgname <- sub("[.][rR]$","",basename(f))
  pkgdir <- file.path(tempdir(),pkgname)
  if(file.exists(pkgdir))unlink(pkgdir,recursive=TRUE)
  rdir <- file.path(pkgdir,"R")
  if(verbose)print(rdir)
  dir.create(rdir,recursive=TRUE)
  sillydir <- system.file("silly",package="inlinedocs")
  tocopy <- file.path(sillydir,c("DESCRIPTION","NAMESPACE"))
  file.copy(tocopy,pkgdir)
  file.copy(f,rdir)
  package.skeleton.dx(pkgdir,parsers)
  cmd <- sprintf("%s CMD check --as-cran %s",
                 file.path(R.home("bin"), "R"),
                 pkgdir)
  if(verbose)cat(cmd,"\n")
  checkLines <- system(cmd,intern=TRUE)
  warnLines <- grep("(WARNING|ERROR)",checkLines,value=TRUE)
  if(length(warnLines)>0){
    print(warnLines)
    stop("ERROR/WARNING encountered in package check!")
  }
}

save.test.result <- function
### For unit tests, this is an easy way of getting a text
### representation of the list result of extract.docs.file.
(f
### R code file with inlinedocs to process with extract.docs.file.
 ){
  .result <- extract.docs.file(f)
  dump(".result",tmp <- tempfile(),control=NULL)
  lines <- readLines(tmp)
  cat(paste(lines,"\n"))
}

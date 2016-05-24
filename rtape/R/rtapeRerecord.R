#' Rerecord the tape.
#'
#' This function reads the objects from one tape, executes a callback function on them and updates them with/appends to the other tape the objects that the callback has returned. 
#'
#' @param FUN Callback function which transforms the objects.
#' @param fNamesIn Name of the tape file to read; if this argument is a vector of several names, function behaves as reading a single tape made of all those tapes joined in a given order. 
#' @param fNameOut Name of the tape to which store the output of filtering; if this file is one of the input files, this file is overwritten with the output; otherwise the output is appended to this tape. This must be a one-element vector.
#' @param moreArgs Additional arguments to \code{FUN}, given as a list.
#' @param skipNULLs If true, all the \code{NULL}s returned by \code{FUN} are not appended to the output. Useful to remove some objects from the tape; see \code{\link{rtapeFilter}} for convenience function to do just this task.
#' @param fileFormatOut File format; should be left default. See \code{\link{guessFileFormat}} and \code{\link{makeFileFormat}} for details.
#' @note Overwriting is NOT realised in place, rather by a creation of a temporary file and then using it to overwrite the filtered tape.
#' @author Miron B. Kursa \email{M.Kursa@@icm.edu.pl}
#' @examples
#' unlink(c('tmp.tape','tmp2.tape'))
#' #Record something
#' for(i in 1:10) rtapeAdd('tmp.tape',i)
#'
#' #Multiply each object by two
#' rtapeRerecord('*','tmp.tape','tmp2.tape',moreArgs=list(2))
#'
#' #Check it out
#' unlist(rtapeAsList('tmp.tape'))->A
#' B<-unlist(rtapeAsList('tmp2.tape'))
#' print(A);print(B)
#' stopifnot(all(A==B/2))
#'
#' #Now do the same in-place:
#' rtapeRerecord('*','tmp.tape',moreArgs=list(2))
#' unlist(rtapeAsList('tmp.tape'))->B2
#' stopifnot(all(A==B2/2))
#'  
#' #Time to clean up
#' unlink(c('tmp.tape','tmp2.tape'))

rtapeRerecord<-function(FUN,fNamesIn,fNameOut=fNamesIn,moreArgs=NULL,skipNULLs=FALSE,fileFormatOut=guessFileFormat(fNameOut)){
 stopifnot(length(fNameOut)==1)
 FUN<-match.fun(FUN)
 if(fNameOut%in%fNamesIn){
  fNameOut->fNameOverwrite
  fNameOut<-sprintf("%s.tmp_%d",fNameOut,Sys.getpid())
 }
 fileFormatOut(fNameOut,open="ab")->conOut
 rtape_apply(fNamesIn,function(x){
  what<-do.call(FUN,c(list(x),moreArgs))
  if(!(skipNULLs & is.null(what))) serialize(what,conOut,ascii=FALSE)
  NULL
 })
 close(conOut)
 if(exists('fNameOverwrite')){
  if(file.exists(fNameOverwrite)){
    file.remove(fNameOverwrite)
    file.rename(fNameOut,fNameOverwrite)
  }else stop("Something wrong with temporary file. Target tape was NOT altered.");
 }
 return(invisible(NULL))
}


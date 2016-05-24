#' Load the whole tape as a list.
#'
#' This function reads are the objects from the tape, in the order they were written on it, and returns them as a list.
#'
#' @param fNames Name of the tape file to read; if this argument is a vector of several names, function behaves as reading a single tape made of all those tapes joined in a given order. 
#' @return A list containing all the objects stored on the tape.
#' @author Miron B. Kursa \email{M.Kursa@@icm.edu.pl}
#' @examples
#' unlink('tmp.tape')
#' #Record something on the tape
#' rtapeAdd('tmp.tape',iris)
#' rtapeAdd('tmp.tape',sin(1:10))
#' #Read whole tape to the list, so we could examine it
#' rtapeAsList('tmp.tape')->stored
#' print(stored)
#' unlink('tmp.tape')

rtapeAsList<-function(fNames){
 rtapeLapply(fNames,function(x) x)
}


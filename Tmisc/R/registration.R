#' Get registration link
#' 
#' Returns a link to register for a course. This reads a text file that I keep 
#' on my website that reveals a URL used to register for a course. I change the 
#' URL on my website, and the result of running the function changes. I use this
#' so that I can give instructions needed to set up for a course, which includes
#' instructions to download this package and run this function. This is used to
#' ensure that participants were able to successfully download and install R and
#' the required packages before coming to the course.
#' 
#' @author Stephen Turner
#' @keywords NA
#'   
#' @return A registration link
#' 
#' @examples
#' registration()
#'   
#' @export
registration <- function() {
    # Where do you keep the registration link?
    myurl <- "http://stephenturner.us/f/Tmisc-registration-link.txt"
    reglink <- scan(myurl, what="char", sep="\n", n=1, quiet=TRUE)
    if (reglink=="") reglink <- NA
    message(paste0("Please visit the following link to register:\n", reglink))
}

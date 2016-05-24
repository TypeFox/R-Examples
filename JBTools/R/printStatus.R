printStatus <- function(
##title<< Print a status message
    string.in ##<< character string: status message to print
)
  ##description<< This function prints a given string as a status report to the screen together
  ## with the current time.
{
    cat(paste(Sys.time(),' : ',string.in,'\n',sep=''))
    ##value<<   Nothing is returned but a message is printed to the screen.
}

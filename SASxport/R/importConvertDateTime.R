##
## Code originally from Frank Harrell's 'Hmisc' library: 
##   http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/Hmisc
## Copied with permission on 2007-08-04
##

#' @importFrom chron times chron

importConvertDateTime <- 
  function(x, type=c('date','time','datetime'), input='sas', form)
{
  type <- match.arg(type)
  input <- match.arg(input)

  adjdays <- 3653   ## 1970-1-1 minus 1960-1-1

  switch(type,
         date = structure(x - adjdays, class='Date'), 
         time = times(x/86400, out.format=c(dates='day mon year', times='h:m:s')),
         datetime = chron((x - adjdays*86400)/86400,
                          out.format=c(dates='day mon year', times='h:m:s'))
         )
}

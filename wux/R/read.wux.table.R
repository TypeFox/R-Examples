
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-07 16:09:44 +0200 (Tue, 07 Apr 2015) $
# $Rev: 339 $
# ----------------------------------------------------------------

read.wux.table <- function(file, ...){
  ## Reads in the WUX data.frame csv file saved on harddisk.
  ## 
  ## Args:
  ##   file: Character filename
  ##   ...: further args used in read.table
  ## 
  ## Returns:
  ##   A wux.df object (basically a data.frame) .
  ##
  ## History:
  ##   2011-09-29 | original code

  data.frame.out <- read.table(file, sep = ";", header = TRUE, ...)
  class(data.frame.out) <- append("wux.df", class(data.frame.out))
  
  return(data.frame.out)
}


                            

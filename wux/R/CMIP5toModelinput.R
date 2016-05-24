CMIP5toModelinput <- function(filedir = NULL,
                              save.to = NULL,
                              verbose = FALSE) {
  ## Creates a modelinput list based on downloaded CMIP5 data by
  ## CMIP5fromESGF. This function executes a python program in the
  ## background.
  ## 
  ## Args:
  ##   filedir: Direcotry where the CMIP5 data are being stored.
  ##            It is crucial for this funcion that the data is
  ##            saved in the same directory strucure as created
  ##            by "CMIP5fromESGF".
  ##   save.to: character. Filename to safe the modelinupt.
  ##   verbose: boolean. For additional information printed on the screen.
  ##
  ## Returns:
  ##   Nothing.
  ##
  ## History:
  ##   2015-04-02 | original code (thm)
  ##
  
  if (is.null(filedir))
    stop("Please provide the directoy of your CMIP5 data storage.")
  if (is.null(save.to))
    stop("Please provide a file for data storage (save.to parameter).")
  
  python.filename <- system.file('exec', 'cmip5_to_wux_modeldict.py', package='wux')
  
  command <- paste("python", python.filename,
                   "-i", filedir,
                   "-o", save.to)
  ## execute python script
  system(command, ignore.stdout = !verbose)
}

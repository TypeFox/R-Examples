xport.namestr <- function(
                          # required arguments
                          var,       # value of variable
                          varName,   # name of variable

                          varNum,    # variable number (starting at 1)
                          varPos,    # record position of varible (starting at 0) 

                          # optional arguments
                          varLength, # variable length
                          varLabel="",  # variable label (uses 'label'
                                     # attribute if present, otherwise
                                     # defaults to R  variable name)

                          fName="",  
                          fLength=0, 
                          fDigits=0,
                          just=c("left","right"),
                          
                          iName="",
                          iLength=0,
                          iDigits=0
                          )
{
  if(is.factor(var))
    var <- as.character(var)
  
  isChar = is.character(var)

  if(missing(varLength))
    if(isChar)
      varLength <- max(nchar(var))
    else
      varLength <- 8

  if( missing(varLabel) || is.null(varLabel) )
    varLabel <- ""    

  just <- match.arg(just)
  if(just=="left")
    justVal <- 0
  else
    justVal <- 1
  
  
  ## force variable name into format permitted by SAS.  Starts with
  ## alpha, alpha, numbers, and underscore permitted.  R's
  ## make.names() function almost does what we want, but allows
  ## periods instead of underscores, so just use make.names() and
  ## replace periods with underscores, :-)
  varName <- gsub("\\.","_", make.names(varName))

  ## Note that the variable name field in the xport file only permits
  ## 8 characters, so names will be truncated.

  .C("fill_namestr",
     isChar = as.integer(isChar),              # Bool: Is this a character varible
     nlng   = as.integer(varLength),           # LENGTH OF VARIABLE IN OBSERVATION
     nvar0  = as.integer(varNum),              # VARNUM
     nname  = toupper(as.character(varName)),  # NAME OF VARIABLE
     nlabel = as.character(varLabel),          # LABEL OF VARIABLE
		   
     nform  = toupper(as.character(fName)),    # NAME OF FORMAT
     nfl    = as.integer(fLength),             # FORMAT FIELD LENGTH OR 0
     nfd    = as.integer(fDigits),             # FORMAT NUMBER OF DECIMALS
     nfj    = as.integer(justVal),             # 0=LEFT JUSTIFICATION, 1=RIGHT JUST
		  
     niform = toupper(as.character(iName)),    # NAME OF INPUT FORMAT
     nifl   = as.integer(iLength),             # INFORMAT LENGTH ATTRIBUTE
     nifd   = as.integer(iDigits),             # INFORMAT NUMBER OF DECIMALS

     npos   = as.integer(varPos),              # POSITION OF VALUE IN OBSERVATION  
     PACKAGE="SASxport"
     )

    .Call("getRawBuffer", PACKAGE="SASxport")
}


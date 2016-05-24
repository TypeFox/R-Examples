"vfonal" <-
function (MICIKE=NULL) {

  #############################################################
  # 
  # TITLE:  	vfonal()
  # AUTHOR: 	SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   	26 NOV, 2003
  # CALLS:  
  # CALLED BY:
  # NEEDS:  	
  # NOTES:  	MUCH LIKE MY MAKEPATH(), BUT NOT AS FLEXIBLE
  #			ALLOWS THE CONSTRUCTION OF CUSTON DECOMPOSITION
  #			PATHS 
  #
  #############################################################


  # RESET DECOMPOSITION PATH TITLE
  CIM <- ""

  # IF DECOMPOSITION PATH IS NOT GIVEN, ENTER INTERACTIVE MODE
  if(is.null(MICIKE)){
    BE <- names(.QND)
    H <- length(BE)
    VFONAL <- matrix(-1,nrow=2,ncol=H,dimnames=list(NULL,BE))
    SZOVEG <- "\n\nHelp for DecoPath matrix\nwrite 0 for NullHipo (typically either Z or X)"
    SZOVEG <- paste(collapse="",SZOVEG,"\nwrite 1 for AltHipo (typically either Y or YandZ)")
    SZOVEG <- paste(collapse="",SZOVEG,"\n\nDecoPath title: ")
    CIM <- readline(SZOVEG)
    assign(".VFONAL",VFONAL,pos=1)
    fix(.VFONAL)
  } 
  # OTHERWISE USE THE SUPPLIED DECOMPOSITION PATH MATRIX
  else assign(".VFONAL",MICIKE,pos=1)

  vfon(BE=.VFONAL)
  
  if(nchar(CIM)!=0) attr(.VFONAL,"cim") <- CIM
  assign(".VFONAL",.VFONAL,pos=1)
}


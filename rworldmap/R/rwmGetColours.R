#' to choose map colours for classified data
#' 
#' Returns a vector of colours based upon the palette specified and number of
#' colours specified.  If colourPalette specifies a number of colours and this
#' is different from numColours, numColours takes precedence and colours are
#' interpolated to make the number fit.
#' 
#' @param colourPalette string describing the colour palette to use, choice of:
#' \enumerate{ \item "palette" for the current palette \item a vector of valid
#' colours, e.g. =c('red','white','blue') or output from RColourBrewer \item one
#' of "heat", "diverging", "white2Black", "black2White", "topo", "rainbow",
#' "terrain", "negpos8", "negpos9" }
#' @param numColours the number of colour categories desired
#' @return A vector specifiying a number of colours.
#' @export
rwmGetColours <- function(colourPalette, numColours)
{

#if invalid option is chosen it's set to palette at end
paletteList <- c("white2Black","black2White","palette","heat","topo","terrain","rainbow","negpos8","negpos9")


#to allow user specified colours
#this could test whether colours are valid
if( length(colourPalette) > 1 ) #only if the colourPalette passed is a vector of length > 1
	{
		#coloursToUse <- colourPalette
		
	  #If the palette has a different number of colours to the number of breaks, interpolation is used.
    #if(length(colourPalette)==numColours) {coloursToUse<-colourPalette}
    if(length(colourPalette)==numColours) {coloursToUse<-as.character(colourPalette)}
    else 
        {
         warning(length(colourPalette)," colours specified and ", numColours," required, using interpolation to calculate colours")
         coloursToUse<-colorRampPalette(colourPalette,alpha=TRUE)(numColours)
         #coloursToUse<-colorRampPalette(colourPalette,space = "Lab")(numColours) #tried this to preserve transparency, didn't work
        }					
	} else
if(colourPalette=="negpos9")
	{
	  #later add that these can cope with diff num cats too
		coloursToUse <- c('darkblue','blue','lightblue','lightgreen','grey','yellow','orange','red','darkred')#'purple')
	} else		
if(colourPalette=="negpos8")
	{
		coloursToUse <- c('darkblue','blue','lightblue','lightgreen','yellow','orange','red','darkred')#'purple')
	} else
	
if(colourPalette=="diverging")
	{
   coloursAbove <- colorRampPalette(c('yellow','red2'))(floor(numColours/2))
   coloursBelow <- colorRampPalette(c('darkblue','lightgreen'))(floor(numColours/2))   			

	if ( numColours%%2 == 0 ) #even
	   {
	    coloursToUse <- c(coloursBelow,coloursAbove)				   
	   } else #i.e. odd
	   {
	    #adding colour for central break
	    coloursToUse <- c(coloursBelow,'white',coloursAbove)	    
	   }		
   
	} else
  	
if(colourPalette=="white2Black")
	{
		white2Black <- colorRampPalette(c(grey(0.8),grey(0.2)))
		coloursToUse <- white2Black(numColours)
	} else

if(colourPalette=="black2White")
	{
		black2White <- colorRampPalette(c(grey(0.2),grey(0.8)))
		coloursToUse <- black2White(numColours)
	} else

if(colourPalette=="palette")
	{	
	  #If the palette has a different number of colours to the number of breaks, interpolation is used.
	  #This makes custom palettes easier to use. You can change numCats in a mapping call with out changing the palette.
	  #It also makes it easier to use custom palettes with pretty or quantiles, where the number of categories is always numCats.
    if(length(palette())==numColours) {coloursToUse<-palette()}
    else 
        {
         warning(length(palette())," colours specified and ", numColours," required, using interpolation to calculate colours")
         coloursToUse<-colorRampPalette(palette())(numColours)
        }
  } else
	
if(colourPalette=="heat") { coloursToUse <- rev(heat.colors(numColours))	
	} else
if(colourPalette=="topo") { coloursToUse <- rev(topo.colors(numColours))	
	} else
if(colourPalette=="terrain") { coloursToUse <- rev(terrain.colors(numColours))	
	} else
if(colourPalette=="rainbow") { coloursToUse <- rev(rainbow(numColours))	
	} else	
	{	
    warning("colourPalette should be set to either a vector of colours or one of :",paste(paletteList,""),"\nsetting to heat colours as default") 
		coloursToUse <- rev(heat.colors(numColours))
	} 	
	

return(coloursToUse)

} #end of rwmGetColours


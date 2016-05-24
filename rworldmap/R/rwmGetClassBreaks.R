#' Internal function to set the numeric values for the breaks between data
#' categories
#' 
#' Sets the values that determine how a vector of continuous data is classified
#' into categories. Called by mapCountryData() and mapGriddedData()
#' 
#' 
#' @param dataColumn the data vector to be classified, must be numeric
#' @param catMethod the method to use to classify the data into categories,
#' choice of "pretty", "fixedWidth", "diverging",
#' "logFixedWidth","quantiles","categorical" or a numeric vector defining
#' breaks
#' @param numCats number of categories to put the data in, may be overidden if
#' not possible under some classification methods
#' @param verbose whether to print information messages to console TRUE/FALSE
#' @param midpoint the midpoint to use if catMethod='diverging', default=0
#' @return A vector specifying the numeric breaks between data categories.
#' @author andy south and matthew staines
#' @seealso The classInt package
#' @keywords dplot
#' @export rwmGetClassBreaks
`rwmGetClassBreaks` <-
function(dataColumn, catMethod, numCats, verbose=TRUE, midpoint=0)
{

#! I should offer option to set min & max (& centre for diverging)

#browser()

functionName <- as.character(sys.call()[[1]])

catMethodList <- c("fixedWidth","diverging","quantiles","pretty","logFixedWidth","categorical")

if( ! catMethod %in% catMethodList)
  {
   warning("classification method should be set to one of :",paste(catMethodList,""),"\nsetting to fixedWidth as default\n")
   catMethod="fixedWidth"
  }

if(catMethod=="fixedWidth")
		{
			#Categorising the data, fixed width intervals
			minVal <- min(dataColumn,na.rm=TRUE)
			maxVal <- max(dataColumn,na.rm=TRUE) 
			
			cutVector <- minVal + ( ((0:numCats)/numCats) * (maxVal-minVal) )				
		} else
		
if(catMethod=="diverging")
		{
			#Categorising the data, fixed width intervals
			minVal <- min(dataColumn,na.rm=TRUE)
			maxVal <- max(dataColumn,na.rm=TRUE) 
			
			#?set the break interval across the whole range
			#?or have options for diff scales above & below
			
			#will probably eventually need to have extra arguments
			above <- abs(maxVal-midpoint)
			below <- abs(midpoint-minVal)
			
			#num categories above or below
			#if numCats is odd this will include .5
      sideCats <- numCats/2
			
			interval <- max(c(above,below)) / sideCats
			
			if ( numCats%%2 == 0 ) #even
			   {
			    fromAbove <- midpoint + interval
			    fromBelow <- midpoint - interval
			   } else #i.e. odd
			   {
			    fromAbove <- midpoint + interval/2
			    fromBelow <- midpoint	- interval/2		    
			   }
			
			cutsAbove <- seq(from=fromAbove, to=midpoint+(sideCats*interval), by=interval)
			cutsBelow <- seq(from=fromBelow, to=midpoint-(sideCats*interval), by=-interval)
			
			if ( numCats%%2 == 0 ) #even
			   {
			    #adding in the midpoint
			    cutVector <- c(rev(cutsBelow),midpoint,cutsAbove)				   
			   } else #i.e. odd
			   {
			    cutVector <- c(rev(cutsBelow),cutsAbove)	    
			   }			
						
		} else	
    	
if(catMethod=="quantiles")
		{
			#Categorising the data, using Quantiles.

   #Using Quantiles will crash if the data contains too many repeats and numCats is high.
   #The break points must be unique. The algorithm below will use numCats if it can.
   #If numCats does not produce unique breakpoints,
   #it will keep reducing the number of quantiles it will use, till unique break points are found.
   #It will also warn if the number of quantiles used was less than asked for.

    testNumCats<-numCats              #The next number of quantiles to try. starts at numCats, and decreases till unique breeakpoints are found.
    uniqueBreaksFlag<-FALSE           #Flags if unique breaks have been found. When TRUE, the while loop stops, and the current value of testNumCats is used to produce quantiles.
    while(uniqueBreaksFlag==FALSE && testNumCats > 0 )
              {
              testQuantiles<-quantile(dataColumn,probs = seq(0, 1, 1/testNumCats),na.rm=TRUE)

              if(length(testQuantiles)==length(unique(testQuantiles)))     #Are the breaks unique?
                    {
                    uniqueBreaksFlag<-TRUE   #Stop looping
                    }else{
                    testNumCats<-testNumCats-1         #Carry on looping, trying one fewer quantile.
                    }
              }
    if(testNumCats!=numCats && verbose )message(paste("You asked for",numCats,"quantiles, only",testNumCats,"could be created in quantiles classification"))  #Warning if the number of quantiles was reduced.

    cutVector <-  quantile(dataColumn, probs=seq(0,1, 1/testNumCats), na.rm=TRUE)

    } else
    
if(catMethod=="pretty")
		{
			#Compute a sequence of about n+1 equally spaced values
      #which cover the range of the values in x.
      #The values are chosen so that they are 1, 2 or 5 times a power of 10.

			cutVector <- pretty(dataColumn, n=numCats)

			#Pretty will choose a number of categories similar to the number of categories asked for.
			#The following code warns when pretty has used a different number of breaks to that which was asked for.

      actualNumberOfBreaks<-length(cutVector)-1
      if(actualNumberOfBreaks!=numCats && verbose ) message(paste("You asked for",numCats,"categories,",actualNumberOfBreaks, "were used due to pretty() classification"))

		} else

# if min = 0 adds 0.01 to avoid problems with zeroes
if ( catMethod=="logFixedWidth") 
    {
         
      # to do for Logs will want to Log the data calc the CutVector then antiLog
      # to get a cutVector that can be directly applied to the data
      
      if (min( dataColumn, na.rm=TRUE ) < 0 ) 
         {stop("negative values in your data cannot be classified using catMethod=logFixedWidth")
          return(FALSE) 
         } else if (min( dataColumn, na.rm=TRUE ) == 0 )
         {
          if (verbose) message("zero values are replaced with NA as they can't be logged in catMethod=logFixedWidth")

          dataColumn[which(dataColumn==0)] <- NA
          #dataColumnLogged <- log(addTo0ForLog+dataColumn)         
          dataColumnLogged <- log(dataColumn) 
         } else
         {
          dataColumnLogged <- log(dataColumn) 
         }
              		
      minVal <- min(dataColumnLogged,na.rm=TRUE)
			maxVal <- max(dataColumnLogged,na.rm=TRUE) 
			maxValNotLogged <- max(dataColumn,na.rm=TRUE)
      			
			#there was a rounding problem with this, that meant that highest value could get excluded
			cutVector <- minVal + ( ((0:numCats)/numCats) * (maxVal-minVal) )	
       			
 			#antilog
      #cutVector <- exp(cutVector) -  exp(log(addTo0ForLog))
      cutVector <- exp(cutVector)

      #to correct potential rounding problem, make sure upper val is equal to max value
      cutVector[length(cutVector)] <- maxValNotLogged
      
    }

if(length(catMethod)==1 && catMethod=="categorical")
		{
			stop(functionName," shouldn't be called when catMethod == 'categorical'")
      return(0)			
		} else 
    {		
     return(cutVector)
		}
} #end of rwmGetClassBreaks


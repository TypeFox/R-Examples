#' Creates a plot device set up for maps
#' 
#' Creates a plot device suited for rworldmap plotting functions.
#' 
#' 
#' @param device Character string which controls the type of plot default.  The
#' default uses your standard plot device. Giving the name of a plotting device
#' function will use that instead. e.g. "pdf", "png", etc.
#' @param rows The number of rows. Default 1
#' @param columns The number of columns. Default 1
#' @param plotOrder Option of 'rows' or 'columns'. For multiple plots whether
#' to plot in row or column order. However, note that addMapLegend can have the
#' effect of reverting order to rows.
#' @param width The width of a single plot. This includes the margins. If you
#' do not specify both width and height, suitable values will be calculated
#' @param height The height of a single plot. This includes the margins. If you
#' do not specify both width and height, suitable values will be calculated
#' @param titleSpace The height in inches of the gap at the plot.
#' @param mai The margin sizes in inches. If titleSpace is given this overrides
#' mai[3].
#' @param mgp As per par(mgp) in the graphics package
#' @param xaxs As per par(xaxs) in the graphics package
#' @param yaxs As per par(yaxs) in the graphics package
#' @param \dots Further arguments to the device function
#' @return Used for the side effect of creating a plot device, and setting
#' graphical parameters for the device.
#' @seealso mapCountryData,mapGridAscii
#' @keywords device
#' @examples
#' 
#' \dontrun{
#' #Basic Usage
#' mapDevice()
#' mapCountryData()
#' 
#' #2 by 2 plot
#' mapDevice(rows=2,columns=2)
#' columns<-c("BIODIVERSITY","EPI","ENVHEALTH","Population2005")
#' for(i in columns){
#'  mapCountryData(nameColumnToPlot=i)
#' }
#' #Creating a pdf that is 5 inches wide
#' mapDevice(device="pdf",width=5,file=tempfile())
#' mapCountryData()
#' dev.off()
#' 
#' }
#' 
#' @export mapDevice
`mapDevice`<-
function(device="dev.new"
                   ,rows=1
                   ,columns=1
                   ,plotOrder="rows"
                   ,width=NULL
                   ,height=NULL
                   ,titleSpace=NULL
                   ,mai=c(0,0,0.2,0)
                   ,mgp=c(0,0,0)
                   ,xaxs="i"
                   ,yaxs="i"
                   ,...)
{
  #aspectRatio = width/height. Unless you specify both width and height,
  #mapDevice will choose an aspect ratio for you, depending on the projection.
  #aspectRatio<-switch(projection,EqualArea=1.7,none=2)
  aspectRatio<-2

  #The margin with the title is only margin you would reguarly adjust,
  #so a shortcut is provided
  if(!is.null(titleSpace)){
   mai[3]<-titleSpace
  }
  
  #If at least one of height and width is NULL, they are set to sensible values
  if(is.null(height) && is.null(width)){
   if (device=='png') width=720 #default width in pixels
   else width=9 #this width works well for powerpoint for a single plot & gets reduced for word    
  }
  if(is.null(width)){
   width=(height-mai[1]-mai[3])*aspectRatio+mai[2]+mai[4]
  }
  if(is.null(height)){
   height=(width-mai[2]-mai[4])/aspectRatio+mai[3]+mai[1]
  }

  #Increase the device size for multiple plots.
  width=width*columns
  height=height*rows
  
  #Create a device  and set parameters.
  do.call(device,c(list(width=width,height=height),list(...)))
  par(mai=mai,mgp=mgp,xaxs=xaxs,yaxs=yaxs)
  
  #2/10/09 allow setting of whether, mfrow or mfcol
  if (plotOrder=='rows') 
     {par(mfrow=c(rows,columns))
     } else 
     {par(mfcol=c(rows,columns))
      print('####columns#####')
     }   
}




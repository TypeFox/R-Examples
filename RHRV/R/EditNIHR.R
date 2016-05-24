EditNIHR <-
function(HRVData,scale=1.0, verbose=NULL) {
#---------------------------------------
# Edits beats interactively
#	Requires tcltk and tkrplot libraries
#---------------------------------------
#	Plots Non-interpolated instantaneous heart rate for manual removing of outliers
#	scale -> allow scaling for small screens


	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Manually editing non-interpolated instantaneous heart rate **\n");
	}
	
	if (is.null(HRVData$Beat$Time)) { 
		cat("   --- ERROR: Beats not present!! ---\n")
		return(HRVData)
	}
	
	if (is.null(HRVData$Beat$niHR)) { 
		cat("   --- ERROR: Non-interpolated heart rate not present!! ---\n")
		return(HRVData)
	}
	
	editFunction <- function(HRVData) {
		
		HRVDataOld <- HRVData
	
		Myhscale <- 2*scale    # Horizontal scaling
		Myvscale <- 1.5*scale    # Vertical scaling
	
		plt <- c()
		usr <- c()
		coords <- c()
		pointsInArea <- c()
		numPointsInArea <- 0
		numCoords <- 0
		numRemovedPoints <- 0
		vectorx <- c()
		vectory <- c()
	

		plotFunction <- function()
		{
			vectorx <<- HRVData$Beat$Time
			vectory <<- HRVData$Beat$niHR
			plot(vectorx,vectory,type="l",xlab="time (sec.)",ylab="HR (beats/min.)",ylim=c(min(vectory),max(vectory)*1.1))
			title(main="Non-interpolated instantaneous heart rate")
		
			if (numCoords==1) {
				points(coords[1],coords[2],pch="+",col="red")
			}
			if (numCoords==2) {
				rect(min(coords[1],coords[3]),min(coords[2],coords[4]),
					max(coords[1],coords[3]),max(coords[2],coords[4]),border="red")
				areaString=paste("No. of selected points: ",numPointsInArea)
				text((coords[1]+coords[3])/2,max(coords[2],coords[4]),areaString,pos=3,col="red")
				points(vectorx[pointsInArea==TRUE],vectory[pointsInArea==TRUE],pch=20,col="red")
			
			}
		
			usr<<-par('usr')
			plt<<-par('plt')		
		}
	
		tt <- tktoplevel()
		tkwm.deiconify(tt)
		tkgrab.set(tt)
		tkfocus(tt)
		tkwm.title(tt,"Outliers removal")
		img <- tkrplot(tt,fun=plotFunction,hscale=Myhscale,vscale=Myvscale)
	
	
		Remove <- function()
		{
			numCoords <<-0
			coords <<- c()
			HRVData$Beat <<- subset(HRVData$Beat, pointsInArea==FALSE)
			numRemovedPoints <<- numRemovedPoints + numPointsInArea
			if (HRVData$Verbose) {
				cat("   Removing ",numPointsInArea," points (",numRemovedPoints," so far)\n",sep="")
			}			
			pointsInArea <<- c()
			numPointsInArea <<- 0
			tkrreplot(img)
		}
	
		Clear <- function()
		{
			numCoords <<- 0
			coords <<- c()
			pointsInArea <<- c()
			numPointsInArea <<- 0
			if (HRVData$Verbose) {
				cat("   Clearing point selection\n")
			}
			tkrreplot(img)
		}
	
		Quit <- function()
		{
			if (HRVData$Verbose) {
				cat("   Manual edition ended... quitting\n")
			}
			if (numRemovedPoints > 0) {
				msg <- paste(numRemovedPoints,"outliers to be removed\nProceed?")
				mbval<- tkmessageBox(title="Confirmation",message=msg,
					type="yesnocancel",icon="question")

				if (tclvalue(mbval)=="no") {
					tkgrab.release(tt)
					tkdestroy(tt)
					HRVData <<- HRVDataOld
				}
				if (tclvalue(mbval)=="yes") {
					tkgrab.release(tt)
					tkdestroy(tt)
				}		
			} else {
				tkdestroy(tt)
			}
		
		
		}
	
		OnLeftClick <- function(x,y)
		{
		
			xClick <- as.numeric(x)
			yClick <- as.numeric(y)
			width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
			height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
			#cat("Width:",width,"\n")
			#cat("Height:",height,"\n")
		
			xMin <- plt[1] * width
			xMax <- plt[2] * width
			yMin <- plt[3] * height
			yMax <- plt[4] * height
			#cat("xMin:",xMin,"\n")
			#cat("xMax:",xMax,"\n")
			#cat("yMin:",yMin,"\n")
			#cat("yMax:",yMax,"\n")		

			rangeX <- usr[2] - usr[1]
			rangeY <- usr[4] - usr[3]
			#cat("X range:",rangeX,"\n")
			#cat("Y range:",rangeY,"\n")
		
			xCoord <- usr[1]+(xClick-xMin)*rangeX/(xMax-xMin)
		 	yCoord <- usr[3]+((height - yClick)-yMin)*rangeY/(yMax-yMin)
		
			if (HRVData$Verbose) {
				cat("   Point clicked: (",xCoord,",",yCoord,")\n",sep="")
			}	
		
			if ((xClick>xMin*0.95)&&(xClick<xMax*1.1)&&(yClick>yMin*0.95)&&(yClick<yMax*1.1)) {
				#cat("Punto en gráfica\n")
				if (numCoords==0) {
					numCoords <<- 1
					#cat("Primer punto\n")
					coords <<- c(xCoord,yCoord)
				} else if (numCoords==1) {
					coords <<- c(coords,xCoord,yCoord)
					numCoords <<- 2
					pointsInArea <<- (
						(vectorx>min(coords[1],coords[3])) &
						(vectorx<max(coords[1],coords[3])) &
						(vectory>min(coords[2],coords[4])) &
						(vectory<max(coords[2],coords[4]))
					)
					numPointsInArea <<- length(pointsInArea[pointsInArea==TRUE])
					if (HRVData$Verbose) {
						cat("   ",numPointsInArea," points found in area\n",sep="")
					}
				}
				tkrreplot(img)
			}
		
		  
		}
	
		buttonremove <- tkbutton(tt,text="Remove outliers",command=Remove)
		buttonclear <- tkbutton(tt,text="Clear",command=Clear)	
		buttonremove2 <- tkbutton(tt,text="End",command=Quit)
	
	
		tkgrid(img,columnspan=3)
		tkgrid(buttonremove,buttonclear,buttonremove2)
		tkbind(img, "<Button-1>",OnLeftClick)
		tkwait.window(tt)
	
		return(HRVData)
	
	}
	
	HRVDataNew = editFunction(HRVData)
	
	return(HRVDataNew)

}


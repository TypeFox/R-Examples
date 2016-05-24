	
############################################################################################
# load BIFIEdata objects when objects are saved as full BIFIEdata
#    objects in one file
load.BIFIEdata <- function(filename, dir =getwd() ){
    d1 <- base::load( file=file.path(dir,filename) )
	objname <- "bdat_temp"
	cdata <- NULL
	miceadds::Reval( paste0("cdata <- " , d1  ,"$cdata" )	, print.string=FALSE  )
	if (cdata){		# if cdata=TRUE
	    l1 <- paste0( d1 , "$wgtrep <- matrix(" , 
					   d1,"$wgtreplist$unique_wgt[ " , d1 , "$wgtreplist$indexwgtrep ] ,
						     nrow= " , d1 , "$N , ncol= " , d1 , "$RR ) ")					
		miceadds::Reval( l1 , print.string=FALSE )						
		miceadds::Reval( paste0( d1 , "$wgtreplist <- NULL" ) , print.string=FALSE)
			}
	# save object in global environment
	# changed this setting: ARb 2014-07-29
#	eval(parse(text = paste(objname, "<<- ", d1)))	
	eval(parse(text = paste(objname, "<- ", d1)))	
    eval( parse(text= paste0( "return( " , objname , ")" ) ) )
    }
############################################################################################	

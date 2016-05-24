

#########################################################################
# miceadds::load.data: load conveniently R objects of different data formats
load.data <- function( filename , type="Rdata" , path=getwd() , 
				spss.default=TRUE , ...){
	#*** the resulting object is dat4!	
	dir <- path
	file <- filename
	
	i1 <- grep.vec( c("Rdata" , "csv" , "csv2" , "table" , "sav" ) , file ,
				"OR" )$x
	if ( length(i1) == 0 ){							
		files <- list.files( dir , filename )	
					} else {
		files <- file
					}
	type1 <- type
	if ( type=="table" ){
		files <- grep.vec( c("dat","txt") , files , "OR" )$x
		type1 <- "dat"
						}		

	files <- grep( gsub("csv2","csv" , type1) , files , value=TRUE)
	file <- max(files)
	cat( paste0( "*** Load " , file , "\n"))

    #*** Rdata objects	
	if (type == "Rdata" ){
		dat4 <- load.Rdata2( filename=file , path=dir )
				}
    #*** csv2 objects
	if (type == "csv2" ){
		dat4 <- utils::read.csv2( file.path(dir,file) , ... )
				}
    #*** csv objects
	if (type == "csv" ){
		dat4 <- utils::read.csv( file.path(dir,file) , ... )
				}
    #*** table objects
	if (type == "table" ){
		dat4 <- utils::read.table( file.path(dir,file) , header=TRUE , ... )
						}
    #*** sav objects (SPSS objects)
	if (type == "sav" ){
	  if ( ! spss.default){
			dat4 <- foreign::read.spss( file.path(dir,file) , ... )
							}
	  if (  spss.default){
			dat4 <- foreign::read.spss( file.path(dir,file) , 
				to.data.frame=TRUE , use.value.labels=FALSE , ... )
							}			
				}				
	return(dat4)
			}
#########################################################################			
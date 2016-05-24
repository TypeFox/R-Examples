

#########################################################################
# miceadds: saving data
save.data <- function( data , filename , type="Rdata" , path=getwd() , 
          row.names=FALSE , na = NULL , suffix = NULL , suffix_space = "__" , 
		  index = FALSE , systime = FALSE ,  ...){
	#*** the resulting object is dat4!	
	dir <- path
	file <- filename	
	if ( ! is.null(suffix) ){
		file <- paste0( file , suffix_space , suffix )
							}
	#*** add index in data frame if requested
	if ( index){ 
		data <- index.dataframe(data , systime=systime)
				}
	
	
	#*** missing handling
#	if ( is.null(na) ){
#		na <- base::switch( type , 
#					"csv" = "" , 
#					"csv2" = "" ,
#					"table" = "." )
#						}
#	type2 <- type
#	if ( type == "csv2" ){ 
#			type2 <- "csv" 
#					}
#	if ( type == "table" ){ 
#			type2 <- "dat" 
#				}
	# i1 <- grep( type2 , file )
	
	file0 <- file
				
    #*** Rdata objects	
	if ( "Rdata" %in% type  ){
	    file <- calc_filename( file=file0 , type="Rdata")
		base::save( data , file= file.path( dir , file ) )
				}
    #*** csv2 objects
	if ("csv2" %in% type  ){
	    if ( is.null(na)){ na <- "" }
		file <- calc_filename( file=file0 , type="csv2")
		utils::write.csv2( data , file.path(dir,file) , row.names=row.names , na=na ,... )
				}
    #*** csv objects
	if ("csv" %in% type ){
	    if ( is.null(na)){ na <- "" }
		file <- calc_filename( file=file0 , type="csv")
		utils::write.csv( data , file.path(dir,file) , row.names=row.names , na=na,  ... )
				}
    #*** table objects
	if ("table" %in% type ){
	    if ( is.null(na)){ na <- "." }
		file <- calc_filename( file=file0 , type="table")
		utils::write.table( data , file.path(dir,file) , na=na , row.names=row.names , ... )
				}
    #*** sav objects (SPSS objects)
	if ( "sav" %in% type ){
	    file <- calc_filename( file=file0 , type="sav")
		data <- sjmisc::set_label( data, lab=attr(data, "variable.labels") )
		sjmisc::write_spss( data , file.path( dir , file ) )		
				}				
			}
#########################################################################			
calc_filename <- function( file , type){
    type2 <- type 
	if ( type == "csv2" ){ 
			type2 <- "csv" 
					}
	if ( type == "table" ){ 
			type2 <- "dat" 
				}
	i1 <- grep( type2 , file )
	i1 <- grep( paste0("\\.",type2) , file , fixed=TRUE)	
	if ( length(i1) == 0 ){
		file <- paste0( file , "." , type2 )
						}
	return(file)
		}
########################################################		
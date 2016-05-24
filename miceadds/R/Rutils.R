

#****************************************
# load R data objects globally
load.Rdata <- function( filename , objname ){ 
    d1 <- load( filename )
    eval( parse( text=paste( objname ,  "<<- " , d1  ) ) )
            }

#**********************************************
# load Rdata objects			
load.Rdata2 <- function( filename , path=getwd() ){
    d1 <- load( file=file.path(path,filename) )
    objname <- "ma01"
	eval(parse(text = paste(objname, "<- ", d1)))	
    eval( parse(text= paste0( "return( " , objname , ")" ) ) )			
					}
			
#****************************************			
# save R data objects			
save.Rdata <- function( dat , name , path = NULL , part.numb = 1000 ){
    if ( is.null(path) ){ path <- getwd() }
    save( dat , file= file.path( path , paste( name , ".Rdata" , sep="") ) )
	M1 <- min( part.numb,nrow(dat) )
	if (M1== part.numb){ kuerzel <- "__PART" } else { kuerzel <- "" }
    write.csv2( dat[ seq(1, M1 ) , ] , 
                file= file.path( path , paste( name , kuerzel , ".csv" , sep="") ) ,
                quote=TRUE , row.names=F , na = "" )
		sink( file.path( path , paste( name , "__STR.Rout" , sep="") ) )	
	cat( paste( Sys.time() , "\n") )
	cat( Rsessinfo() , "\n\n")	
	cat( file.path( path , paste( name , ".Rdata" , sep="") ) , "\n\n")
	cat( utils::str( dat , list.len = ncol(dat) + 50 ))
		sink()
	dfr <- data.frame( "column" = 1:(ncol(dat)) , "variable" = colnames(dat) )
    write.csv2( dfr, 
                file= file.path( path , paste( name , "__VARIABLES.csv" , sep="") ) ,
                quote=TRUE , row.names=F , na = "" )	
                }
#****************************************


#**********************************************************
# function for sourcing all files within a folder
source.all <- function( path , grepstring= "\\.R" , print.source=TRUE ){
#        files <- list.files( path , grepstring ) 
        files <- list.files( path  ) 
		files <- grep.vec(  grepstring , files , "OR")$x
#		files <- files[ grep( "\\.R" , files ) ]
        for ( ff in files ){ 
			source( file.path( path , ff ) ) 
			if ( print.source ){ 
				cat( paste( "*** source" , ff ) , "\n") 
				utils::flush.console()
						}
				}
            }
#**********************************************************


#---------------------------------------------------------------------------
# stringtype <- "*" ; label <- "X->Y" ; rep.N <- 1 # replications of "\n"
output.format1 <- function( stringtype , label , rep.N = 1 , stringlength = 70){
    h1 <-  paste( rep("\n" ,  rep.N) , collapse="") 
    cat(h1)
    cat( paste( paste( rep( stringtype , stringlength ) , collapse ="") ) )
    cat(h1)
    cat( label )
    cat(h1)
	cat("\n")
    }
#---------------------------------------------------------------------------





#**************************************************
# evaluates an R function
Reval <- function( Rstring , print.string=TRUE){ 
    if (print.string){ cat( paste( "R>" , Rstring ) , "\n"  ) }
	eval.parent( parse( text = paste( Rstring )) , n=1 )	
	# .GlobalEnv
	# or eval.parent(expr, n = 1)
            }
#**************************************************
Revalpr <- function( Rstring , print.string=TRUE){
    if (print.string){ cat( paste( "R> print(" , Rstring , ")") , "\n"  ) }
	eval.parent( parse( text = paste( "print(",Rstring,")" )) , n=1 )	
		}
		
		
Revalprstr <- function( Rstring , print.string=TRUE){
    if (print.string){ cat( paste( "R> print(str(" , Rstring , "))") , "\n"  ) }
	eval.parent( parse( text = paste( "print(str(",Rstring,"))" )) , n=1 )	
		}

##################################################
# vector version of grep
grep.vec <- function( pattern.vec , x , operator="AND"){
	x0 <- x
	xv <- NULL
    for (vv in 1:(length(pattern.vec) ) ){
	if (operator == "AND"){
        x <- x[ grep( pattern.vec[vv] , x ) ]
					} else {
			xv <- union( xv ,x0[ grep( pattern.vec[vv] , x0 ) ] )
						}
                }
	if (operator!="AND"){ x <- xv }
		index.x <- which( x0 %in% x )
		res <- list( "x" = x , "index.x" = index.x )
        return(res)
        }
##################################################	
				
				
####################################################
# adds an index to a data frame
index.dataframe <- function (data, systime = FALSE) {
    data <- as.data.frame(data)
    data$index <- NULL
	data$file_created <- NULL
    data1 <- data.frame(index = seq(1, nrow(data)))
    if (systime) {
        data1$file_created <- paste0(" ", substring(Sys.time(), 1, 20)) 
		}
	data1 <- cbind( data1 , data )
    return(data1)
	}
########################################################	
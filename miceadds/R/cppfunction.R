cxxfunction.copy <- function( cppfct , name ){
#    code1 <- paste( readLines( codefile ) , collapse="\n      ")
#    fx <- cxxfunction( signature.input , code1 , plugin = plugin , verbose=verbose )
#    g1 <-  getDynLib(cppfct)
    g1 <-  inline::getDynLib(cppfct)
    cppname <- gsub( "\\.dll" , "\\.cpp" ,  g1[["path"]] )
    h1 <- readLines( cppname )
    tempname <- g1[["name"]]
    h1 <- gsub( tempname , name , h1 )
	h1 <- c( paste0( "//  Code created: " , Sys.time() ) , "" , h1 )
	name1 <- paste0( tolower(name) , ".cpp" )
    writeLines( h1 , name1 )
	crlrem( filename1=name1 , filename2=name1 )
            }

######################################################
# remove line endings
crlrem <- function( filename1 , filename2 ){
    filename <- filename1
    con <- file(filename, "rb")
    bin <- base::readBin(con, raw(), 100000)
    bin <- bin[ which(bin != "0d") ]
    close(con)
    Sys.sleep(1)
    con <- file(filename2, "wb")
    base::writeBin(bin, con)
    close(con)
        }
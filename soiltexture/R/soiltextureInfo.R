
# soiltextureInfo ====================================================

soiltextureInfo <- function(# Display and / or export system and package version information
### Display and / or export system and package version information. 
###  Can be used to provide an overview of the system and the R 
###  packages that were used to produce some calculations, thus 
###  improving the traceability of that work in the long run.

##seealso<< The base functions that were used internally to compile 
##  the information: \code{\link[base]{Sys.time}}, 
##  \code{\link[base]{Sys.info}}, \code{\link[base]{version}}, 
##  \code{\link[base]{.packages}}, 
##  \code{\link[utils]{installed.packages}}, 
##  \code{\link[tools]{package_dependencies}}. See also the 
##  \code{MD5} file in each package directory (and 
##  \code{\link[tools]{md5sum}} for generating these MD5 checksums).

 file = NULL, 
###  Single character string. Name of the text file (with or without 
###  its path) in which the information will be exported. If 
###  \code{NULL} (default), information are not exported.

    verbose  = TRUE, 
###  Single logical value. If \code{TRUE}, information are displayed 
###  on the screen.

    depends  = FALSE, 
###  Single logical value. If \code{TRUE}, information on packages 
###  dependencies are also displayed, in the same way

    md5      = TRUE,
###  Single logical value. If \code{TRUE}, the package MD5 checksums 
###  are returned too

    packages = "soiltexture" 
###  Single character string. Name of the package whose information 
###  must be returned.

){  
    #   depends <- TRUE; md5 <- TRUE; packages <- "soiltexture"; verbose <- TRUE
    
    tmpFile <- tempfile() 
    con     <- file( description = tmpFile, open = "w+" ) 
    
    #   In case the function crashes, the file will be closed
    on.exit( { close( con ); unlink( tmpFile ) } )
    
    cat2 <- function( ..., .file = con, sep = " " ){ 
        cat( ..., file = .file, sep = sep ) 
    }   
    
    cap2 <- function( ..., .file = con, .append = TRUE ){ 
        capture.output( ..., file = .file, append = .append ) 
    }   
    
    cat2( "INFORMATION ON SYSTEM AND PACKAGE(S) VERSION(S)\n" ) 
    cat2( "===============================================\n\n" ) 
    
    cat2( sprintf( "Date and time: %s\n\n", Sys.time() ) ) 
    
    cat2( "System info and R version\n" ) 
    cat2( "-------------------------\n\n" ) 
    
    #   System information
    cat2( "System info:\n" ) 
    
    sysInfo <- Sys.info() 
    sysInfo <- data.frame( 
        "info"  = names( sysInfo ), 
        "value" = as.character( sysInfo ), 
        stringsAsFactors = FALSE 
    )   
    
    cap2( sysInfo ); rm( sysInfo )
    cat2( "\n" ) 
    
    #   R version
    cat2( "R version:\n" ) 
    
    rversion <- data.frame( 
        "info"  = names( version ), 
        "value" = as.character( version ), 
        stringsAsFactors = FALSE 
    )   
    
    cap2( rversion ); rm( rversion )
    cat2( "\n" ) 
    
    #   R version
    cat2( "Loaded packages:\n" ) 
    
    loadedPackages    <- .packages() 
    
    # require( "utils" ) 
    installedPackages <- utils::installed.packages() 
    
    loadedPackages    <- installedPackages[ sort( loadedPackages ), 
        "Version" ] 
    loadedPackages    <- data.frame( 
        "package"  = names( loadedPackages ), 
        "version"  = as.character( loadedPackages ), 
        stringsAsFactors = FALSE 
    )   
    
    cap2( loadedPackages ) 
    cat2( "\n" ) 
    
    #   Package information
    # if( length( package ) > 1 ){ 
        # stop( "length( package ) > 1" )
    # }   
    
    cat2( sprintf( 
        "Information for packages: %s\n\n", 
        paste( packages, collapse = "; " )
    ) ) 
    
    
    testPack <- packages %in% installedPackages[, "Package" ]
    if( !all( testPack ) ){ 
        stop( sprintf( 
            "Can't find some package(s): %s", 
            paste( packages[ !testPack ], collapse = "; " )
        ) ) 
    };  rm( testPack ) 
    
    
    if( depends ){ 
        #   Find dependencies
        # require( "tools" ) 
        
        dep <- unlist( lapply( 
            X   = packages, 
            FUN = function(X){ 
                out <- tools::package_dependencies( 
                    packages  = packages, 
                    db        = installedPackages, 
                    recursive = TRUE ) 
                
                return( unlist( out ) )
            }   
        ) ) 
        dep <- sort( unique( dep ) ) 
        
        #   Remove base packages
        xPriority <- installedPackages[ dep, "Priority" ]
        xPriority[ is.na( xPriority ) ] <- ""
        dep <- dep[ xPriority != "base" ]
        rm( xPriority )
        dep <- dep[ !(dep %in% packages) ] 
        
        if( length( dep ) > 0 ){ 
            cat2( "Dependencies (except base packages):\n" ) 
            for( p in dep ){
                cat2( sprintf( "*   %s\n", p ) )
            }   
            cat2( "\n" ) 
            
            #   Add to the list of packages
            packages <- c( packages, dep )
        }else{ 
            cat2( "Dependencies (except base packages): none\n\n" ) 
        }   
        
        rm( dep )
    }   
    
    
    for( p in packages ){ 
        cat2( sprintf( "Package information for %s\n", p ) ) 
        cat2(          "-----------------------\n\n" ) 
        
        #   Package version
        cat2( sprintf( "Package version: %s\n\n", 
            installedPackages[ p, "Version" ] ) ) 
        
        #   Files in packages directory
        packFiles <- list.files( system.file( package = p ) )
        
        svnRev <- c( "SVN_VERSION", "SVN_REVISON", "REVISION" ) 
        svnRev <- svnRev[ svnRev %in% packFiles ]
        
        if( length( svnRev ) > 0 ){ 
            svnRev <- svnRev[ 1L ] 
            svnRev <- system.file( "SVN_VERSION", package = p )
            svnRev <- readLines( con = svnRev ) 
            cat2( "SVN revision: ", svnRev, "\n\n", sep = "" )
        };  rm( svnRev )
        
        
        if( "MD5" %in% packFiles ){ 
            md5  <- system.file( "MD5", package = p )
            .md5 <- readLines( con = md5 )
            cat2( "Package MD5: ", md5, "\n\n", sep = "" )
            cat2( .md5, sep = "\n" ) 
            cat2( "\n\n" )
            
            rm( md5, .md5 )
            
        }else{ 
            cat2( "No MD5 checksum found\n\n" )
        }   
        
        rm( packFiles )
        
    };  rm( p )
    
    close( con ); on.exit() 
    info <- readLines( con = tmpFile ) 
    
    if( verbose ){ cat( info, sep = "\n" ) }
    
    if( !is.null( file ) ){ writeLines( text = info, con = file ) }
    
    return( invisible( info ) ) 
###  Invisibly returns the information as a vector of character 
###  strings
}   

#   soiltextureInfo( file = "soiltextureInfo.txt", depends = TRUE ) 
#   soiltextureInfo( packages = "macrolegions" ) 



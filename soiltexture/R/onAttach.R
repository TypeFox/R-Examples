
.onAttach <- function(# Internal. Message displayed when loading the package.
 libname, 
 pkgname  
){  
    #   Try to load triangles and languages with special characters
    tryRes <- try( lang.par2 <- readRDS( file = system.file( 
        "languages.rds", package = pkgname ) ) ) 

    tryRes2 <- try( triPar2 <- readRDS( file = system.file( 
        "triangles.rds", package = pkgname ) ) )
    
    
    #   Extend language parameter if the polish triangle was loaded
    #   successfully
    if( !("try-error" %in% c( class( tryRes ), class( tryRes2 ) ) ) ){ 
        #   Save languages with special characters
        lang.par <- TT.get( "lang.par" ) 
        
        lang.par <- rbind( 
            lang.par, 
            lang.par2 
        )   
        
        TT.set( "lang.par" = lang.par ) 
        
        
        #   Save triangles with special characters
        for( i in 1:length( triPar2 ) ){ 
            do.call( what = "TT.add", args = triPar2[ i ] )
        }   
        
    }else{ 
        packageStartupMessage( "soiltexture: Languages and triangles with special characters could not be loaded" )
    }   
    
    
    # Welcome message
    if( interactive() ){ 
        svnVersion <- system.file( "SVN_VERSION", package = pkgname ) 
        
        if( svnVersion != "" ){ 
            svnVersion <- readLines( con = svnVersion )[ 1L ] 
            svnVersion <- sprintf( "(svn revision: %s)", svnVersion ) 
        }else{ 
            svnVersion <- "(svn revision: ?)" 
        }   
        
        msg <- sprintf( 
            "%s %s %s. For help type: help(pack='%s')", 
            pkgname, 
            as.character( packageVersion( pkgname ) ), 
            svnVersion, 
            pkgname ) 
        
        packageStartupMessage( msg ) 
    }   
}   



# .select.list2 
# ==========================================================

.select.list2 <- function( # Wrapper around 'menu' with error handling
### Wrapper around 'menu' with error handling

##seealso<< \code{\link[utils]{select.list}}, on which 
##  \code{.select.list2} is based.

##keywords<< internal 

    title = NULL, 
###  See \code{\link[utils]{select.list}}.

    choices, 
###  See \code{\link[utils]{select.list}}.

    graphics = FALSE, 
###  See \code{\link[utils]{select.list}}.

    preselect = NULL, 
###  See \code{\link[utils]{select.list}}.

    error = "You haven't chosen anything", 
###  Single character string. Error message to be displayed if 
###  the user does not chose any item (code 0).

    multi = FALSE
###  Single logical. If \code{TRUE}, then multiple choices are 
###  allowed.
){  
    #   Index of choices
    choicesNum <- 1:length(choices) 
    
    #   Label choice index
    names( choicesNum ) <- choices 
    
    # library( "utils" )
    
    mRes <- select.list( 
        title       = title,
        choices     = choices, 
        preselect   = preselect, 
        multiple    = multi, 
        graphics    = graphics 
    )   
    
    ## Error handling:
    if( length(mRes) == 0 ){ 
        stop( error ) 
    }   
    
    mRes <- choicesNum[ mRes ] 
    names( mRes ) <- NULL 
    
    if( any( is.na( mRes ) ) ){ 
        stop( "Wrong value(s) chosen" )
    }   
    
    return( mRes ) 
###  The \bold{index} of the user's choice (index from 
###  \code{choices}).
}   

    # .select.list2(
        # title     = "Please chose", 
        # choices   = c( "A", "B" ) 
    # )   



# .chooseTextFiles
# ==========================================================

.chooseTextFiles <- function( # Pop-up a menu to choose one or several text file(s) from the file system.
### Pop-up a menu to choose one or several text file(s) from the file system.

##keywords<< internal 

    caption = "Select a text file", 
###   See \code{\link[utils]{choose.files}} or 
###   \code{\link[tcltk]{tk_choose.files}}.

    multi   = FALSE 
###   See \code{\link[utils]{choose.files}} or 
###   \code{\link[tcltk]{tk_choose.files}}.

){  
    if( !interactive() ){ 
        stop( "'.chooseTextFiles' can only be used in interactive mode" )
    }   
    
    
    ## Create a template of file extension to be read:
    filterz <- matrix( 
        data  = c( 
            "Text file (*.txt)", "*.txt", 
            "CSV file (*.csv)", "*.csv" ), 
        nrow  = 2, 
        ncol  = 2, 
        byrow = TRUE  
    )   
    rownames( filterz ) <- c( "txt", "csv" ) 
    
    ## Pop-up a menu to choose the bin file to be 
    ## imported
    if( exists( "choose.files", where = "package:utils" ) ){ 
        # fun <- get( "choose.files" ) 
        
        .file <- utils::choose.files(
            # default = "", # , "*.bin"
            caption = caption, 
            multi   = multi, 
            filters = filterz 
        )   
        
    }else{ 
        .file <-tcltk::tk_choose.files(
            # default = "", # , "*.bin"
            caption = caption, 
            multi   = multi, 
            filters = filterz 
        )   
    }   
    
    if( length( .file ) == 0 ){ 
        stop( "No file was chosen" ) 
    }   
    
    return( .file ) 
### Returns the full path to the selected file.
}   



# .read.table.menu
# ==========================================================

.read.table.menu <- function( # Text-based menu for importing a text table (typically a CSV file)
### Text-based menu for importing a text table (typically a 
###  CSV file).

##keywords<< internal 

    graphics = FALSE, 
###  See \code{\link[utils]{select.list}}.

    stringsAsFactors = FALSE,
###  See \code{\link{read.table}}.
    
    ... 
###  Additional parameters passed to 
### \code{\link[utils]{read.table}}.

){  
    message( "==== Select the file to import ====" ) 
    # ==============================================
    
    message( "* Provide a text file (.txt or csv) containing soil texture data" ) 
    message( "* This must be a text file containing tabular data" ) 
    message( "* The file *must* contain the following headers: CLAY, SILT, SAND" ) 
    message( "* The column order does not matter, and other variables can be present" ) 
    
    f <- .chooseTextFiles(
        caption  = "* Select a text file with soil texture data", 
        multi    = FALSE  
    )   
    
    
    
    message( "==== File format ====" ) 
    # ================================
    
    #   Field separator character
    sep <- c( 
        "Comma (,)"           = ",", 
        "Semi-colon (;)"      = ";", 
        "Tabulation (\t)"     = "\t", 
        "Single space ( )"    = " ", 
        "Multiple space (  )" = "" 
    )   
    
    sep <- sep[ .select.list2(
        title     = "* Choose the field/column separator character", 
        choices   = names( sep ), 
        graphics  = graphics, 
        preselect = names( sep )[ 1L ], 
        error     = "You haven't chosen anything", 
        multi     = FALSE
    ) ] 
    
    
    #   'decimal points' character
    dec <- c( 
        "Dot (.)"    = ".", 
        "Comma  (,)" = ","  
    )   
    
    dec <- dec[ .select.list2(
        title     = "* Choose the 'decimal points' character", 
        choices   = names( dec ), 
        graphics  = graphics, 
        preselect = names( dec )[ 1L ], 
        error     = "You haven't chosen anything", 
        multi     = FALSE
    ) ] 
    
    
    #   'decimal points' character
    fileEncoding <- c( 
        "Internal (current locale)" = "native.enc", 
        "UTF-8"                     = "UTF-8"  
    )   
    
    fileEncoding <- fileEncoding[ .select.list2(
        title     = "* Choose the encoding of the file", 
        choices   = names( fileEncoding ), 
        graphics  = graphics, 
        preselect = names( fileEncoding )[ 1L ], 
        error     = "You haven't chosen anything", 
        multi     = FALSE
    ) ] 
    
    
    message( "==== Importing the file ====" ) 
    # =======================================
    
    message( sprintf( "* File: %s", f ) ) 
    
    dat <- read.table( file = f, sep = sep, dec = dec, 
        stringsAsFactors = stringsAsFactors, 
        fileEncoding = fileEncoding, header = TRUE )
    
    message( sprintf( "* Dimensions: %s columns, %s rows", 
        ncol(dat), nrow(dat) ) )  
    
    message( sprintf( "* Columns: %s", 
        paste( colnames( dat ), collapse = "," ) ) )  
    
    
    #   Save attributes
    attr( x = dat, which = "sep" ) <- sep 
    attr( x = dat, which = "dec" ) <- dec 
    attr( x = dat, which = "fileEncoding" ) <- fileEncoding 
    
    return( dat ) 
###  A \code{\link{data.frame}} with the imported soil texture 
###  data (column names \code{CLAY}, \code{SILT} and \code{SAND}, 
###  along with any additional column present in the file 
###  imported).
}   

# library( "soiltexture" )
# .read.table.menu()



# soiltexture_gui
# ==========================================================

soiltexture_gui <- function( # Text-based menu for plotting and classifying soil texture data
### Text-based menu for plotting and classifying soil texture 
### data. 
### 
### If you simply want to obtain a figure with an 
### empty soil texture triangle, just call 
### \code{soiltexture_gui}() and follow the instructions.
### 
### If you want to a figure with your own soil texture data 
### on top of a texture triangle, you must first prepare 
### a tabular text file containing your texture data, as 
### \code{.txt} or \code{.csv}. Such a file can be prepared 
### with \code{MS Excel} or \code{Libre Office}, and exported 
### as CSV ("CSV (comma delimited) (*.csv)" or "CSV (MS-DOS) 
### (*.csv)" for example). The table \bold{must} contain 
### headers (column names) and it \bold{must} the following 
### columns and headers: \code{CLAY}, \code{SILT} and 
### \code{SAND}. Other columns are allowed and will be ignored. 
### In the texture data file, each row represent a record 
### (a sample) and each column a variable.
### 
### You will be asked about the format of this text file, in 
### particular about the field / column separator (it can be 
### commas, semi-colons, tabulations or (multiple) spaces) 
### and the decimal mark (comma or dot). The file encoding 
### can be either the native encoding of the computer, or 
### UTF-8 (without BOM).
### 
### The sum of the texture of each row must be either 1 (if 
### expressed as a fraction) or 100 (if expressed as a 
### percentage). You will be asked about the unit. Only small 
### divergences from 1 or 100 are allowed, but you will be 
### asked if you want to normalise your data beforehand, so 
### larger divergences are possible.
###
### You will also be asked which texture classification system 
### you want to use (FAO, USDA, etc.). It is possible to 
### plot a texture triangle without texture classification.
### 
### Finally, if you have chosen a texture classification system, 
### \code{soiltexture_gui} can classify each record according 
### to this classification system and 
### \bold{return you the texture class of each record}, 
### as a CSV text file.
###
### The texture triangle is show to you with R default 
### graphical device, and you can choose to export a 
### PNG figure of the resulting texture triangle (between 
### 512 and 2048 pixel width/height, depending on what you 
### chose).
### 

    main = NULL, 
###  Single character string. Main title of the texture 
###  diagram. Set to \code{NA} to obtain a a slightly bigger 
###  figure, with no title. See \code{\link[soiltexture]{TT.plot}}.

    graphics = FALSE, 
###  See \code{\link[utils]{select.list}}.

    ...
###  Additional parameters passed to 
###  \code{soiltexture:::.read.table.menu} or 
###  (subsequently) to \code{\link[utils]{read.table}}.

){  
    message( "+-------------------------------+" ) 
    message( "|    The Soil Texture Wizard    |" ) 
    message( "+-------------------------------+" ) 
    
    message( "Text-based interface for plotting and classifying soil texture data" ) 
    
    if( interactive() ){ 
        #   Field separator character
        fileImport <- c( 
            "Import soil texture data"      = TRUE, 
            "Plot an empty texture diagram" = FALSE 
        )   
        
        fileImport <- fileImport[ .select.list2(
            title     = "* Do you want to:", 
            choices   = names( fileImport ), 
            graphics  = graphics, 
            preselect = names( fileImport )[ 1L ], 
            error     = "You haven't chosen anything", 
            multi     = FALSE
        ) ] 
        
        # Import the soil texture file
        # ============================
        
        if( fileImport ){
            dat <- .read.table.menu( graphics = graphics, ... )
            
            fileEncoding <- attr( x = dat, which = "fileEncoding" ) 
            dec          <- attr( x = dat, which = "dec" ) 
            sep          <- attr( x = dat, which = "sep" ) 
            
            
            # Data control
            # ============
            
            css <- c( "CLAY", "SILT", "SAND" )
            
            #   Test if the required column are present:
            testCol <- css %in% colnames( dat )
            
            if( any( !testCol ) ){
                stop( sprintf( 
                    "Some required column(s) is(are) missing: %s", 
                    paste( css[ !testCol ], collapse = "; " ) 
                ) ) 
            };  rm( testCol )
            
            
            #   Test that the required columns are numeric:
            for( cl in css ){
                if( !is.numeric( dat[, cl ] ) ){
                    stop( 
                        sprintf( 
                            "Wrong format for column %s. Expected 'numeric' or 'integer', got %s. ", 
                            cl, paste( dat[, cl ], collapse = "; " ) 
                        ),  
                        "Can be a problem with field separator or decimal mark"
                    ) 
                }   
            };  rm( cl )
            
            
            
            message( "==== Data normalisation ====" ) 
            # =======================================
            
            #   Unit for soil texture
            text.sum <- c( 
                "Percentage (0-100%)" = 100, 
                "Fraction  (0-1)"     = 1  
            )   
            
            text.sum <- text.sum[ .select.list2(
                title     = "* Choose the unit for soil texture data", 
                choices   = names( text.sum ), 
                graphics  = graphics, 
                preselect = names( text.sum )[ 1L ], 
                error     = "You haven't chosen anything", 
                multi     = FALSE
            ) ] 
            
            message( sprintf( 
                "* For each row, the sum CLAY + SILT + SAND *must* be equal to %s%s", 
                text.sum, ifelse( text.sum == 100, "%", "" ) 
            ) ) 
            
            #   Unit for soil texture
            normalise <- c( 
                "Yes" = TRUE, 
                "No"  = FALSE 
            )   
            
            normalise <- normalise[ .select.list2(
                title     = sprintf( 
                    "* Should the Clay+Silt+Sand sums be normalised to %s%s", 
                    text.sum, ifelse( text.sum == 100, "%", "" ) 
                ),  
                choices   = names( normalise ), 
                graphics  = graphics, 
                preselect = names( normalise )[ 1L ], 
                error     = "You haven't chosen anything", 
                multi     = FALSE
            ) ] 
            
            if( normalise ){
                message(  "* Normalising the texture data" )
                
                dat[, css ] <- TT.normalise.sum( 
                    tri.data = dat[, css ], 
                    text.sum = text.sum ) 
            }   
        }else{
            dat <- NULL # No data was imported
        }   
        
        
        
        message( "==== Classification system ====" ) 
        # ==========================================
        
        #   Create a list of possible texture triangles 
        class.sys.list <- c(
            #"none"             = "none",
            # "FAO"             = "FAO50.TT", 
            "HYPRES"            = "HYPRES.TT", 
            "USDA (US)"         = "USDA.TT",
            "Aisne (France)"    = "FR.AISNE.TT",
            "GEPPA (France)"    = "FR.GEPPA.TT",
            "BK94 (Germany)"    = "DE.BK94.TT",
            "SEA74 (Germany)"   = "DE.SEA74.TT", 
            "TGL85 (Germany)"   = "DE.TGL85.TT", 
            "SSEW (UK)"         = "UK.SSEW.TT",
            #"Australia"        = "AU.TT", 
            "Australia"         = "AU2.TT",
            "Belgium"           = "BE.TT",
            "Canadian (fr)"     = "CA.FR.TT", 
            "Canadian (en)"     = "CA.EN.TT", 
            "ISSS"              = "ISSS.TT", 
            "Romania"           = "ROM.TT", 
            "Poland"            = "PL.TT", 
            "Brasil (1996)"     = "BRASIL.TT",
            "Brasil (2013)"     = "SiBCS13.TT", 
            "USDA 1911 (US)"    = "USDA1911" 
        )   
        
        #   Eliminates diagrams that are not found in current 
        #   version of soiltexture
        class.sys.list <- class.sys.list[ class.sys.list %in% names( TT.get() ) ]
        
        #   Add class none
        class.sys.list <- c( "none" = "none", class.sys.list ) 
        
        class.sys.list <- class.sys.list[ .select.list2(
            title     = "* Choose the texture classification system",  
            choices   = names( class.sys.list ), 
            graphics  = graphics, 
            preselect = names( class.sys.list )[ 1L ], 
            error     = "You haven't chosen anything", 
            multi     = FALSE
        ) ] 
        
        
        
        message( "==== Control plot ====" ) 
        # =================================
        
        plotTexture <- function( 
            .main      = main, 
            .tri.data  = dat, 
            .class.sys = class.sys.list, 
            ... 
        ){  
            TT.plot( tri.data = .tri.data, class.sys = .class.sys, 
                main = .main, family.op = "serif", ... )
        }   
        
        plotTexture( cex.axis = 1, cex.lab = 1, cex.main = 1 ) 
        
        
        
        message( "==== Export PNG figure ====" ) 
        # =================================
        
        #   Export figure of 
        exportFigure <- c( 
            "Yes" = TRUE, 
            "No"  = FALSE 
        )   
        
        exportFigure <- exportFigure[ .select.list2(
            title     = "* Export a PNG figure of the diagram?",  
            choices   = names( exportFigure ), 
            graphics  = graphics, 
            preselect = names( exportFigure )[ 1L ], 
            error     = "You haven't chosen anything", 
            multi     = FALSE
        ) ] 
        
        
        if( exportFigure ){
            message( "* Note: The figure will look a bit different the default graphical device" )
            
            
            #   Find out the path of the temporary directory
            d <- tempdir() 
            
            
            #   Choice of figure size
            .size <- c( 
                "512 px"  = 512, 
                "1024 px" = 1024, 
                "2048 px" = 2048
            )   
            
            .size <- .size[ .select.list2(
                title     = "* Choose the size of the PNG figure:",  
                choices   = names( .size ), 
                graphics  = graphics, 
                preselect = names( .size )[ 1L ], 
                error     = "You haven't chosen anything", 
                multi     = FALSE
            ) ] 
            
            #   Find out the best resolution
            if( .size == 512 ){
                .res <- 50
            }else if( .size == 1024 ){
                .res <- 100
            }else if( .size == 2048 ){
                .res <- 150
            }else{
                .res <- 50
            }   
            
            
            #   Clean classification system name
            cs <- gsub( x = class.sys.list, pattern = ".TT", 
                replacement = "", fixed = TRUE ) 
            cs <- gsub( x = cs, pattern = ".", replacement = "-", 
                fixed = TRUE ) 
            
            #   File name and path
            filename <- sprintf( 
                "texture_%s_%s.png", 
                cs, 
                format( Sys.time(), "%Y%m%d-%H%M%S" ) 
            )
            
            filename <- file.path( d, filename )
            
            
            png(
                filename = filename, 
                width    = .size, 
                height   = .size, 
                res      = .res 
            )   
            
                plotTexture() 
            
            dev.off() 
                
            message( sprintf( "* Figure exported in: %s", filename ) )
        }   
        
        
        
        if( (!is.null( dat )) & (class.sys.list != "none") ){
            message( "==== Classify texture data ====" ) 
            # ==========================================
            
            #   Export figure of 
            classifData <- c( 
                "Yes" = TRUE, 
                "No"  = FALSE 
            )   
            
            classifData <- classifData[ .select.list2(
                title     = "* Classify the soil texture data?",  
                choices   = names( classifData ), 
                graphics  = graphics, 
                preselect = names( classifData )[ 1L ], 
                error     = "You haven't chosen anything", 
                multi     = FALSE
            ) ] 
            
            if( classifData ){
                #   Find out the path of the temporary directory
                d <- tempdir() 
                
                out1 <- TT.points.in.classes( 
                    tri.data  = dat, 
                    class.sys = class.sys.list, 
                    PiC.type  = "l" ) 
                
                cn <- colnames( out1 )
                
                out2 <- unlist( lapply(
                    X   = 1:nrow(out1), 
                    FUN = function(i){
                        out <- cn[ out1[ i, ] ]
                        if( length( out ) == 0 ){
                            out <- ""
                        }else{ 
                            out <- out[ 1L ] # 1st value is taken
                        }   
                        
                        return( out )
                    }   
                ) ) 
                
                dat <- data.frame( 
                    dat,          # Original data
                    "firstClass" = out2, # Simple classification
                    out1                 # full classification 
                )   
                
                rm( out2, out1 )
                
                
                #   File name and path
                filename <- sprintf( 
                    "texture_%s_%s.csv", 
                    cs, 
                    format( Sys.time(), "%Y%m%d-%H%M%S" ) 
                )
                
                filename <- file.path( d, filename )
                
                write.table( 
                    x     = dat, 
                    file  = filename, 
                    quote = TRUE, 
                    sep   = sep, 
                    dec   = dec, 
                    fileEncoding = fileEncoding, 
                    row.names    = FALSE, 
                    col.names    = TRUE ) 
                
                message( sprintf( "* Classification exported in: %s", filename ) )
                message( "* Note: The exported table has the same format as the imported one (sep, dec, encoding)" ) 
            }   
        }else{
            classifData <- FALSE 
        }   
        
        
        
        if( exportFigure | classifData ){
            # shell.exec( tempdir() ) 
            browseURL( url = tempdir() ) 
        }   
    }else{ # Not interactive
        dat <- NULL 
        
        message( "R is not run interactively. Operation cancelled. See ?interactive" ) 
    }   
    
    return( invisible( dat ) ) 
###  Either \code{NULL} if no texture data was imported, 
###  or a \code{\link{data.frame}} (if texture data was 
###  imported). The texture classification is also returned 
###  (when the user asked for a texture classification).
}   


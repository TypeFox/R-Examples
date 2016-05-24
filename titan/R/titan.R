# TITAN
# TITration ANalyser to interpolate concentrations from Mass Spectrometry data
# created using R 2.4.0 www.cran.r-project.org
# (c) TS Price 2006


.packageName <- "titan"

titan <-
function(
    data =       NULL,
    trace =      TRUE,
    widget =     TRUE,
    dataFile =   "",
    outFile =    "",
    pdfFile =    "",
    flagRaw =    FALSE,
    flagFitted = FALSE,
    freqLo =     .05,
    freqHi =     .95,
    reg =        "least.squares",
    term =       "linear",
    sel =        "wald",
    alpha =      .05,
    rx0 =        NULL,
    rx1 =        NULL,
    gene0 =      NULL,
    gene1 =      NULL,
    R =          1000,
    seed =       0,
    ciConf =     .95,
    ciType =     "all"
)
{

titan.version <- "TITAN v.1.0-13: Titration Analysis"

# Trim call from error report
trim.error <-
function( err ) { sub( ".*:\\s*(.*?)\\s*$", "\\1", err, perl = TRUE ) }




###########################################################################





log10it <-
function( x ) { log10( x / ( 1 - x ) ) }

exp10it <-
function( y ) { x <- 10^y; return( x / ( 1 + x ) ) }

# degree 1 polynomial
# d1 = ( x - p$alpha[ 1 ] ) / sqrt( p$norm2[ 3 ] )
# => x - p$alpha[ 1 ] - sqrt( p$norm2[ 3 ] ) * d1 == 0
#
# degree 2 polynomial
# d2 = ( ( x - p$alpha[ 1 ] ) * ( x - p$alpha[ 2 ] ) - p$norm2[ 3 ] / p$norm2[ 2 ] ) / sqrt( p$norm2[ 4 ] )
# => x ^ 2 - ( p$alpha[ 1 ] + p$alpha[ 2 ] ) * x + p$alpha[ 1 ] * p$alpha[ 2 ] - p$norm2[ 3 ] / p$norm2[ 2 ] - sqrt( p$norm2[ 4 ] ) * d2 == 0

# equivalent to poly( x, degree =  2, coefs =  p )[ , 1 ]
polyL <-
function( x, p ) ( x - p$alpha[ 1 ] ) / sqrt( p$norm2[ 3 ] )

# equivalent to poly( x, degree =  2, coefs =  p )[ , 2 ]
polyQ <-
function( x, p ) ( ( x - p$alpha[ 1 ] ) * ( x - p$alpha[ 2 ] ) -
    p$norm2[ 3 ] / p$norm2[ 2 ] ) / sqrt( p$norm2[ 4 ] )

titan.widget <-
function( widget, ... )
{
    tkblank <- function( tt ) tklabel( tt, text =  "     " )

    tkblankline <- function( tt ) tkgrid( tkblank( tt ) )

    up <-
    function( str, index )
    {
        str[ index ] <- toupper( str[ index ] )
        return( str )
    }

    checkButtons <-
    function( tt, value, label, vmax = NULL )
    {
        cb <- cbVal <- cbLabel <- list()
        len <- length( value )
        if ( is.null( vmax ) )
        {
            vmax <- len
        }
        valueNames <- names( value )
        value <- as.character( as.numeric( value ) )
        for ( i in 1:len ) {
            cbVal[[ i ]] <- tclVar( value[ i ] )
            cb[[ i ]] <- tkcheckbutton( tt, variable = cbVal[[ i ]] )
            cbLabel[[ i ]] <- tklabel( tt, text = valueNames[ i ] )
        }
        # print check buttons in columns of `nr' buttons
        nr <- min( len, vmax )
        for ( i in 1:nr )
        {
            args <- list()
            args[[ 1 ]] <- tkblank( tt )
            if ( i == 1 )
            {
                args[[ 2 ]] <- tklabel( tt, text=  label )
            }
            else
            {
                args[[ 2 ]] <- tkblank( tt )
            }
            args[[ 3 ]] <- tkblank( tt )
            for ( j in 0:( ( len - 1 ) %/% nr ) )
            {
                if ( ( p <- i + j * nr ) > len )
                {
                    next
                }
                args[[ j * 4 + 4 ]] <- cb[[ p ]]
                args[[ j * 4 + 5 ]] <- cbLabel[[ p ]]
                args[[ j * 4 + 6 ]] <- tkblank( tt )
                args[[ j * 4 + 7 ]] <- tkblank( tt )
            }
            args$sticky <-  "w"
            do.call( "tkgrid", args )
        }
        for ( i in 1:len )
        {
            tkgrid.configure( cb[[ i ]], sticky = "e" )
        }
        return( list( cb = cb, value = cbVal, label = cbLabel ) )
    }

    checkButtons.nCols <-
    function( tt, value, label, vmax = NULL )
    {
        cb <- cbVal <- cbLabel <- list()
        len <- dim( value )[ 1 ]
        n <- dim( value )[ 2 ]
        if ( is.null( vmax ) )
        {
            vmax <- len
        }
        valueNames <- rownames( value )
        colNames <- colnames( value )
        for ( j in 1:len )
        {
            cb[[ j ]] <- list()
            cbVal[[ j ]] <- list()
        }
        for ( j in 1:len )
        {
            for ( k in 1:n )
            {
                cbVal[[ j ]][[ k ]] <- tclVar( as.character( as.numeric(
                value[ j, k ] ) ) )
                cb[[ j ]][[ k ]] <- tkcheckbutton( tt,
                    variable = cbVal[[ j ]][[ k ]] )
            }
            cbLabel[[ j ]] <- tklabel( tt, text = valueNames[ j ] )
        }

        # print column headings
        if ( !is.null( colnames( value ) ) )
        {
            args <- list()
            for ( j in 0:( ( len - 1 ) %/% vmax ) )
            {
                m <- j * ( n + 2 ) + 2
                args[[ m - 1 ]] <- tkblank( tt )
                args[[ m ]] <- tkblank( tt )
                for ( k in 1:n )
                {
                    args[[ m + k ]] <- tklabel( tt, text = colNames[ k ] )
                }
                args[[ m + k + 1 ]] <- tkblank( tt )
                args[[ m + k + 2 ]] <- tkblank( tt )
            }
            do.call( "tkgrid", args )
        }
        # print check buttons in columns of `nr' buttons
        nr <- min( len, vmax )
        for ( i in 1:nr )
        {
            args <- list()
            args[[ 1 ]] <- tkblank( tt )
            if ( i == 1 )
            {
                args[[ 2 ]] <- tklabel( tt, text = label )
            }
            else
            {
                args[[ 2 ]] <- tkblank( tt )
            }
            for ( j in 0:( ( len - 1 ) %/% nr ) )
            {
                if ( ( p <- i + j * nr ) > len )
                {
                    next
                }
                m <- j * ( n + 2 ) + 2
                for ( k in 1:n )
                {
                    args[[ m + k ]] <- cb[[ p ]][[ k ]]
                }
                args[[ m + k + 1 ]] <- cbLabel[[ p ]]
                args[[ m + k + 2 ]] <- tkblank( tt )
            }
            args$sticky <- "w"
            do.call( "tkgrid", args )
        }
        return( list( cb = cb, value = cbVal, label = cbLabel ) )
    }

    checkButtonValues <-
    function( cb, cbNames, col = NULL )
    {
        len <- length( cb$value )
        res <- logical( len )
        names( res ) <- cbNames
        if ( is.null( col ) )
        {
            for ( i in 1:len )
            {
                res[ i ] <- as.logical( as.numeric( tclvalue( cb$value[[ i ]] ) ) )
            }
        }
        else
        {
            for ( i in 1:len )
            {
                res[ i ] <- as.logical( as.numeric( tclvalue( cb$value[[ i ]][[ col ]] ) ) )
            }
        }
        return( res )
    }

    # First input widget
    widget1 <-
    function( opt, data )
    {
        # initialise variables
        regStr <- c( "Least squares polynomial", "Robust polynomial", "Least squares natural spline" )
        termStr <- list()
        termStr[[ 1 ]] <- termStr[[ 2 ]] <- c( "Linear", "Parallel Linear", "Quadratic" )
        termStr[[ 3 ]] <- paste( 2:4, "degrees of freedom" )
        selStr <- list( c( "Backwards by Wald test", "Stepwise by AIC" ), c( "Backwards by Wald test", "" ), c( "Stepwise by AIC", "" ) )
        regOptions <- tclVar( up( regStr, opt$reg ) )
        termOptions <- tclVar( up( termStr[[ 1 ]], opt$term ) )
        selOptions <- tclVar( up( selStr[[ 1 ]], opt$sel ) )
        regOpt <- c( "least.squares", "robust", "spline" )
        termOpt <- list()
        termOpt[[ 1 ]] <- termOpt[[ 2 ]] <- make.names( tolower( termStr[[ 1 ]] ) )
        termOpt[[ 3 ]] <- paste( 2:4, "df", sep = "." )
        selOpt <- list( c( "wald", "aic" ), c( "wald", "" ), c( "aic", "" ) )

        reg <- match( opt$reg, regOpt )
        tclvalue( regOptions  ) <- up( regStr, reg )
        tclvalue( termOptions ) <- up( termStr[[ reg ]], match( opt$term, termOpt[[ reg ]] ) )
        tclvalue( selOptions  ) <- up( selStr[[ reg ]], match( opt$sel, selOpt[[ reg ]] ) )

        # set up widget
        eventTcl <- tclVar( "NULL" )
        dataFile <- tclVar( opt$dataFile )
        outFile <- tclVar( opt$outFile )
        pdfFile <- tclVar( opt$pdfFile )
        flagRawTcl <- tclVar( opt$flagRaw )
        flagFittedTcl <- tclVar( opt$flagFitted )
        freqLoTcl <- tclVar( opt$freqLo )
        freqHiTcl <- tclVar( opt$freqHi )
        alphaTcl <- tclVar( opt$alpha )
        tt <- tktoplevel()
        tkwm.title( tt, titan.version )
        cbFlagRaw <- tkcheckbutton( tt, variable = flagRawTcl )
        cbFlagFitted <- tkcheckbutton( tt, variable = flagFittedTcl )
        entryFreqLo <- tkentry( tt, width = "8", textvariable = freqLoTcl )
        entryFreqHi <- tkentry( tt, width = "8", textvariable = freqHiTcl )
        tkblankline( tt )
        if ( is.null( data ) )
        {
            blank1 <- tkblank( tt )
            entry1 <- tkentry( tt, width = "75", textvariable = dataFile )
            onBrowse1 <- function() { tclvalue( dataFile ) <- tclvalue( tkgetOpenFile( initialdir = getwd() ) ) }
            browseBut1 <- tkbutton( tt, text = "  Browse  ", command = onBrowse1 )
            tkgrid(
                tkblank( tt ),
                tklabel( tt, text = "Enter MS data file name:" ),
                tkblank( tt ),
                entry1,
                browseBut1,
                blank1,
                sticky = "w"
            )
            tkgrid.configure( entry1, columnspan = 3 )
            tkgrid.configure( browseBut1, column = 7 )
            tkgrid.configure( blank1, column = 8 )
        }
        else
        {
            blank1 <- tkblank( tt )
            tkgrid(
                tkblank( tt ),
                tklabel( tt, text = "MS data:" ),
                tkblank( tt ),
                tklabel( tt, text = "Supplied to function" ),
                blank1,
                sticky = "w"
            )
            tkgrid.configure( blank1, column = 8 )
        }
        tkblankline( tt )
        onBrowse2 <- function() { tclvalue( outFile ) <- tclvalue( tkgetSaveFile( initialdir = getwd() ) ) }
        browseBut2 <- tkbutton( tt, text = "  Browse  ", command = onBrowse2 )
        entry2 <- tkentry( tt, width = "75", textvariable = outFile )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Enter text output file name:" ),
            tkblank( tt ),
            entry2,
            browseBut2,
            sticky = "w"
        )
        tkgrid.configure( entry2, columnspan = 3 )
        tkgrid.configure( browseBut2, column = 7 )
        tkblankline( tt )
        onBrowse3 <- function() { tclvalue( pdfFile ) <- tclvalue( tkgetSaveFile( initialdir = getwd() ) ) }
        browseBut3 <- tkbutton( tt, text = "  Browse  ", command = onBrowse3 )
        entry3 <- tkentry( tt, width = "75", textvariable = pdfFile )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Enter graphical output file name:" ),
            tkblank( tt ),
            entry3,
            browseBut3,
            sticky = "w"
        )
        tkgrid.configure( entry3, columnspan = 3 )
        tkgrid.configure( browseBut3, column = 7 )
        tkblankline( tt )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Flag" ),
            cbFlagRaw,
            tklabel( tt, text = "Data" ),
            entryFreqLo,
            tklabel( tt, text = "Minimum" ),
            sticky = "w"
        )
        tkgrid(
            tkblank( tt ),
            tkblank( tt ),
            cbFlagFitted,
            tklabel( tt, text = "Fitted values" ),
            entryFreqHi,
            tklabel( tt, text = "Maximum" ),
            sticky = "w"
        )
        tkgrid.configure( entryFreqLo, sticky = "e" )
        tkblankline( tt )
        regList <- tklistbox(
            tt,
            height = 3,
            width = 35,
            selectmode = "single",
            background = "white",
            activestyle = "dotbox",
            relief = "sunken",
            listvariable = regOptions
        )
        termList <- tklistbox(
            tt,
            height = 3,
            width = 30,
            selectmode = "single",
            background = "white",
            activestyle = "dotbox",
            relief = "sunken",
            listvariable = termOptions
        )
        selList <- tklistbox(
            tt,
            height = 2,
            width = 35,
            selectmode = "single",
            background = "white",
            activestyle = "dotbox",
            relief = "sunken",
            listvariable = selOptions
        )
        alphaEntry <-  tkentry(
            tt,
            width = "8",
            relief = "sunken",
            textvariable = alphaTcl,
            background = "white"
        )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Regression method:" ),
            tkblank( tt ),
            regList,
            tkblank( tt ),
            termList,
            sticky = "w"
        )
        tkgrid(
            tklabel( tt, text = "Alpha:" ),
            sticky = "w",
            column = 5
        )
        okBut <- tkbutton(
            tt,
            text = "  OK  ",
            command = function() { tclvalue( eventTcl ) <- "OK" }
        )
        cancelBut <- tkbutton(
            tt,
            text="  Cancel  ",
            command = function() { tclvalue( eventTcl ) <- "Cancel" }
        )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Selection method:" ),
            tkblank( tt ),
            selList,
            tkblank( tt ),
            alphaEntry,
            okBut,
            cancelBut,
            sticky = "w"
        )
        tkblankline( tt )
        tkfocus( tt )

        # bind events
        tkbind( tt, "<Return>",  function() tclvalue( eventTcl ) <- "OK" )      ## hit return = OK
        tkbind( tt, "<Destroy>", function() tclvalue( eventTcl ) <- "Cancel" )  ## destroy = cancel
        tkbind( tt, "<Escape>",  function() tclvalue( eventTcl ) <- "Cancel" )  ## hit escape = cancel
        tkbind( regList,  "<ButtonRelease-1>", function() { tclvalue( eventTcl ) <- "reg"  } )
        tkbind( termList, "<ButtonRelease-1>", function() { tclvalue( eventTcl ) <- "term" } )
        tkbind( selList,  "<ButtonRelease-1>", function() { tclvalue( eventTcl ) <- "sel"  } )

        # event loop
        repeat
        {
            tkwait.variable( eventTcl )
            eventVal <- tclvalue( eventTcl )
            if ( eventVal == "OK" )
            {
                # check frequency range
                opt$freqLo <- as.numeric( tclvalue( freqLoTcl ) )
                opt$freqHi <- as.numeric( tclvalue( freqHiTcl ) )
                if ( opt$freqLo < 0 || opt$freqLo >= opt$freqHi || opt$freqHi > 1 )
                {
                    tkmessageBox(
                        message = "Invalid range for frequency",
                        type = "ok",
                        icon = "error",
                        title = titan.version,
                        parent=tt
                    )
                    next
                }
                # check value of alpha
                alpha <- as.numeric( tclvalue( alphaTcl ) )
                if ( ( opt$reg != regOpt[ 3 ] && opt$sel == selOpt[ 1 ] )
                    && ( opt$alpha <= 0 || opt$alpha > 1 ) )
                    {
                    tkmessageBox(
                        message = "Invalid value for alpha",
                        type = "ok",
                        icon = "error",
                        title = titan.version,
                        parent = tt
                    )
                    next
                }
                tkdestroy( tt )
                if ( opt$reg == regOpt[ 3 ] )
                {
                    opt$sel <- selOpt[ 2 ]
                }
                break
            }
            if ( eventVal == "Cancel" )
            {
                tkdestroy( tt )
                return( NULL )
            }
            if ( eventVal == "reg" )
            {
                temp <- regOpt[ as.numeric( tkcurselection( regList ) ) + 1 ]
                if ( opt$reg != temp )
                {
                    reg <- match( temp, regOpt )
                    if ( temp == regOpt[ 3 ] || opt$reg == regOpt[ 3 ] )
                    {
                        opt$term <- termOpt[[ reg ]][ 1 ]
                    }
                    opt$reg <- temp
                    if ( opt$reg == regOpt[ 3 ] ) {
                        tkgrid(
                            tklabel( tt, text = "          " ),
                            row = 11,
                            column = 5,
                            sticky = "w"
                        )
                        tkgrid(
                            tklabel( tt, text = "                    " ),
                            row = 12,
                            column = 5,
                            sticky = "w"
                        )
                    }
                    else {
                        alphaEntry <- tkentry(
                            tt,
                            width = "8",
                            relief = "sunken",
                            textvariable = alphaTcl,
                            background = "white"
                        )
                        tkgrid(
                        tklabel( tt, text = "Alpha:" ),
                            row = 11,
                            column = 5,
                            sticky = "w"
                        )
                        tkgrid(
                            alphaEntry,
                            row = 12,
                            column = 5,
                            sticky = "w"
                        )
                    }
                    opt$sel <- selOpt[[ reg ]][ 1 ]
                }
            }
            if ( eventVal == "term" )
            {
                opt$term <- termOpt[[ reg ]][ as.numeric( tkcurselection( termList ) ) + 1 ]
            }
            if ( eventVal == "sel" )
            {
                temp <- selOpt[[ reg ]][ as.numeric( tkcurselection( selList ) ) + 1 ]
                if ( opt$sel != temp )
                {
                    opt$sel <- temp
                    if ( opt$sel == "wald" )
                    {
                        alphaEntry <- tkentry(
                            tt,
                            width = "8",
                            relief = "sunken",
                            textvariable = alphaTcl,
                            background = "white"
                        )
                        tkgrid(
                            tklabel( tt, text = "Alpha:" ),
                            row = 11,
                            column = 5,
                            sticky = "w"
                        )
                        tkgrid(
                            alphaEntry,
                            row = 12,
                            column = 5,
                            sticky = "w"
                        )
                    }
                    else
                    {
                        tkgrid(
                            tklabel( tt, text = "          " ),
                            row = 11,
                            column = 5,
                            sticky = "w"
                        )
                        tkgrid(
                            tklabel( tt, text = "                    " ),
                            row = 12,
                            column = 5,
                            sticky = "w"
                        )
                    }
                }
            }
            if ( eventVal == "reg" || eventVal == "term" || eventVal == "sel" )
            {
                tclvalue( regOptions  ) <- up( regStr, reg )
                tclvalue( termOptions ) <- up( termStr[[ reg ]], match( opt$term, termOpt[[ reg ]] ) )
                tclvalue( selOptions  ) <- up( selStr[[ reg ]], match( opt$sel, selOpt[[ reg ]] ) )
            }
        }
        opt$dataFile <- tclvalue( dataFile )
        opt$outFile <- tclvalue( outFile )
        opt$pdfFile <- tclvalue( pdfFile )
        opt$flagRaw <- as.logical( as.numeric( tclvalue( flagRawTcl ) ) )
        opt$flagFitted <- as.logical( as.numeric( tclvalue( flagFittedTcl ) ) )
        return( opt )
    }


    # Second input widget
    widget2 <-
    function( opt, gene, rx )
    {
        star <-
        function( x, n )
        {
            x <- paste( " ", x )
            x[ n + 1 ] <- paste( "*", sub( "  ", "", x[ n + 1 ] ) )
            return( x )
        }

        # initialise variables
        if ( is.null( opt$rx0 ) )
        {
            opt$rx0 <- rx[ 1 ]
        }
        if ( is.null( opt$rx1 ) )
        {
            opt$rx1 <- rx
        }
        if ( is.null( opt$gene0 ) )
        {
            opt$gene0 <- character( 0 )
            if ( is.null( opt$gene1 ) )
            {
                opt$gene1 <- gene
            }
        }
        else if ( is.null( opt$gene1 ) )
        {
            opt$gene1 <- gene[ -match( opt$gene0, gene ) ]
        }
        if ( any( is.na( match( opt$rx0, rx ) ) ) )
        {
            stop( "Unknown treatment in rx0", call. = FALSE )
        }
        if ( any( is.na( match( opt$rx1, rx ) ) ) )
        {
            stop( "Unknown treatment in rx1", call. = FALSE )
        }
        if ( any( is.na( match( opt$gene0, gene ) ) ) )
        {
            stop( "Unknown gene in gene0", call. = FALSE )
        }
        if ( any( is.na( match( opt$gene1, gene ) ) ) )
        {
            stop( "Unknown gene in gene1", call. = FALSE )
        }
        rx0 <- match( opt$rx0, rx ) - 1
        nrx <- length( rx )
        ngene <- length( gene )
        rxNames <- logical( nrx )
        rxNames[ match( opt$rx1, rx ) ] <- TRUE
        names( rxNames ) <- rx
        if ( nrx > 1 )
        {
            geneNames <- matrix( FALSE, nr = ngene, nc = 2 )
            geneNames[ match( opt$gene0, gene ), 2 ] <- TRUE
            colnames( geneNames ) <- c( "T", "H" )
        }
        else
        {
            geneNames <- matrix( FALSE, nr = ngene, nc = 1 )
            colnames( geneNames ) <- c( "Genes" )
        }
        geneNames[ match( opt$gene1, gene ), 1 ] <- TRUE
        rownames( geneNames ) <- gene
        nc <- function( x ) ceiling( sqrt( x / 5 ) )
        nr <- function( x, nc ) ceiling( x / nc )
        nCols <- max( nc( nrx - 1), nc( ngene ) )
        if ( is.null( opt$ciConf ) )
        {
            opt$ciConf <- c( 0.99, 0.95, 0.90 )
            ciConf <- c( FALSE, TRUE, FALSE )
            ciConfNames <- c( "99%", "95%", "90%" )
        }
        else
        {
            ciConf <- rep( TRUE, length( opt$ciConf ) )
            ciConfNames <- paste( opt$ciConf * 100, "%", sep = "" )
        }
        ciTypeNames <- c( "norm", "basic", "perc", "bca" )    ## not "stud" since we don't know the variances
        names( ciConf ) <- ciConfNames
        if ( is.null( opt$ciType ) || opt$ciType == "all" )
        {
            ciType <- rep( TRUE, 4 )
        }
        else
        {
            ciType <- rep( FALSE, 4 )
            ciType[ match( opt$ciType, ciTypeNames ) ] <- TRUE
        }
        names( ciType ) <- ciTypeNames
        rTcl <- tclVar( opt$R )
        if ( is.null( opt$seed ) )
        {
            seedTcl <- tclVar( as.character( ( runif( 1 ) * 99999999 ) %/% 1 + 1 ) )
        }
        else
        {
            seedTcl <- tclVar( opt$seed )
        }
        eventTcl <- tclVar( "NULL" )

        # set up widget
        tt <- tktoplevel()
        tkwm.title( tt, titan.version )
        if ( nrx > 1 )
        {
            rx0Options <- tclVar( star( rx, rx0 ) )
            args <- list(
                tt,
                height = min( nrx, 5 ),
                width = 20,
                selectmode = "single",
                background = "white",
                activestyle = "dotbox",
                relief = "sunken",
                listvariable = rx0Options
            )
            rx0List <- do.call( "tklistbox", args )
            tkselection.set( rx0List, rx0 )
            tkblankline( tt )
            tkgrid(
                tkblank( tt ),
                tklabel( tt, text = "Baseline treatment:" ),
                tkblank( tt ),
                tkblank( tt ),
                rx0List,
                sticky = "w"
            )
            if ( nrx > 5 )
            {
                tkgrid.configure( rx0List, column = 5, row = 1, rowspan = 5, sticky = "nsw" )
            }
            tkblankline( tt )
            rx1Cb <- checkButtons( tt, rxNames, "Comparison treatments:", nr( nrx - 1, nCols ) )
            tkconfigure( rx1Cb$cb[[ rx0 + 1 ]], foreground = "white" )
            tkconfigure( rx1Cb$label[[ rx0 + 1 ]], foreground = "gray" )
        }
        else
        {
            rx0 <- NULL
            rx1 <- rx
        }
        tkblankline( tt )
        if ( nrx > 1 )
        {
            gene1Cb <- checkButtons.nCols(
                tt,
                geneNames,
                "Select Test and \nHousekeeping genes:",
                nr( ngene, nCols )
            )
        }
        else
        {
            gene1Cb <- checkButtons.nCols(
                tt,
                geneNames,
                "Select genes:",
                nr( ngene, nCols )
            )
        }
        tkblankline( tt )
        ciConfCb <- checkButtons( tt, ciConf, "Confidence intervals:" )
        tkblankline( tt )
        ciTypeCb <- checkButtons( tt, ciType, "Confidence interval types:" )
        tkblankline( tt )
        rEnt <- tkentry( tt, width = "20", textvariable = rTcl )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Number of bootstrap replicates:" ),
            tkblank( tt ),
            tkblank( tt ),
            rEnt,
            sticky = "w"
        )
        tkblankline( tt )
        seedEnt <- tkentry( tt, width = "20", textvariable = seedTcl )
        tkgrid(
            tkblank( tt ),
            tklabel( tt, text = "Seed for random number generator:" ),
            tkblank( tt ),
            tkblank( tt ),
            seedEnt,
            sticky = "w"
        )
        okBut <- tkbutton(
            tt, text = "  OK  ",
            command = function() { tclvalue( eventTcl ) <- "OK" }
        )
        cancelBut <- tkbutton(
            tt,
            text = "  Back  ",
            command = function() { tclvalue( eventTcl ) <- "Cancel" }
        )
        tkblankline( tt )
        tkgrid(
            tkblank( tt ),
            tkblank( tt ),
            tkblank( tt ),
            okBut,
            cancelBut,
            sticky = "w"
        )
        tkblankline( tt )
        tkfocus( tt )

        # bind events
        tkbind( tt, "<Return>",  function() tclvalue( eventTcl ) <- "OK" )      ## hit return = OK
        tkbind( tt, "<Destroy>", function() tclvalue( eventTcl ) <- "Cancel" )  ## destroy = cancel
        tkbind( tt, "<Escape>",  function() tclvalue( eventTcl ) <- "Cancel" )  ## hit escape = cancel
        if ( nrx > 1 )
        {
            tkbind( rx0List, "<ButtonRelease-1>", function() tclvalue( eventTcl ) <- "rx0" )
        }
        for ( i in 1:ngene )
        {
            func <- eval( parse( text = paste( "function() tclvalue( eventTcl ) <- \"gT", i, "\"", sep = "" ) ) )
            args <- list( gene1Cb$cb[[ i ]][[ 1 ]], "<ButtonRelease-1>", func )
            do.call( "tkbind", args )
            if ( nrx > 1 )
            {
                func <- eval( parse( text = paste( "function() tclvalue( eventTcl ) <- \"gH", i, "\"", sep = "" ) ) )
                args <- list( gene1Cb$cb[[ i ]][[ 2 ]], "<ButtonRelease-1>", func )
                do.call( "tkbind", args )
            }
        }

        # event loop
        repeat
        {
            tkwait.variable( eventTcl )
            eventVal <- tclvalue( eventTcl )
            if ( eventVal == "OK" )
            {
                tkdestroy( tt )
                if ( nrx > 1 )
                {
                    opt$rx0 <- rx[ rx0 + 1 ]
                    rx1 <- checkButtonValues( rx1Cb, rx )
                    rx1[ rx0 + 1 ] <- FALSE
                    opt$rx1 <- rx[ rx1 ]
                }
                else
                {
                    opt$rx0 <- rx[ 1 ]
                    opt$rx1 <- NULL
                }
                opt$gene1 <- gene[ checkButtonValues( gene1Cb, gene, 1 ) ]
                if ( nrx > 1 )
                {
                    opt$gene0 <- gene[ checkButtonValues( gene1Cb, gene, 2 ) ]
                }
                opt$ciConf <- as.numeric( sub( "%", "", ciConfNames ) )[ checkButtonValues( ciConfCb, ciConfNames ) ] / 100
                opt$ciType <- ciTypeNames[ checkButtonValues( ciTypeCb, ciTypeNames ) ]
                opt$R <- as.numeric( tclvalue( rTcl ) )
                opt$seed <- as.numeric( tclvalue( seedTcl ) )
                return( opt )
            }
            if ( eventVal == "Cancel" )
            {
                tkdestroy( tt )
                return( NULL )
            }
            if ( eventVal == "rx0" && nrx > 1 )
            {
                tkconfigure( rx1Cb$cb[[ rx0 + 1 ]], foreground = "black" )
                tkconfigure( rx1Cb$label[[ rx0 + 1 ]], foreground = "black" )
                rx0 <- as.numeric( tkcurselection( rx0List ) )
                tkconfigure( rx1Cb$cb[[ rx0 + 1 ]], foreground = "white" )
                tkconfigure( rx1Cb$label[[ rx0 + 1 ]], foreground = "gray" )
                tclvalue( rx0Options ) <- star( rx, rx0 )
                }
            if ( length( grep( "gT", eventVal ) ) > 0 )
            {
                i <- as.numeric( sub( "gT", "", eventVal ) )
                if ( tclvalue( gene1Cb$value[[ i ]][[ 1 ]] ) == "1" )
                {
                    tclvalue( gene1Cb$value[[ i ]][[ 2 ]] ) <- "0"
                }
            }
            if ( length( grep( "gH", eventVal ) ) > 0 )
            {
                i <- as.numeric( sub( "gH", "", eventVal ) )
                if ( tclvalue( gene1Cb$value[[ i ]][[ 2 ]] ) == "1" )
                {
                    tclvalue( gene1Cb$value[[ i ]][[ 1 ]] ) <- "0"
                }
            }
        }
    }

    if ( widget == 1 )
    {
        widget1( ... )
    }
    else if (widget == 2 )
    {
        widget2( ... )
    }

}


# Widget input sequence
titan.input <-
function(
    data =       NULL,
    trace =      TRUE,
    widget =     TRUE,
    dataFile =   "",
    outFile =    "",
    pdfFile =    "",
    flagRaw =    FALSE,
    flagFitted = FALSE,
    freqLo =     .05,
    freqHi =     .95,
    reg =        "least.squares",
    term =       "linear",
    sel =        "wald",
    alpha =      .05,
    rx0 =        NULL,
    rx1 =        NULL,
    gene0 =      NULL,
    gene1 =      NULL,
    R =          1000,
    seed =       0,
    ciConf =     .95,
    ciType =     "all"
)
{
    get.data <-
    function( opt, data )
    {
        # get data
        if ( is.null( data ) )
        {
            dat <- try( titan.load.file( opt$dataFile ), silent = TRUE )
            if ( inherits( dat, "try-error" ) )
            {
                stop( paste( trim.error( dat ), "in MS data file" ),
                    call. = FALSE )
            }
        }
        else
        {
            dat <- data
        }
        datNames <- c( "FREQUENCY", "GENE", "TREATMENT", "RX", "COMPETITOR", "CONCENTRATION", "FLAG" )
        names( dat ) <- toupper( sub( ".*?(\\w+).*", "\\1", names( dat ), perl = TRUE ) )
        datCols <- match( 1:7, pmatch( names( dat ), datNames ) )
        if ( any( is.na( datCols[ 1:2 ] ) ) )
        {
            stop( "MS data does not contain column ", paste( datNames[ 1:2 ][ is.na( datCols[ 1:2 ] ) ], collapse = " or " ), call. = FALSE )
        }
        if ( all( is.na( datCols[ 3:4 ] ) ) )
        {
            stop( "MS data does not contain column ", "TREATMENT", call. = FALSE )
        }
        datCols[ 3 ] <- datCols[ 3:4 ][ !is.na( datCols[ 3:4 ] ) ][ 1 ]
        datCols <- datCols[ -4 ]
        if ( all( is.na( datCols[ 4:5 ] ) ) )
        {
            stop( "MS data does not contain column ", "COMPETITOR CONCENTRATION", call. = FALSE )
        }
        datCols[ 4 ] <- datCols[ 4:5 ][ !is.na( datCols[ 4:5 ] ) ][ 1 ]
        datCols <- datCols[ -5 ]
        if ( is.na( datCols[ 5 ] ) )
        {
            dat <- data.frame( dat[ , datCols[ 1:4 ], drop = FALSE ], FLAG = 0 )
        }
        else
        {
            dat <- data.frame( dat[ , datCols, drop = FALSE ] )
        }
        names( dat ) <- c( "freq", "gene", "rx", "comp.conc", "flag" )

        # factor variables
        dat$gene <- factor( dat$gene )
        dat$rx <- factor( dat$rx )

        # dependent variable = log ratio of test signal to competitor signal
        dat$y <- log10it( dat$freq )

        # independent variable = log of competitor concentration
        dat$x <- log10( dat$comp.conc )

        # flag extreme data points
        dat$flag[ ( dat$freq <= 0 ) | ( dat$freq >= 1 ) ] <- 1

        if ( sum( dat$flag == 0 ) == 0 )
        {
            stop( "No valid data points", call. = FALSE )
        }

    return( dat )
    }

    data.levels <-
    function( dat )
    {
        # create dummy variable matrix from factor
        dummy <- function( f )
        {
            len <- length( f )
            levels <- levels( f )
            nlevels <- nlevels( f )
            dummy <- matrix( nr = len, nc = nlevels )
            colnames( dummy ) <- levels
            for ( i in 1:nlevels )
                dummy[ , i ] <- as.numeric( f == levels[ i ] )
            return( dummy )
        }

        generx <- interaction( dat$rx, dat$gene, drop = TRUE )
        lev.generx <- levels( interaction( dat$rx, dat$gene[ order( as.character( dat$gene ) ) ] ) )
        dat$generx <- factor( lev.generx[ match( generx, lev.generx ) ], levels = lev.generx )

        # create dummy variables for factors
        dat$genei <- dummy( dat$gene )
        dat$rxi <- dummy( dat$rx )
        dat$generxi <- dummy( dat$generx )

        # rx contrasts
        if ( nlevels( dat$rx ) > 1 )
        {
            contr.rx <- contr.treatment( levels( dat$rx ) )
            ncontr.rx <- dim( contr.rx )[ 2 ]
            dat$contr.rxi <- matrix( nr = length( dat$x ), nc = ncontr.rx )
            for ( j in 1:ncontr.rx )
            {
                dat$contr.rxi[ , j ] <- ( dat$rxi %*% contr.rx[ , j ] )
            }
        }
        else
        {
            contr.rx <- NULL
            ncontr.rx <- 0
            dat$contr.rxi <- NULL
        }
        return(
            list(
                dat = dat,
                levels = list(
                    gene = levels( dat$gene ),
                    rx = levels( dat$rx ),
                    generx = levels( dat$generx ),
                    ngene = nlevels( dat$gene ),
                    nrx = nlevels( dat$rx ),
                    ngenerx = nlevels( dat$generx ),
                    ncontr.rx = ncontr.rx
                )
            )
        )
    }

    flipAndFlag <-
    function( data, opt, trace )
    {
        # do a quick regression to flip data which are 'upside down'
        # i.e. wrong frequency captured from sequenom file
        # robust regression, parallel linear, alpha = .05

        # first get rid of 0's and 1's ( y == -Inf and y == Inf )
        data$dat$flag[ ( data$dat$freq <= 0 ) | ( data$dat$freq >= 1 ) ] <- 1

        opt.check <- opt
        opt.check$reg <- "robust"
        opt.check$term <- "parallel.linear"
        opt.check$sel <- "wald"
        opt.check$alpha <- 0
        input.check <- list( data = data, opt = opt.check )
        cf.check <- titan.create.formulae( input.check )
        if ( length( cf.check$gene.valid ) == 0 ) {
            stop( "No genes have valid data", call. = FALSE )
        }
        reg.check <- titan.regression( input.check, cf.check, FALSE )
        for ( gene in seq( data$levels$ngene )[ cf.check$gene.valid ] )
        {
            if ( titan.find.coefs( reg.check, gene, 1, data$levels )$L > 0 )
            {
                data$dat$y[ data$dat$gene == data$levels$gene[ gene ] ] <- - data$dat$y[ data$dat$gene == data$levels$gene[ gene ] ]
            }
        }
        if ( opt$flagFitted )
        {
        # iteratively flag extreme data points and refit regression
        # to ensure that freqLo <= fitted freq <= freqHi

        # get highest and lowest competitor concentrations
        # and flag the most extreme fitted values ( if necessary )
            fittedMax <- log10it( opt$freqHi )
            fittedMin <- log10it( opt$freqLo )
            if ( trace )
            {
                cat( "\nRemoving data points outside valid range of fitted values" )
            }
            repeat
            {
                if ( trace )
                {
                    cat( "." )
                }
                input <- list( data = data, opt = opt )
                cf <- titan.create.formulae( list( data = data, opt = opt ) )
                if ( length( cf$gene.valid ) == 0 || all( cf$generx.valid == 0 ) )
                {
                    break
                }
                reg <- titan.regression( input, cf, FALSE )
                newFlag <- data$dat$flag
                moreFlagged <- 0
                # take fitted range FOR EACH GENE / TREATMENT
                for ( g in seq( data$levels$ngene ) )
                {
                    for ( r in seq( data$levels$nrx ) )
                    {
                        if ( !cf$generx.valid[ g, r ] )
                        {
                            next
                        }
                        notFlagged <- reg$call$subset
                        sel <- ( data$dat$gene[ notFlagged ] == data$levels$gene[ g ] ) & ( data$dat$rx[ notFlagged ] == data$levels$rx[ r ] )
                        fitted.range <- range( reg$fitted[ sel ] )
                        ##cat( "Fitted range: ", fitted.range[ 1 ], ", ", fitted.range[ 2 ], "\n", sep = "" )
                        if ( fitted.range[ 1 ] >= fittedMin && fitted.range[ 2 ] <= fittedMax )
                        {
                            next
                        }
                        extreme <- which.max( abs( fitted.range ) )
                        if ( extreme == 1 && fitted.range[ 1 ] < fittedMin )
                        {
                            newFlag[ notFlagged ][ sel ][ reg$fitted[ sel ] <= fitted.range[ 1 ] + .Machine$double.eps * 5 ] <- 1
                            moreFlagged <- 1
                            next
                        }
                        if ( extreme == 2 && fitted.range[ 2 ] > fittedMax )
                        {
                            newFlag[ notFlagged ][ sel ][ reg$fitted[ sel ] >= fitted.range[ 2 ] - .Machine$double.eps * 5 ] <- 1
                            moreFlagged <- 1
                            next
                        }
                        if ( fitted.range[ 1 ] < fittedMin )
                        {
                            newFlag[ notFlagged ][ sel ][ reg$fitted[ sel ] <= fitted.range[ 1 ] + .Machine$double.eps * 5 ] <- 1
                            moreFlagged <- 1
                            next
                        }
                        if ( fitted.range[ 2 ] > fittedMax )
                        {
                            newFlag[ notFlagged ][ sel ][ reg$fitted[ sel ] >= fitted.range[ 2 ] - .Machine$double.eps * 5 ] <- 1
                            moreFlagged <- 1
                            next
                        }
                    } # next rx
                } # next gene
                data$dat$flag <- newFlag
                if ( moreFlagged == 0 )
                {
                    break
                }
            }
        }

        if ( opt$flagRaw )
        {
            # flag extreme data points
            data$dat$flag[ ( data$dat$freq < opt$freqLo ) | ( data$dat$freq > opt$freqHi ) ] <- 1
        }

        if ( trace )
        {
            cat( "\n" )
        }
        if ( sum( data$dat$flag == 0 ) == 0 )
        {
            stop( "No valid data points: try a wider range of valid values", call. = FALSE )
        }

        return( data )
    }

    opt.init <- list(
        dataFile =   dataFile,
        outFile =    outFile,
        pdfFile =    pdfFile,
        flagRaw =    flagRaw,
        flagFitted = flagFitted,
        freqLo =     freqLo,
        freqHi =     freqHi,
        reg =        reg,
        term =       term,
        sel =        sel,
        alpha =      alpha,
        rx0 =        rx0,
        rx1 =        rx1,
        gene0 =      gene0,
        gene1 =      gene1,
        R =          R,
        seed =       seed,
        ciConf =     ciConf,
        ciType =     ciType
    )
    if ( widget && interactive() )
    {
        widget <- 1
        repeat
        {
            # input filenames
            if ( widget == 1 )
            {
                opt <- try( titan.widget( 1, opt.init, data ) )
                if ( is.null( opt ) )
                {
                    stop( "Input cancelled", call. = FALSE )
                }
                if ( inherits( opt, "try-error" ) )
                {
                    stop( "Input error", call. = FALSE )
                }
                opt.init <- opt
                widget <- 2
            }
            # select baseline rx
            # select bootstrap function and parameters
            # then preprocess data, run regression and find roots
            if ( widget == 2 )
            {
                dat <- try( get.data( opt, data ), silent = TRUE )
                if ( inherits( dat, "try-error" ) )
                {
                    err <- tkmessageBox(
                        message = dat,
                        type = "retrycancel",
                        default = "retry",
                        icon = "error",
                        title = titan.version
                    )
                    if ( tclvalue( err ) == "cancel" )
                    {
                        stop( trim.error( dat ), call. = FALSE )
                    }
                    widget <- 1
                    next
                }
                opt <- titan.widget( 2, opt.init, levels( dat$gene ), levels( dat$rx ) )
                if ( is.null( opt ) )
                {
                    widget <- 1
                }
                else
                {
                    break
                }
            }
        }
    }
    else
    {
        opt <- opt.init
    }

    if ( !is.null( data ) )
    {
        opt$dataFile <- ""
    }
    dat <- get.data( opt, data )
    if ( opt$reg == "spline" && is.numeric( opt$term ) )
    {
        opt$term <- paste( opt$term, "df", sep = "." )
    }
    rx <- levels( dat$rx )
    gene <- levels( dat$gene )
    if ( is.null( opt$rx0 ) )
    {
        opt$rx0 <- rx[ 1 ]
    }
    if ( is.null( opt$rx1 ) )
    {
        opt$rx1 <- rx[ -match( opt$rx0, rx ) ]
    }
    if ( is.null( opt$gene0 ) )
    {
        opt$gene0 <- character( 0 )
        if ( is.null( opt$gene1 ) )
        {
            opt$gene1 <- gene
        }
    }
    else if ( is.null( opt$gene1 ) )
    {
        opt$gene1 <- gene[ -match( opt$gene0, gene ) ]
    }
    if ( any( is.na( match( opt$rx0, rx ) ) ) )
    {
        stop( "Unknown treatment in rx0", call. = FALSE )
    }
    if ( any( is.na( match( opt$rx1, rx ) ) ) )
    {
        stop( "Unknown treatment in rx1", call. = FALSE )
    }
    if ( any( is.na( match( opt$gene0, gene ) ) ) )
    {
        stop( "Unknown gene in gene0", call. = FALSE )
    }
    if ( any( is.na( match( opt$gene1, gene ) ) ) )
    {
        stop( "Unknown gene in gene1", call. = FALSE )
    }
    if ( length( opt$gene1 ) == 0 )
    {
        stop( "No test genes selected", call. = FALSE )
    }
    ##if ( length( opt$rx1 ) == 0 )
    ##{
    ##    stop( "No test treatments selected", call. = FALSE )
    ##}

    # preprocess data
    if ( nlevels( dat$rx ) > 1 )
    {
        dat$rx <- relevel( dat$rx, match( opt$rx0, levels( dat$rx ) ) )
    }
    return(
        list(
            data = titan.process.data(
                flipAndFlag(
                    titan.process.data( data.levels( dat ) ),
                    opt,
                    trace
                )
            ),
            opt = opt
        )
    )
}

# Load data from text file
titan.load.file <-
function( filename = NULL, ... )
{
    if ( is.null( filename ) )
    {
        stop( "no file name" )
    }
    dat <- try( read.table( file = filename, header = TRUE, sep = "\t", ... ), silent = TRUE )
    if ( inherits( dat, "try-error" ) )
    {
        stop( trim.error( dat ) )
    }
    names( dat ) <- toupper( names( dat ) )
    return( dat )
}

titan.process.data <-
function( data )
{
    notFlagged <- ( data$dat$flag == 0 )

    # mean of logx by gene / rx
    generx.mean <- numeric( data$levels$ngenerx )
    names( generx.mean ) <- data$levels$generx
    is.na( generx.mean ) <- TRUE
    for ( i in 1:data$levels$ngenerx )
    {
        generx.mean[ i ] <- mean( data$dat$x[ notFlagged & data$dat$generxi[ , data$levels$generx[ i ] ] ] )
    }
    data$generx.mean <- generx.mean

    # take orthogonal polynomials with degree 2
    # for the log concentrations centred within gene / rx groups
    # because of the different
    # ranges of concentrations in the gene / rx groups,
    # the polynomials are not orthogonal within these groups
    # ( but it still helps )

    poly <- poly( data$dat$x[ notFlagged ] - generx.mean[ as.character( data$dat$generx[ notFlagged ] ) ], degree = 2 )
    data$poly.coefs <- attr( poly, "coefs" )
    data$dat$polyx1 <- polyL( data$dat$x - generx.mean[ data$dat$generx ], data$poly.coefs )
    data$dat$polyx2 <- polyQ( data$dat$x - generx.mean[ data$dat$generx ], data$poly.coefs )

    return( data )
}





###########################################################################




# linear model formulae
# the idea is to nest linear and quadratic x
# and their interactions with rx terms
# within dummy-coded gene effects
titan.create.formulae <-
function( input )
{
    gene.from.generx <-
    function( generx, gene )
    {
        match = 0
        for ( i in seq( length( gene ) ) )
        {
            if ( length( grep( paste( "\\.", gene[ i ], sep = "" ), generx, perl = TRUE ) ) > 0 )
            {
                match = i
                break
            }
        }
        return( match )
    }

    rx1 <- match( input$opt$rx1, input$data$levels$rx )
    rx1 <- rx1[ !is.na( rx1 ) ]
    genes <- c(
        match( input$opt$gene0, input$data$levels$gene ),
        match( input$opt$gene1, input$data$levels$gene )
    )
    genes <- genes[ !is.na( genes ) ]
    nrx1 <- length( rx1 )
    ngenes <- length( genes )

    # NB at first: generx.valid = # of distinct data points by gene / rx
    #       later: generx.valid = flag for whether there is a valid regression line for each gene / rx
    generx.valid <- matrix( NA, nr = input$data$levels$ngene, nc = input$data$levels$ncontr.rx + 1 )
    valid <- ( input$data$dat$flag == 0 )
    for ( i in genes )
    {
        generx.valid[ i, 1 ] <- length(
            unique(
                input$data$dat$polyx1[ valid ][ input$data$dat$genei[ valid, i ] & input$data$dat$rxi[ valid, 1 ] ]
            )
        )
        for ( j in rx1 )
        {
            generx.valid[ i, j ] <- length(
                unique(
                    input$data$dat$polyx1[ valid ][ input$data$dat$genei[ valid, i ] & input$data$dat$contr.rxi[ valid, j - 1 ] ]
                )
            )
        }
    }

    lower <- "y ~ -1 +"

    # Natural Spline
    if ( input$opt$reg == "spline" ) {
        regdf <- as.numeric( gsub( "\\D", "", input$opt$term, perl = TRUE ) )

        # Convert generx.valid to logical
        generx.valid[ is.na( generx.valid ) ] <- 0
        generx.valid[ generx.valid < regdf ]  <- 0
        generx.valid[ generx.valid >= regdf ] <- 1
        gene.valid <- rowSums( generx.valid )
        gene.valid <- seq( input$data$levels$ngene )[ as.logical( gene.valid ) ]

        generxi.valid <- numeric( input$data$levels$ngenerx )
        for ( i in seq( input$data$levels$ngenerx ) )
        {
            generxi.valid[ i ] <- length(
                unique(
                    input$data$dat$polyx1[ valid ][ as.logical( input$data$dat$generxi[ valid, i ] ) ]
                )
            ) >= regdf
            generxi.valid[ i ] <- generxi.valid[ i ] && gene.from.generx( input$data$levels$generx[ i ], input$data$levels$gene ) %in% gene.valid
        }
        generxi.valid <- seq( input$data$levels$ngenerx )[ as.logical( generxi.valid ) ]
        for ( i in generxi.valid )
        {
            lower <- paste(
                lower,
                "generxi[ ,", i, "] +"
            )
        }
        upper <- lower
        for ( i in generxi.valid )
        {
            upper <- paste(
                upper,
                "generxi[ ,", i, "]:ns( polyx1, df =", regdf, ") +"
            )
        }
        baseline.rx <- parallel.rx <- upper
    } # if spline

    else
    {
        gene.valid <- seq( input$data$levels$ngene )[ generx.valid[ , 1 ] > 1 ]
        gene.valid <- gene.valid[ !is.na( gene.valid ) ]
        for ( i in gene.valid )
        {
            lower <- paste(
                lower,
                "genei[ ,", i, "] +"
            )
        }
        baseline.rx <- lower
        for ( i in gene.valid )
        {
            baseline.rx <- paste(
                baseline.rx,
                "genei[ ,", i, "]:polyx1 +"
            )
        }
        parallel.rx <- baseline.rx
        for ( i in gene.valid )
        {
            if ( nrx1 > 0 )
            {
                for ( j in seq( input$data$levels$ncontr.rx ) )
                {
                    if ( !is.na( generx.valid[ i, j + 1 ] ) && generx.valid[ i, j + 1 ] > 1 )
                    {
                        lower <- paste(
                            lower,
                            "genei[ ,", i, "]:contr.rxi[ ,", j, "] +"
                        )
                        parallel.rx <- paste(
                            parallel.rx,
                            "genei[ ,", i, "]:contr.rxi[ ,", j, "] +"
                        )
                    }
                }
            }
            upper <- parallel.rx
        }

        # Parallel Linear
        if ( input$opt$term == "parallel.linear" && input$data$levels$ncontr.rx > 0 )
        {
            upper <- parallel.rx
        }

        # Linear
        else if ( input$opt$term == "linear" && input$data$levels$ncontr.rx > 0 )
        {
            upper <- baseline.rx
            for ( i in gene.valid )
            {
                if ( nrx1 > 0 )
                {
                    for ( j in seq( input$data$levels$ncontr.rx ) )
                    {
                        if ( !is.na( generx.valid[ i, j + 1 ] ) && generx.valid[ i, j + 1 ] > 1 )
                        {
                            upper <- paste(
                                upper,
                                "genei[ ,", i, "]:contr.rxi[ ,", j, "] +",
                                "genei[ ,", i, "]:polyx1:contr.rxi[ ,", j, "] +"
                            )
                        }
                    }
                }
            }
        }

        # Quadratic
        else if ( input$opt$term == "quadratic" )
        {
            upper <- parallel.rx <- baseline.rx
            for ( i in gene.valid )
            {
                if ( !is.na( generx.valid[ i, j + 1 ] ) && generx.valid[ i, j + 1 ] > 2 )
                {
                    baseline.rx <- paste(
                        baseline.rx,
                        "genei[ ,", i, "]:polyx2 +"
                    )
                    parallel.rx <- paste(
                        parallel.rx,
                        "genei[ ,", i, "]:polyx2 +"
                    )
                    upper <- paste(
                        upper,
                        "genei[ ,", i, "]:polyx2 +"
                    )
                }
                if ( nrx1 > 0 )
                {
                    for ( j in seq( input$data$levels$ncontr.rx ) )
                    {
                        if ( !is.na( generx.valid[ i, j + 1 ] ) && generx.valid[ i, j + 1 ] > 1 )
                        {
                            parallel.rx <- paste(
                                parallel.rx,
                                "genei[ ,", i, "]:contr.rxi[ ,", j, "] +"
                            )
                            upper <- paste(
                                upper,
                                "genei[ ,", i, "]:contr.rxi[ ,", j, "] +",
                                "genei[ ,", i, "]:polyx1:contr.rxi[ ,", j, "] +"
                            )
                            if ( generx.valid[ i, j + 1 ] > 2 )
                            {
                                upper <- paste(
                                    upper,
                                    "genei[ ,", i, "]:polyx2:contr.rxi[ ,", j, "] +"
                                )
                            }
                        }
                    }
                }
            }
        }

        # Convert generx.valid to logical
        if ( length( gene.valid ) == 0 )
        {
            generx.valid[ , ] <- 0
        }
        else
        {
            generx.valid[ -gene.valid, ]          <- 0
            generx.valid[ is.na( generx.valid ) ] <- 0
            generx.valid[ generx.valid <= 1 ]     <- 0
            generx.valid[ generx.valid > 1 ]      <- 1
        }

    } # not spline

    upper <- as.formula( substr( upper, 1, nchar( upper ) - 1 ) )
    baseline.rx <- as.formula( substr( baseline.rx, 1, nchar( baseline.rx ) - 1 ) )
    parallel.rx <- as.formula( substr( parallel.rx, 1, nchar( parallel.rx ) - 1 ) )
    lower <- as.formula( substr( lower, 1, nchar( lower ) - 1 ) )
    return(
        list(
            form = list(
                upper = upper,
                parallel.rx = parallel.rx,
                baseline.rx = baseline.rx,
                lower = lower
            ),
            generx.valid = generx.valid,
            gene.valid = gene.valid
        )
    )
}

# run regression
titan.regression <-
function( input, cf, trace = TRUE )
{
    dropWald <-
    function( reg, alpha, scope, trace )
    {
        if ( alpha == 0 )
        {
            return( reg )
        }
        r <- summary( reg, correlation = FALSE )
        repeat
        {
            ds <- drop.scope( reg$terms, scope )
            if ( length( ds ) == 0 )
            {
                break
            }
            coef <- coef( r )[ rownames( coef( r ) ) %in% ds, , drop = FALSE ]
            drop <- rownames( coef )[ which.min( abs( coef[ , 3 ] ) ) ]
            if ( pt( -abs( coef( r )[ rownames( coef( r ) ) == drop, 3 ] ), df = r$df[ 2 ] ) < alpha / 2 )
            {
                break
            }
            new.formula <- update.formula( formula( reg$call ), as.formula( paste( "~ . -", ds[ drop == ds ] ) ) )
            reg <- update( reg, formula = new.formula )
            r <- summary( reg, correlation = FALSE )
            if ( trace )
            {
                cat( "Drop term", titan.regsub( ds[ drop == ds ], input$data$levels ), "\n" )
                r2 <- r
                r2$call$formula <- titan.regformula( r2$call$formula, input$data$levels )
                rownames( r2$coefficients ) <- titan.regsub( rownames( coef( r2 ) ), input$data$levels )
                r2$call$data <- "data"
                r2$call$subseet <- "subset"
                ##print( r2 )
            }
        }
        return( reg )
    }

    if ( length( cf$gene.valid ) == 0 )
    {
        warning( "no genes have valid data" )
        return( NULL )
    }
    if ( trace )
    {
        cat( "Performing regression" )
    }
    method <- c( "lm", "rlm" )[ ( input$opt$reg == "robust" ) + 1 ]
    args <- list(
        formula = cf$form$upper,
        data = input$data$dat,
        subset = ( ( input$data$dat$flag == 0 ) & rowSums( input$data$dat$genei[ , cf$gene.valid, drop = FALSE ] ) )
    )
    if ( input$opt$reg == "robust" )
    {
        args <- c( args, method = "M", psi = "psi.huber" )
    }
    reg <- try( do.call( method, args ), silent = FALSE )
    if ( inherits( reg, "try-error" ) )
    {
        tkmessageBox(
            message = "Regression failed",
            type = "ok",
            icon = "error",
            title = titan.version
        )
        stop( sub( "Error *(.*)", "\\1", reg, perl = TRUE ), call. = FALSE )
    }
    if ( trace )
    {
        cat( "\n" )
        cat( c( "Least squares", "Robust" )[ ( input$opt$reg == "robust" ) + 1 ] )
        cat( " regression with " )
        if ( input$opt$reg != "spline" )
        {
            cat( sub( "\\.", " ", input$opt$term, perl = TRUE ) )
            cat( " terms\n" )
            cat(
                c(
                    paste( "Backwards selection by Wald test (alpha = ", input$opt$alpha, ")\n", sep = "" ),
                    "Stepwise selection by AIC\n"
                )[ ( input$opt$sel == "aic" ) + 1 ]
            )
        }
        else
        {
            cat( "natural spline terms (" )
            cat( sub( "\\.", " ", input$opt$term, perl = TRUE ) )
            cat( ")\nStepwise selection by AIC\n" )
        }
        r <- summary( reg, correlation = FALSE )
        r$call <- titan.regformula( r$call$formula, input$data$levels )
        rownames( r$coefficients ) <- titan.regsub( rownames( coef( r ) ), input$data$levels )
        print( r )
    }

    # Stepwise selection by AIC
    if ( input$opt$sel == "aic" ) ## && cf$form$upper != cf$form$parallel.rx )
    {
        reg <- stepAIC( reg, scope = list( lower = cf$form$parallel.rx, upper = cf$form$upper ), trace = 0 )
        ## reg <- stepAIC( reg, scope = list( lower = cf$form$lower ),
        ##     trace = 0 )
        reg$anova$Step <- titan.regsub( as.character( reg$anova$Step ), input$data$levels )
        if ( trace )
        {
            cat( "Stepwise Model Path\nAnalysis of Deviance Table\n\n" )
            print.data.frame( reg$anova )
            r <- summary( reg, correlation = FALSE )
            r$call <- titan.regformula( r$call$formula, input$data$levels )
            rownames( r$coefficients ) <- titan.regsub( rownames( coef( r ) ), input$data$levels )
            ## print( r )
        }
    }

    # Backwards selection by Wald test
    if ( input$opt$sel == "wald" ) ## && cf$form$upper != cf$form$parallel.rx )
    {
        ## reg <- dropWald( reg, input$opt$alpha, scope = cf$form$lower, trace )
        reg <- dropWald( reg, input$opt$alpha, scope = cf$form$parallel.rx, trace )
    }

    # Results per gene
    gene.res <- matrix( NA, 0, 4 )
    dat <- reg$call$data[ reg$call$subset, ]
    for ( gene in cf$gene.valid )
    {
        sel <- ( dat$gene == input$data$levels$gene[ gene ] )
        means.x.rx <- tapply( dat$y[ sel ], list( dat$x[ sel ], dat$rx[ sel ] ), mean )
        means.row <- match( dat$x[ sel ], rownames( means.x.rx ) )
        means.col <- match( dat$rx[ sel ], colnames( means.x.rx ) )
        means <- as.vector( means.x.rx )[ means.row + ( means.col - 1 ) * dim( means.x.rx )[ 1 ] ]
        rss <- sum( reg$residuals[ sel ] ^ 2 )
        mss <- sum( reg$fitted[ sel ] ^ 2 )
        ss <- rss + mss
        rss.res <- sum( ( dat$y[ sel ] - means ) ^ 2 )
        rss.lof <- sum( ( means - reg$fitted[ sel ] ) ^ 2 )
        gene.res <- rbind( gene.res, c( ss, rss.lof, rss.res, mss ) )
    }
    # Overall results
    gene.res <- rbind( gene.res, colSums( gene.res ) )
    dimnames( gene.res ) <- list(
        c( input$data$levels$gene[ cf$gene.valid ], "Total" ),
        c( "Sum Sq", "Lack of Fit","Residuals","R-Squared" )
    )
    gene.res[ , 2:4 ] <- gene.res[ , 2:4 ] / gene.res[ , 1 ]

    # Results per treatment per gene
    generx.res <- array( NA, dim = c( length( cf$gene.valid ) + 1, length( input$opt$rx0 ) + length( input$opt$rx1 ), 4 ) )
    for ( g in seq( length( cf$gene.valid ) ) )
    {
        for ( r in seq( length( input$opt$rx0 ) + length( input$opt$rx1 ) ) )
        {
            sel <- (
                ( dat$gene == input$data$levels$gene[ cf$gene.valid[ g ] ] ) &
                ( dat$rx == c( input$opt$rx0, input$opt$rx1 )[ r ] )
            )
            rss <- sum( reg$residuals[ sel ] ^ 2 )
            mss <- sum( reg$fitted[ sel ] ^ 2 )
            ss <- rss + mss
            generx.res[ g, r, ] <- c( ss, rss, mss, NA )
        }
    }
    # Overall results
    generx.res[ length( cf$gene.valid ) + 1, , ] <- 0
    generx.res[ length( cf$gene.valid ) + 1, , ] <- apply( generx.res, 2:3, sum )
    dimnames( generx.res ) <- list(
        c( input$data$levels$gene[ cf$gene.valid ], "Total" ),
        c( input$opt$rx0, input$opt$rx1 ),
        c( "Sum Sq", "Resid SS", "Mean SS", "R-Squared" )
    )
    generx.res[ , , 4 ] <- generx.res[ , , 3 ] / generx.res[ , , 1 ]

    reg$gene.res <- gene.res
    reg$generx.res <- generx.res
    reg$levels$gene <- input$dat$levels$gene
    reg$levels$rx <- input$dat$levels$rx
    reg$levels$generx <- input$dat$levels$generx

    # print final model
    if ( trace )
    {
        cat( "\nFinal Model:\n" )
        titan.printreg( reg )
    }

    return( reg )
}




###########################################################################




titan.find.coefs <-
function( reg, gene, rx, levels )
{
    if ( any( is.na( coef( reg ) ) ) )
    {
        ##warning( "NA in regression coefficients", call. = FALSE )
        stop( "NA in regression coefficients" )
        return( list( I = NA, L = NA, Q = NA ) )
    }
    terms <- names( coef( reg ) )
    g <- gr <- gri <- grep( paste( "genei\\[, ", gene, "]", sep = "" ), terms )
    r <- grep( "contr.rxi", terms[ g ] )
    if ( length( r ) > 0 )
    {
        gr <- g[ -r ]                               # baseline rx
    }
    if ( rx > 1 )
    {
        gr2 <- grep( paste( "contr.rxi\\[, ", rx - 1, "]", sep = "" ), terms[ g ] )
        if ( length( gr2 ) > 0 )
        {
            gr <- c( gr, g[ gr2 ] )
        }
    }
    grx <- grep( "polyx", terms[ gr ] )
    if ( length( grx ) > 0 )
    {
        gri <- gr[ -grx ]
    }
    I <- coef( reg )[ gri ]
    L <- coef( reg )[ gr[ grep( "polyx1", terms[ gr ] ) ] ]
    Q <- coef( reg )[ gr[ grep( "polyx2", terms[ gr ] ) ] ]
    if ( length( L ) == 0 ) L <- 0
    if ( length( Q ) == 0 ) Q <- 0
    if ( length( I ) == 0 )
    {
        warning(
            "regression not possible for gene ",
            levels$gene[ gene ], ", rx ",
            levels$rx[ rx ],
            call. = FALSE
        )
    }
    else if ( all( ( L == 0 ) & ( Q == 0 ) ) )
    {
         warning(
            "no regression slope estimated for gene ",
             levels$gene[ gene ], ", rx ",
             levels$rx[ rx ],
             call. = FALSE
         )
    }
    return( list( I = sum( I ), L = sum( L ), Q = sum( Q ) ) )
}

# interpolate value for x at y
titan.interpolate.line <-
function( data, reg, y, gene, rx, trace = FALSE )
{
    pred.y <-
    function( data, reg, gene, rx, x, ... )
    {
        nx <- length( x )
        sel <- ( 1:dim( data$dat )[ 1 ] )[ data$dat$gene == data$levels$gene[ gene ] & data$dat$rx == data$levels$rx[ rx ] ][ 1 ]
        newdata <- data.frame(
            polyx1 = polyL( x, data$poly.coefs ),
            polyx2 = polyQ( x, data$poly.coefs )
        )
        newdata$genei <- data$dat$genei[ rep( sel, nx ), , drop = FALSE ]
        newdata$contr.rxi <- data$dat$contr.rx[ rep( sel, nx ), , drop = FALSE ]
        newdata$generxi <- data$dat$generxi[ rep( sel, nx ), , drop = FALSE ]
        return( predict( reg, newdata, ... ) )
    }

    interpolate.polynomial <-
    function( y, coefs, p, trace = FALSE )
    {
        # the solution of a quadratic curve
        # y == Q * d2 + L * d1 + I
        # is a quadratic of the form
        # a * x^2 + b * x + c == 0
        # where
        # a = Q / sqrt( p$norm2[4] )
        # b = - Q * ( p$alpha[1] + p$alpha[2] ) / sqrt( p$norm2[4] ) + L / sqrt( p$norm2[3] )
        # c = Q * ( p$alpha[1] * p$alpha[2] - p$norm2[3] / p$norm2[2] ) / sqrt( p$norm2[4] ) - L * p$alpha[1] / sqrt( p$norm2[3] ) + I - y
        # with roots
        # ( - b +- sqrt( b^2 - 4 * a * c ) ) / ( 2 * a )
        # which are real when b^2 >= 4 * a * c
        # only one root when b^2 == 4 * a * c

        ##if ( coefs$Q == 0 && coefs$L == 0 )
        ##    warning( "Regression slope is constant", call. = FALSE )
        a <- coefs$Q / sqrt( p$norm2[4] )
        b <- - coefs$Q * ( p$alpha[1] + p$alpha[2] ) / sqrt( p$norm2[4] ) + coefs$L / sqrt( p$norm2[3] )
        c <- coefs$Q * ( p$alpha[1] * p$alpha[2] - p$norm2[3] / p$norm2[2] ) / sqrt( p$norm2[4] ) - coefs$L * p$alpha[1] / sqrt( p$norm2[3] ) + coefs$I - y
        d <- b^2 - 4 * a * c
        if ( trace ) {
            cat( "\ta", a, "\n\tb", b, "\n\tc", c, "\n\td", d, "\n" )
        }
        if ( coefs$Q == 0 )
        {
            return( - c / b )
        }
        if ( d < 0 )
        {
            ##warning( "Regression slope never reaches y == ", y, call. = FALSE )
            return( NA )
        }
        ##if ( d == 0 )
        ##    warning( "Regression slope does not cross y == ", y, call. = FALSE )
        roots <- ( - b + ( c( 1, -1 ) * sqrt( d ) ) ) / ( 2 * a )
        if ( trace)
        {
            cat( "quadratic roots", roots, "\n" )
        }
        # always take the root where the slope is negative
        # this is the lower root if Q > 0
        # and the upper root if Q < 0
        return( ifelse( coefs$Q < 0, max( roots ), min( roots ) ) )
    }

    # find root of spline curve by bisection
    # cutting x.range into successive deciles
    # NB this method may not find the root
    # if the curve is not a monotonic decreasing function of x
    # and will not detect the presence of multiple roots
    interpolate.spline <-
    function( data, y, x.range, reg, gene, rx, decreasing = TRUE, maxiter = 10, epsilon = 1e-5, trace = FALSE )
    {
        if ( decreasing )
        {
            xr <- range( x.range )
        }
        else
        {
            xr <- range( x.range )[ 2:1 ]
        }
        xseq <- seq( xr[ 1 ], xr[ 2 ], length = 11 )
        yseq <- pred.y( data, reg, gene, rx, xseq )
        if ( all( diff( sign( yseq ) ) >= 0 ) )
        {
            # expand range once only
            xseq <- seq( xr[ 1 ] - 3 * diff( xr ), xr[ 2 ] + 3 * diff( xr ), length = 11 )
            yseq <- pred.y( data, reg, gene, rx, xseq )
            if ( all( diff( sign( yseq ) ) >= 0 ) )
            {
                warning(
                    "cannot find root for gene ",
                    data$levels$gene[ gene ],
                    ", rx ",
                    data$levels$rx[ rx ],
                    call. = FALSE
                )
                return( NA )
            }
        }
        i <- which.min( diff( sign( yseq ) ) )
        xr <- xseq[ i:( i + 1 ) ]
        conv <- FALSE
        for ( iter in 1:maxiter )
        {
            ##cat( "xr =", xr, "\n" )
            xseq <- seq( xr[ 1 ], xr[ 2 ], length = 11 )
            yseq <- pred.y( data, reg, gene, rx, xseq )
            ##for ( i in 1:11 )
            ##    cat( "\tx", xseq[ i ], "\ty", yseq[ i ], "\n" )
            i <- which.min( diff( sign( yseq ) ) )
            yp <- yseq[ i:( i + 1 ) ]
            xr <- xseq[ i:( i + 1 ) ]
            if ( abs( diff( yp ) ) / ( .1 + abs( mean( yp ) ) ) < epsilon )
            {
                conv <- TRUE
                break
            }
        }
        if ( !conv )
        {
            warning(
                "bisection did not converge for gene",
                data$levels$gene[ gene ],
                ", rx ",
                data$levels$rx[ rx ],
                call. = FALSE
            )
        }
        return( mean( xr ) )
    }

    x.offset <- data$generx.mean[ ( gene - 1 ) * data$levels$nrx + rx ]
    if ( length( grep( "\\:ns\\(",
        attr( reg$terms, "term.labels" ) ) ) > 0 )
        {
        r <- try(
            interpolate.spline(
                data,
                0,
                range( data$dat$x - x.offset ),
                reg,
                gene,
                rx,
                trace = trace
            ),
            silent = FALSE
        )
    }
    else {
        coefs <- titan.find.coefs( reg, gene, rx, data$levels )
        if ( trace )
        {
            cat( "\n\tI", coefs$I )
            cat( "\n\tL", coefs$L )
            cat( "\n\tQ", coefs$Q, "\n" )
        }
        if ( any( is.na( coefs ) ) || length( coefs$I ) == 0 )
        {
            return( NA )
        }
        r <- try(
            interpolate.polynomial(
                0,
                coefs = coefs,
                p = data$poly.coefs,
                trace = trace
            ),
            silent = FALSE
        )
        if ( is.na( r ) )
        {
            warning(
                "Regression slope never reaches y == ", y,
                " for gene ", data$levels$gene[ gene ], ", rx ",
                data$levels$rx[ rx ], "\n",
                call. = FALSE
            )
        }
    }
    if ( inherits( r, "try-error" ) )
    {
        r[ 1 ] <- paste(
            trim.error( r[ 1 ] ), "Error occurred for gene ",
            data$levels$gene[ gene ], ", rx ",
            data$levels$rx[ rx ], "\n",
            sep = ""
        )
        stop( r, call. = FALSE )
    }

    return( r + x.offset )
}

titan.interpolate <-
function( input, reg, gene.valid, generx.valid, trace = TRUE )
{
    rx1 <- match( input$opt$rx1, input$data$levels$rx )
    rx1 <- rx1[ !is.na( rx1 ) ]
    nrx1 <- length( rx1 )
    roots <- matrix( NA, nr = length( gene.valid ),
        nc = 1 + nrx1 )
    dimnames( roots ) <- list( input$data$levels$gene[ gene.valid ], input$data$levels$rx[ c( 1, rx1 ) ] )
    for ( i in seq( length( gene.valid ) ) )
    {
        for( j in 0:nrx1 )
        {
            rx <- 1
            if ( j > 0 )
            {
                rx <- rx1[ j ]
            }
            root <- try(
                titan.interpolate.line( input$data, reg, 0, gene.valid[ i ], rx, trace = FALSE ),
                silent = FALSE
            )
            if  ( !inherits( root, "try-error" ) )
            {
                roots[ i, j + 1 ] <- root
            }
        }
    }
    log10fold <- NULL
    gene0 <- match( input$opt$gene0, input$data$levels$gene[ gene.valid ] )
    gene1 <- match( input$opt$gene1, input$data$levels$gene[ gene.valid ] )
    gene0 <- gene0[ !is.na( gene0 ) ]
    gene1 <- gene1[ !is.na( gene1 ) ]
    ngene0 <- length( gene0 )
    ngene1 <- length( gene1 )
    if ( nrx1 > 0 )
    {
        if ( ngene1 )
        {
            log10fold <- matrix( NA, nr = ngene1, nc = nrx1 )
            dimnames( log10fold ) <- list(
                input$data$levels$gene[ gene.valid ][ gene1 ],
                input$data$levels$rx[ rx1 ]
            )
            if ( ngene0 )
            {
                log10fold.g0 <- matrix( NA, nr = ngene0, nc = nrx1 )
                for ( g in seq( ngene0 ) )
                {
                    for ( r in seq( nrx1 ) )
                    {
                        log10fold.g0[ g, r ] <- ( roots[ gene0[ g ], rx1[ r ] ] - roots[ gene0[ g ], 1 ] )
                    }
                }
            }
            for ( g in seq( ngene1 ) )
            {
                for ( r in seq( nrx1 ) )
                {
                    log10fold[ g, r ] <- ( roots[ gene1[ g ], rx1[ r ] ] - roots[ gene1[ g ], 1 ] )
                }
            }
            if ( ngene0 )
            {
                log10fold <- log10fold - matrix( colMeans( log10fold.g0 ), nr = ngene1, nc = nrx1, byrow = TRUE )
            }
        }
    }
    interpolation <- list(
        roots = roots,
        log10fold = log10fold,
        rx0 = input$opt$rx0,
        rx1 = input$data$levels$rx[ rx1 ],
        gene0 = input$data$levels$gene[ gene.valid ][ gene0 ]
    )
    if ( trace )
    {
        titan.print.interpolation( interpolation )
    }
    if ( length( gene0 ) != length( input$opt$gene0 ) )
    {
        cat( "Warning: invalid data for one or more housekeeping genes\n" )
    }
    return( interpolation )
}





###########################################################################





# graph the data

titan.plot <-
function( input, reg, cf, trace = TRUE, ... )
{
    titan.plot.internal <-
    function( input, reg, cf, cex = 0.8, ok = FALSE, ... )
    {
        for ( gene in cf$gene.valid )
        {
            gsel <- titan.plot.gene(
                gene, input$data,
                reg,
                rx.valid = cf$generx.valid[ gene, ],
                input$opt$freqLo,
                input$opt$freqHi,
                cex,
                ...
            )
            if ( ok && interactive() && gene != cf$gene.valid[ length( cf$gene.valid ) ] )
            {
                msg <- tkmessageBox(
                    message = "Review the next gene",
                    type = "ok",
                    icon = "info",
                    title = titan.version
                )
            }
        }
    }

    if ( input$opt$pdfFile != "" )
    {
        pdf( file = input$opt$pdfFile, height = 6, width = 8, onefile = TRUE )
        titan.plot.internal( input, reg, cf, cex = 0.8, ok = FALSE, ... )
        dev.off()
    }
    else if ( trace )
    {
        cols <- ceiling( sqrt( length( cf$gene.valid ) ) )
        rows <- ceiling( length( cf$gene.valid ) / cols )
        par( mfrow = c( rows, cols ) )
        titan.plot.internal( input, reg, cf, cex = 0.8, ok = TRUE, ... )
    }
}

titan.plot.gene <-
function( gene, data, reg, rx.valid, freqLo, freqHi, cex, ... )
{
    titan.plot.generx <-
    function( gene, rx, data, reg, xr, colour, se = TRUE, trace = FALSE )
    {
        len.p <- 2 * dim( data$dat )[ 2 ]
        sel <- ( data$dat$gene == data$levels$gene[ gene ] & data$dat$rx == data$levels$rx[ rx ] )
        x <- data$dat$x[ sel ] - data$generx.mean[ data$dat$generx[ sel ] ]
        new.x <- seq( min( x ), max( x ), length = len.p )
        newdata <- data.frame(
            polyx1 = polyL( new.x, data$poly.coefs ),
            polyx2 = polyQ( new.x, data$poly.coefs )
        )
        sel1 <- which.max( sel )
        newdata$genei <- data$dat$genei[ rep( sel1, len.p ), , drop = FALSE ]
        newdata$contr.rxi <- data$dat$contr.rxi[ rep( sel1, len.p ), , drop = FALSE ]
        newdata$generxi <- data$dat$generxi[ rep( sel1, len.p ), , drop = FALSE ]
        reg.p <- predict( reg, newdata, se.fit = se )
        new.x <- new.x + data$generx.mean[ rep( data$dat$generx[ sel1 ], len.p ) ]
        panel.lines( new.x, reg.p$fit, col = colour )
        if ( inherits( reg, "rlm" ) )
        {
            df.res <- summary( reg )$df[ 2 ]
        }
        else
        {
            df.res <- reg$df.res
        }
        k <- qt( .975, df = df.res )
        if ( se )
        {
            panel.lines( new.x, reg.p$fit + k * reg.p$se.fit, col = colour, lty = 3 )
            panel.lines( new.x, reg.p$fit - k * reg.p$se.fit, col = colour, lty = 3 )
        }
        p <- c( 4, 97:119, 121:122, 65:87, 88:89 )[ as.numeric( data$dat$sample[sel] ) * ( data$dat$flag[sel] == 0 ) + 1 ]
        panel.points(
            data$dat$x[ sel ],
            data$dat$y[ sel ],
            col = colour,
            pch = 4 - 3 * ( data$dat$flag[ sel ] == 0 )
        )
    }

    len <- dim( data$dat )[ 2 ]
    gsel <- ( data$dat$gene == data$levels$gene[ gene ] )
    xr <- range( data$dat$x[ gsel ] )
    yr <- range( data$dat$y[ gsel ][ is.finite( data$dat$y[ gsel ] ) ] )
    yr <- mean(yr) + diff(yr) * 0.5 * ( 1.04 * c( -1, 1 ) )
    xl <- xr[ 1 ]
    yl <- c( 4, 1 ) %*% yr / 5
    par( mar = c( 5, 4, 4, 5 ) + 0.1 )
    rx.valid <- as.logical( rx.valid )
    cols <- numeric( data$levels$nrx )
    cols[ rx.valid ] <- seq( sum( rx.valid ) )
    if ( sum( rx.valid ) )
    {
        t <- paste( "                         ", c( data$levels$rx[ rx.valid ], "flagged spot" ) )
        plot.new()
        titan.xyplot <- xyplot(
            dat$y ~ dat$x,
            data = data,
            subset = gsel,
            xlab = list( "log10 competitor concentration", cex = cex ),
            ylab = list( "log10 ( test signal / competitor signal )", cex = cex ),
            ylim = yr,
            panel = function()
            {
                for ( rx in seq( data$levels$nrx )[ rx.valid ] )
                {
                    titan.plot.generx( gene, rx, data, reg, xr, cols[ rx ] )
                }
                panel.abline( h = 0, col = 8, lty = 1 )
                panel.abline( h = log10it( freqLo ), col = 8, lty = 2 )
                panel.abline( h = log10it( freqHi ), col = 8, lty = 2 )
            },
            auto.key = TRUE,
            key = list(
                text =  list(
                    t,
                    cex = cex
                ),
                lines = list(
                    col = c( seq( sum( rx.valid ) ), 1 ),
                    pch = c( rep( 1, sum( rx.valid ) ), 4 ),
                    lty = rep( 1, sum( rx.valid ) + 1 ),
                    type = "o",
                    size = 3
                ),
                space = "right",
                divide = 1
            ),
            main = data$levels$gene[ gene ],
            scales = list(
                y = list( col = 0, tck = 0 ),
                ...
            )
        )
        plot( titan.xyplot )
        trellis.focus( "panel", column = 1, row = 1, clip.off = TRUE, highlight = FALSE )
        panel.axis( "left", outside = TRUE, check.overlap = TRUE )
        y2ticks <- c( -1, .01, .05, seq( .1, .9, by = .1 ), .95, .99 )
        y2seq <- y2ticks[
            ( which.max( y2ticks[ exp10it( yr[ 1 ] ) > y2ticks ] ) + 1 ):
            ( which.max( y2ticks[ exp10it( yr[ 2 ] ) >= y2ticks ] ) ) ]
        xlim <- current.panel.limits()$xlim
        ylim <- current.panel.limits()$ylim
        panel.axis( "right", outside = TRUE, at = ( log10it(y2seq) - yr[1] ) / diff(yr) * diff(ylim) + ylim[1], labels = y2seq, check.overlap = TRUE )
        panel.text( xr[2] + diff(xr) * 0.2, mean(yr), "Frequency = test signal / ( test signal + competitor signal )", cex = cex, srt = -90 )
        trellis.unfocus()
    }
    return( gsel )
}

titan.manualFlag <-
function( input, reg, cf, ... )
{
    if ( !interactive() )
    {
        return( input$data$dat$flag )
    }
    ok <- tkmessageBox(
        message = paste(
            "Left click mouse\tto change the flag\n\n",
            "Right click mouse\nand <Stop>\tto review the next gene",
            sep = ""
        ),
        type = "ok",
        icon = "info",
        title = titan.version
    )
    par( mfrow = c( 1, 1 ) )
    for ( gene in cf$gene.valid )
    {
        gsel <- titan.plot.gene(
            gene,
            input$data,
            reg,
            rx.valid = cf$generx.valid[ gene, ],
            input$opt$freqLo,
            input$opt$freqHi,
            cex = 0.8,
            ...
        )
        trellis.focus( "panel", column = 1, row = 1, clip.off = TRUE, highlight = FALSE )
        i <- which( gsel )[
            panel.identify(
                input$data$dat$x[ gsel ],
                input$data$dat$y[ gsel ],
                labels = c( "X","O" )[ input$data$dat$flag[ gsel ] + 1 ],
                offset = -.3,
                threshold = 18
            )
        ]
        input$data$dat$flag[ i ] <- 1 - input$data$dat$flag[ i ]
    }
    if (
        input$opt$dataFile != "" &&
        tclvalue(
            tkmessageBox(
                message = "Save flags?",
                type = "yesno",
                default = "yes",
                icon = "question",
                title = titan.version
            )
        )  == "yes"
    )
    {
        data <- try( titan.load.file( input$opt$dataFile ), silent = TRUE )
        if ( inherits( data, "try-error" ) )
        {
            stop( paste( trim.error( data ), "in MS data file" ), call. = FALSE )
        }
        if ( all( is.na( pmatch( names( data ), "FLAG" ) ) ) )
        {
            data$FLAG <- NA
        }
        flagCol <- which.max( pmatch( names( data ), "FLAG" ) )
        data[ , flagCol ] <- input$data$dat$flag
        write.table( data, file = input$opt$dataFile, sep = "\t", row.names = FALSE )
    }
    return( input$data$dat$flag )
}





###########################################################################





titan.bootstrap <-
function( input, reg, gene.valid, trace = TRUE )
{
    # Calculate bootstrap CI
    # by taking random sample of residuals from regression with replacement
    # and adding them to the fitted values
    # then redoing the regression
    # and calculating the 'bootstrap' roots

    if ( input$opt$R == 0 )
    {
        return( NULL )
    }
    if ( length( input$opt$gene1 ) == 0 )
    {
        warning( "No test genes selected" )
        return( NULL )
    }
    if ( length( input$opt$rx1 ) == 0 )
    {
        return( titan.bootstrap.baselineonly( input, reg, gene.valid, trace ) )
    }
    gene0 <- input$data$levels$gene[ match( input$opt$gene0, input$data$levels$gene ) ]
    gene1 <- input$data$levels$gene[ match( input$opt$gene1, input$data$levels$gene ) ]
    rx1 <- input$data$levels$rx[ match( input$opt$rx1, input$data$levels$rx ) ]
    nrx1 <- length( rx1 )
    if ( length( gene1 ) == 0 )
    {
        warning( "No test genes with valid data" )
        return( NULL )
    }
    if ( length( rx1 ) == 0 )
    {
        warning( "No comparison treatments with valid data" )
        return( NULL )
    }
    if ( trace )
    {
        cat( "\n\nRunning bootstrap\n")
    }

    # bootstrap statistic function
    boot.fun <-
    function( boot.data, indices, reg, orig.data, opt )
    {
        boot.data$y <- reg$fitted + reg$residuals[ indices ]
        boot.reg <- update( reg, data = boot.data, subset = NULL )
        orig.data$dat$y <- boot.data$y
        rx1 <- seq( orig.data$levels$rx )[ match( opt$rx1, orig.data$levels$rx ) ]
        gene0 <- seq( orig.data$levels$ngene )[ match( opt$gene0, orig.data$levels$gene ) ]
        gene1 <- seq( orig.data$levels$ngene )[ match( opt$gene1, orig.data$levels$gene ) ]
        nrx1 <- length( rx1 )
        ngene0 <- length( gene0 )
        ngene1 <- length( gene1 )
        theta.star.g0 <- matrix( NA, nr = ngene0, nc = nrx1 )
        theta.star.g1 <- matrix( NA, nr = ngene1, nc = nrx1 )
        if ( ngene0 )
        {
            for ( g in seq( ngene0 ) )
            {
                boot.r0 <- try(
                    titan.interpolate.line(
                        orig.data,
                        boot.reg,
                        0,
                        gene0[ g ],
                        1,
                        trace = FALSE
                    )
                )
                for ( r in seq( nrx1 ) )
                {
                    boot.rx <- try(
                        titan.interpolate.line(
                            orig.data,
                            boot.reg,
                            0,
                            gene0[ g ],
                            rx1[ r ],
                            trace = FALSE
                        )
                    )
                    if ( !inherits( boot.r0, "try-error" ) && !inherits( boot.rx, "try-error" ) )
                    {
                        theta.star.g0[ g, r ] <- boot.rx - boot.r0
                    }
                }
            }
        }
        for ( g in seq( ngene1 ) )
        {
            boot.r0 <- try(
                titan.interpolate.line(
                    orig.data,
                    boot.reg,
                    0,
                    gene1[ g ],
                    1,
                    trace = FALSE
                )
            )
            for ( r in seq( nrx1 ) )
            {
                boot.rx <- try(
                    titan.interpolate.line(
                        orig.data,
                        boot.reg,
                        0,
                        gene1[ g ],
                        rx1[ r ],
                        trace = FALSE
                    )
                )
                if ( !inherits( boot.r0, "try-error" ) && !inherits( boot.rx, "try-error" ) )
                {
                    theta.star.g1[ g, r ] <- boot.rx - boot.r0
                }
            }
        }
        ##cat( "theta.star =\n", theta.star, "\n" )
        if ( ngene0 )
        {
            theta.star.g1 <- theta.star.g1 - matrix( colMeans( theta.star.g0 ), nc = nrx1, nr = ngene1, byrow = TRUE )
        }
        return( as.vector( theta.star.g1 ) )
    }

    # restrict to subset of data
    orig.data <- input$data
    boot.data <- orig.data$dat <- reg$call$data[ reg$call$subset, ]

    set.seed( input$opt$seed )
    boot.res <- boot(
        boot.data,
        boot.fun,
        input$opt$R,
        orig.data = orig.data,
        reg = reg,
        opt = input$opt
    )
    boot.res$gene0 <- gene0
    boot.res$gene1 <- gene1
    boot.res$rx1 <- rx1

    # print bootstrap results
    if ( trace )
    {
        print.titanboot( boot.res, input$opt )
    }
    else
    {
        for ( i in 1:dim( boot.res$t )[ 2 ] )
        {
            if ( nMiss <- sum( is.na( boot.res$t[ , i ] ) ) )
            {
                warning(
                    nMiss, " missing values for gene ",
                    gene1[ i ],
                    call. = FALSE
                )
            }
        }
    }

    # bootstrap confidence intervals
    boot.res.ci <- NULL
    if ( length( input$opt$ciConf ) && length( input$opt$ciType ) )
    {
        if ( !any( is.na( boot.res$t ) ) )
        {
            n <- dim( boot.res$t )[ 2 ]
            boot.res.ci <- list(
                R = input$opt$R,
                t0 = NULL,
                call = NULL,
                normal = NULL,
                basic = NULL,
                percent = NULL,
                bca = NULL,
                gene = NULL,
                rx = NULL
            )
            if ( trace )
            {
                cat( "\nBOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS\n" )
            }
            for ( i in 1:n )
            {
                boot.res.ci$gene[ i ] <- gene1[ i %/% nrx1 ]
                boot.res.ci$rx[ i ] <- rx1[ i %% nrx1 + 1 ]
                if ( trace )
                {
                    cat( "\nGene", boot.res.ci$gene[ i ] )
                    cat( "\nRx  ", boot.res.ci$rx[ i ], "\n" )
                }
                boot.ci.old <- NULL
                boot.ci <- try(
                    boot.ci(
                        boot.res,
                        conf = input$opt$ciConf,
                        type = input$opt$ciType,
                        index = i
                    ),
                    silent = TRUE
                )
                if ( inherits( boot.ci, "try-error" ) && !is.na( match( "bca", input$opt$ciType ) ) )
                {
                    boot.ci.old <- boot.ci
                    ciType <- input$opt$ciType[ -match( "bca", input$opt$ciType ) ]
                    boot.ci <- try(
                        boot.ci(
                            boot.res,
                            conf = input$opt$ciConf,
                            type = ciType,
                            index = i
                        ),
                        silent = TRUE
                    )
                }
                if ( inherits( boot.ci, "try-error" ) )
                {
                    stop( boot.ci, call. = FALSE )
                }
                if (
                    inherits( boot.ci.old, "try-error" ) &&
                    length( grep( "extreme order statistics used as endpoints", boot.ci.old[ 1 ], ignore.case = TRUE ) ) > 0
                )
                {
                    warning(
                        "extreme order statistics used as endpoints\n",
                        "BCa confidence intervals could not be calculated: ",
                        "try increasing the number of bootstrap replicates",
                        call. = FALSE
                    )
                }
                else if ( inherits( boot.ci.old, "try-error" ) )
                {
                    warning(
                        "BCa confidence intervals could not be calculated",
                        call. = FALSE
                    )
                }
                if ( trace && !is.null( boot.ci ) )
                {
                    o <- options( digits = 7 )
                    cat( "\nFold change:" )
                    titan.print.bootci(
                        boot.ci,
                        h = function( x ) 10 ^ x,
                        hdot = function( x ) log( 10 ) * 10 ^ x
                    )
                    cat( "\nLog10 fold change:" )
                    titan.print.bootci( boot.ci )
                    options( digits = o$digits )
                }
                boot.res.ci$t0[ i ] <- boot.ci$t0
                boot.res.ci$normal[[ i ]] <- boot.ci$normal
                boot.res.ci$basic[[ i ]] <- boot.ci$basic
                boot.res.ci$percent[[ i ]] <- boot.ci$percent
                boot.res.ci$bca[[ i ]] <- boot.ci$bca
            }
        class( boot.res.ci ) <- "titanbootci"
        }
        else
        {
            warning(
                "Confidence intervals could not be calculated due to missing values",
                call. = FALSE
            )
        }
    }
    return( list( boot = boot.res, bootci = boot.res.ci ) )
}


titan.bootstrap.baselineonly <-
function( input, reg, gene.valid, trace = TRUE )
{
    gene0 <- input$data$levels$gene[ match( input$opt$gene0, input$data$levels$gene ) ]
    gene1 <- input$data$levels$gene[ match (input$opt$gene1, input$data$levels$gene ) ]
    ngene1 <- length( gene1 )
    if ( ngene1 < 2 )
    {
        warning( "Too few test genes with valid data" )
        return( NULL )
    }
    rx1 <- input$data$levels$rx[ match( input$opt$rx1, input$data$levels$rx ) ]
    nrx1 <- length( rx1 )
    if ( trace )
    {
        cat( "\n\nRunning bootstrap\n" )
    }

    boot.fun <-
    function( boot.data, indices, reg, orig.data, opt )
    {
        boot.data$y <- reg$fitted + reg$residuals[ indices ]
        boot.reg <- update( reg, data = boot.data, subset = NULL )
        orig.data$dat$y <- boot.data$y
        gene1 <- seq( orig.data$levels$ngene )[ match(opt$gene1, orig.data$levels$gene ) ]
        ngene1 <- length( gene1 )
        theta.star.g1 <- matrix( NA, nr = ( ngene1 - 1 ) * ngene1 / 2, nc = 1 )
        g1 <- 1
        g2 <- 1
        for ( g in seq( ( ngene1 - 1 ) * ngene1 / 2 ) )
        {
            g2 <- g2 + 1
            if ( g2 > ngene1 )
            {
                g1 <- g1 + 1
                g2 <- g1 + 1
            }
            boot.r0g1 <- try(
                titan.interpolate.line(
                    orig.data,
                    boot.reg,
                    0,
                    gene1[ g1 ],
                    1,
                    trace = FALSE
                )
            )
            boot.r0g2 <- try(
                titan.interpolate.line(
                    orig.data,
                    boot.reg,
                    0,
                    gene1[ g2 ],
                    1,
                    trace = FALSE
                )
            )
            if ( !inherits( boot.r0g1, "try-error" ) && !inherits( boot.r0g2, "try-error" ) )
            {
                theta.star.g1[ g, 1 ] <- boot.r0g1 - boot.r0g2
            }
        }
        return( as.vector( theta.star.g1 ) )
    }

    orig.data <- input$data
    boot.data <- orig.data$dat <- reg$call$data[ reg$call$subset, ]
    set.seed( input$opt$seed )
    boot.res <- boot(
        boot.data,
        boot.fun,
        input$opt$R,
        orig.data = orig.data,
        reg = reg,
        opt = input$opt
    )
    boot.res$gene0 <- gene0
    boot.res$gene1 <- gene1
    if (trace)
    {
        print.titanboot( boot.res, input$opt )
    }
    else
    {
        g1 <- 1
        g2 <- 1
        for ( i in 1:dim( boot.res$t )[ 2 ] )
        {
            g2 <- g2 + 1
            if ( g2 > ngene1 )
            {
                g1 <- g1 + 1
                g2 <- g1 + 1
            }
            if ( nMiss <- sum( is.na( boot.res$t[ , i ] ) ) )
            {
                if ( nrx1 > 0 )
                {
                    warning(
                        nMiss, " missing values for gene ", gene1[ i ],
                        call. = FALSE
                    )
                }
                else
                {
                    warning(
                        nMiss, " missing values for gene ", gene1[ g1 ], " / gene ", gene1[ g2 ],
                        call. = FALSE
                    )
                }
            }
        }
    }

    boot.res.ci <- NULL
    if ( length(input$opt$ciConf ) && length( input$opt$ciType ) )
    {
        if ( !any( is.na( boot.res$t ) ) )
        {
            n <- dim(boot.res$t)[ 2 ]
            boot.res.ci <- list(
                R = input$opt$R,
                t0 = NULL,
                call = NULL,
                normal = NULL,
                basic = NULL,
                percent = NULL,
                bca = NULL,
                gene1 = NULL,
                gene2 = NULL
            )
            if (trace)
            {
                cat( "\nBOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS\n" )
            }
            g1 <- 1
            g2 <- 1
            for ( i in 1:n )
            {
                g2 <- g2 + 1
                if ( g2 > ngene1 )
                {
                    g1 <- g1 + 1
                    g2 <- g1 + 1
                }
                boot.res.ci$gene1[ i ] <- gene1[ g1 ]
                boot.res.ci$gene2[ i ] <- gene1[ g2 ]
                if (trace)
                {
                    cat(
                        "\n",
                        "Gene",
                        boot.res.ci$gene1[ i ],
                        "/ Gene",
                        boot.res.ci$gene2[i],
                        "\n"
                    )
                }
                boot.ci.old <- NULL
                boot.ci <- try(
                    boot.ci(
                        boot.res,
                        conf = input$opt$ciConf,
                        type = input$opt$ciType,
                        index = i
                    ),
                    silent = TRUE
                )
                if ( inherits( boot.ci, "try-error" ) && !is.na( match( "bca", input$opt$ciType ) ) )
                {
                    boot.ci.old <- boot.ci
                    ciType <- input$opt$ciType[ -match( "bca", input$opt$ciType ) ]
                    boot.ci <- try(
                        boot.ci(
                            boot.res,
                            conf = input$opt$ciConf,
                            type = ciType,
                            index = i
                        ),
                        silent = TRUE
                    )
                }
                if ( inherits( boot.ci, "try-error" ) )
                {
                    stop( boot.ci, call. = FALSE )
                }
                if (
                    inherits( boot.ci.old, "try-error" ) &&
                    length( grep( "extreme order statistics used as endpoints", boot.ci.old[1], ignore.case = TRUE ) ) > 0
                )
                {
                    warning(
                        "extreme order statistics used as endpoints\n",
                        "BCa confidence intervals could not be calculated: ",
                        "try increasing the number of bootstrap replicates",
                        call. = FALSE
                    )
                }
                else if ( inherits( boot.ci.old, "try-error" ) )
                {
                    warning(
                        "BCa confidence intervals could not be calculated",
                        call. = FALSE
                    )
                }
                if ( trace && !is.null( boot.ci ) )
                {
                    o <- options( digits = 7 )
                    cat( "\nFold change:" )
                    titan.print.bootci(
                        boot.ci,
                        h = function( x ) { 10^x },
                        hdot = function( x ) { log( 10 ) * 10^x }
                    )
                    cat( "\nLog10 fold change:" )
                    titan.print.bootci( boot.ci )
                    options( digits = o$digits )
                }
                boot.res.ci$t0[ i ] <- boot.ci$t0
                boot.res.ci$normal[[ i ]] <- boot.ci$normal
                boot.res.ci$basic[[ i ]] <- boot.ci$basic
                boot.res.ci$percent[[ i ]] <- boot.ci$percent
                boot.res.ci$bca[[ i ]] <- boot.ci$bca
            } # next i
            class( boot.res.ci ) <- "titanbootci"
        }
        else
        {
            warning(
                "Confidence intervals could not be calculated due to missing values",
                call. = FALSE
            )
        }
    }
    return( list( boot = boot.res, bootci = boot.res.ci ) )
}





###########################################################################




# print method

print.titan <-
function( x, ... )
{
    titan.printreg( x$reg )
    titan.print.interpolation( x$interpolation )
    if ( !is.null( x$boot ) )
    {
        print.titanboot( x$boot, x$opt )
        if ( !is.null( x$bootci ) )
        {
            print.titanbootci( x$bootci )
        }
    }
}

titan.printreg <-
function( reg )
{
    r <- summary( reg, correlation = FALSE )
    if ( !is.null( coef( r ) ) )
    {
        rownames( r$coefficients ) <- titan.regsub( rownames( coef( r ) ), reg$levels )
    }
    r$call <- titan.regformula( reg$call$formula, reg$levels )
    print( r )
    o <- options()$digits
    options( digits = 4 )
    cat( "\nRegression results per gene:\n" )
    print( reg$gene.res )
    cat( "\nR-Squared per gene per treatment:\n" )
    print( reg$generx.res[ , , "R-Squared" ] )
    options( digits = o )
    cat( "\n" )
}

titan.regsub <-
function( s, levels )
{
    for ( g in seq( length( levels$gene ) ) )
    {
        s <- gsub( paste( "genei\\[, ", g, "]", sep = "" ), levels$gene[ g ], s )
    }
    for ( r in seq( length( levels$rx ) ) )
    {
        s <- gsub( paste( "contr.rxi\\[, ", r - 1, "]", sep = "" ), levels$rx[ r ], s )
    }
    for ( gr in seq( length( levels$generx ) ) )
    {
        s <- gsub( paste( "generxi\\[, ", gr, "]", sep = "" ), levels$generx[ gr ], s )
    }
    s <- gsub( "polyx1", "conc", s )
    s <- gsub( "polyx2", "conc^2", s )
    ## s <- gsub( "ns", "spline", s )
    return( s )
}

titan.regformula <-
function( f, levels )
{
    f <- as.character( f )
    f[3] <- titan.regsub( f[3], levels )
    return( paste( f[2], f[1], f[3] ) )
}

titan.print.interpolation <-
function( interpolation )
{
    cat( "Interpolated Concentrations\n" )
    print( 10 ^ interpolation$roots )
    cat( "\nLog10 Interpolated Concentrations\n" )
    print( interpolation$roots )
    if ( !is.null( interpolation$log10fold ) ) {
        cat( "\nFold Change\n" )
        print( 10 ^ interpolation$log10fold )
        cat( "\nLog10 Fold Change\n" )
        print( interpolation$log10fold )
        cat( "\nFold changes calculated relative to baseline rx: ")
        cat( interpolation$rx0 )
        if ( length( interpolation$gene0 ) )
        {
            cat( "\nFold changes adjusted for housekeeping (control) genes: ")
            cat( paste( interpolation$gene0 ) )
        }
    }
    cat( "\n" )
}

print.titanboot <-
function( boot.res, opt )
{
    gene0 <- boot.res$gene0
    gene1 <- boot.res$gene1
    ngene <- length( gene0 ) + length( gene1 )
    rx1 <- boot.res$rx1
    nrx1 <- length( rx1 )
    ngene1 <- length( gene1 )
    len <- length( boot.res$t0 )
    boot.mat <- matrix( NA, nr = len, nc = 7 )
    boot.mat[ , 1 ] <- 10 ^ boot.res$t0
    boot.mat[ , 2 ] <- boot.res$t0
    boot.mat[ , 3 ] <- colMeans( boot.res$t, na.rm = TRUE ) - boot.res$t0
    for ( i in 1:len )
    {
        ti <- na.omit( boot.res$t[ , i ] )
        if ( length( ti ) > 1 )
        {
            boot.mat[ i, 4 ] <- sd( ti )
        }
        else
        {
            is.na( boot.mat[ i, 4 ] ) <- TRUE
        }
    }
    boot.mat[ , 5 ] <- 1 - ( rowSums( sign( t( boot.res$t ) ) == sign( boot.res$t0 ) ) / opt$R )
    if ( nrx1 > 0 )
    {
        boot.mat[ , 6 ] <- ( seq( len ) - 1 ) %% ngene1 + 1
        boot.mat[ , 7 ] <- ( seq( len ) - 1 ) %/% ngene1 + 1
        colnames( boot.mat ) <- c( "Fold  ", "Log10 Fold", "Bias  ", "Std.Error", "p ", "Gene", "Rx" )
        rn <- paste( "Gene ", gene1[ boot.mat[ , 6 ] ], ", Rx ", rx1[ boot.mat[ , 7 ] ], sep = "" )
    }
    else
    {
        g1 <- 1
        g2 <- 1
        for ( i in 1:len )
        {
            g2 <- g2 + 1
            if ( g2 > ngene1 )
            {
              g1 <- g1 + 1
              g2 <- g1 + 1
            }
            boot.mat[ i, 6 ] <- g1
            boot.mat[ i, 7 ] <- g2
        }
        colnames( boot.mat ) <- c( "Fold  ", "Log10 Fold", "Bias  ", "Std.Error", "p ", "Gene", " / Gene" )
        rn <- paste( "Gene ", gene1[ boot.mat[ , 6 ] ], " / Gene ", gene1[ boot.mat[ , 7 ] ], sep = "" )
    }
    rownames( boot.mat ) <- rn
    cat( "\nORDINARY NONPARAMETRIC BOOTSTRAP\n" )
    cat( "Based on", opt$R, "bootstrap replicates\n" )
    cat( "Seed for random number generator: " )
    cat( opt$seed )
    cat( "\n\n" )
    cat( "Bootstrap statistics:\n" )
    print( boot.mat[ , 1:5 ] )
    cat( "\np-value is one-sided = Pr( fold change differs " )
    cat( "from 1 and is in the expected direction )\n" )
    for ( i in 1:dim( boot.res$t )[ 2 ] )
    {
        if ( nMiss <- sum( is.na( boot.res$t[ , i ] ) ) )
        {
            if ( nrx1 > 0 )
            {
                warning(
                    nMiss, " missing values for gene ",
                    boot.res$gene1[ ( i - 1 ) %/% nrx1 + 1 ],
                    call. = FALSE
                )
            }
            else
            {
                warning(
                    nMiss,
                    " missing values for gene ",
                    boot.res$gene1[ boot.mat[ i, 6 ] ],
                    " / gene ",
                    boot.res$gene1[ boot.mat[ i, 7 ] ],
                    call. = FALSE
                )
            }
        }
    }
}

print.titanbootci <-
function( x, ... )
{
    if ( ( len <- length( x ) ) < 6 )
    {
        return
    }
    cat( "\nBOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS\n" )
    o <- options( digits = 7 )
    for ( i in 1:length( x$t0 ) )
    {
    if (is.null(x$rx))
    {
        cat( "\n", "Gene ", x$gene1[ i ], " / Gene ", x$gene2[ i ] )
    }
    else
    {
        cat( "\nGene", x$gene[ i ], "\nRx  ", x$rx[ i ] )
    }
    cat( "\n\nFold change:" )
        boot.ci <- list( R = x$R, call = NULL, t0 = x$t0[ i ] )
        k <- 4
        for ( j in 4:( len - 2 ) )
        {
            if ( !is.null( x[[ j ]] ) )
            {
                boot.ci[[ k ]] <- x[[ j ]][[ i ]]
                names( boot.ci )[ k ] <- names( x )[ j ]
                k <- k + 1
            }
        }
        titan.print.bootci(
            boot.ci,
            h = function( x ) 10 ^ x,
            hdot = function( x ) log( 10 ) * 10 ^ x
        )
        cat( "\nLog10 fold change:" )
        titan.print.bootci( boot.ci )
    }
    options( digits = o$digits )
}

    # hacked `print.bootci' method
    titan.print.bootci <-
    function( x, hinv = NULL, ... )
    {
    #
    # Print the output from boot.ci
    #
        ci.out <- x
        cl <- ci.out$call
        ntypes <- length( ci.out ) - 3
        nints <- nrow( ci.out[[ 4 ]] )
        t0 <- ci.out$t0
        if ( !is.null( hinv ) )
        {
            t0 <- hinv( t0 )
        }
    # Find the number of decimal places which should be used
        o <- getOption( "digits" )
        digs <- ceiling( log10( abs( t0 ) ) )
        if ( digs <= 0 ) digs <- o
        else if ( digs >= o ) digs <- 0
        else digs <- o - digs
        intlabs <- NULL
        basrg <- strg <- perg <- bcarg <- NULL
        sp <- paste( rep(  " ", o  ), collapse = ""  )
        if ( !is.null( ci.out$normal ) )
        {
            intlabs <- c( intlabs, paste( sp, "  Normal   ", sp, sep ="" ) )
        }
        if ( !is.null( ci.out$basic ) )
        {
            intlabs <- c( intlabs, paste( sp, "  Basic    ", sp, sep = "" ) )
            basrg <- range( ci.out$basic[ , 2:3 ] )
        }
        if ( !is.null( ci.out$student ) )
        {
            intlabs <- c( intlabs, paste( sp, "Studentized", sp, sep = "" ) )
            strg <- range( ci.out$student[ , 2:3 ] )
        }
        if ( !is.null( ci.out$percent ) )
        {
            intlabs <- c( intlabs, paste( sp, " Percentile", sp, sep = "" ) )
            perg <- range( ci.out$percent[ , 2:3 ] )
        }
        if ( !is.null( ci.out$bca ) )
        {
            intlabs <- c( intlabs, paste( sp, "   BCa     ", sp, sep = "" ) )
            bcarg <- range( ci.out$bca[ , 2:3 ] )
        }
        level <- 100 * ci.out[[ 4 ]][ , 1 ]
        if ( ntypes == 4 )
        {
            n1 <- n2 <- 2
        }
        else if ( ntypes == 5 )
        {
            n1 <- 3
            n2 <- 2
        }
        else {
            n1 <- ntypes
            n2 <- 0
        }
        ints1 <- matrix( NA, nints, 2 * n1 + 1 )
        ints1[ , 1 ] <- level
        n0 <- 4
    # Re-organize the intervals and coerce them into character data
        for ( i in n0:( n0 + n1 - 1 ) )
        {
            j <- c( 2 * i - 6, 2 * i - 5 )
            nc <- ncol( ci.out[[ i ]] )
            nc <- c( nc - 1, nc )
            if ( is.null( hinv ) )
            {
                ints1[ , j ] <- ci.out[[ i ]][ , nc ]
           }
           else
           {
                ints1[ , j ] <- hinv( ci.out[[ i ]][ , nc ] )
            }
        }
        n0 <- 4 + n1
        ints1 <- format( round( ints1,digs ) )
        ints1[ , 1 ] <- paste( "\n", level, "%  ", sep = "" )
        ints1[ , 2*( 1:n1 ) ] <- paste( "( ",ints1[ , 2 * ( 1:n1 ) ], ",", sep = "" )
        ints1[ , 2*( 1:n1 ) + 1 ] <- paste( ints1[ , 2*( 1:n1 ) + 1 ]," )  " )
        if ( n2 > 0 )
        {
            ints2 <- matrix( NA, nints, 2 * n2 + 1 )
            ints2[ , 1 ] <- level
            j <- c( 2, 3 )
            for ( i in n0:( n0 + n2 - 1 ) )
            {
                if ( is.null( hinv ) )
                {
                    ints2[ ,j ] <- ci.out[[ i ]][ , c( 4, 5 ) ]
                }
                else
                {
                    ints2[ ,j ] <- hinv( ci.out[[ i ]][ , c( 4, 5 ) ] )
                }
                j <- j + 2
            }
            ints2 <- format( round( ints2,digs ) )
            ints2[ , 1 ] <- paste( "\n", level, "%  ", sep = "" )
            ints2[ , 2*( 1:n2 ) ] <- paste( "( ", ints2[ , 2 * ( 1:n2 ) ], ",", sep = "" )
            ints2[ , 2*( 1:n2 ) + 1 ] <- paste( ints2[ , 2 * ( 1:n2 ) + 1 ], " )  " )
        }
        R <- ci.out$R
    # Print the intervals
        cat( "\nLevel", intlabs[ 1:n1 ] )
        cat( t( ints1 ) )
        cat( "\n" )
        if ( n2 > 0 )
        {
            cat( "\nLevel", intlabs[ ( n1 + 1 ):( n1 + n2 ) ] )
            cat( t( ints2 ) )
            cat( "\n" )
        }
    # Print any warnings about extreme values.
        if ( !is.null( basrg ) )
        {
            if ( ( basrg[ 1 ] <= 1 ) || ( basrg[ 2 ] >= R ) )
            {
                cat( "Warning: Basic Intervals used Extreme Quantiles\n" )
            }
            if ( ( basrg[ 1 ] <= 10 ) || ( basrg[ 2 ] >= R-9 ) )
            {
                cat( "Some basic intervals may be unstable\n" )
            }
        }
        if ( !is.null( strg ) )
        {
            if ( ( strg[ 1 ] <= 1 ) || ( strg[ 2 ] >= R ) )
            {
                cat( "Warning: Studentized Intervals used Extreme Quantiles\n" )
            }
            if ( ( strg[ 1 ] <= 10 ) || ( strg[ 2 ] >= R-9 ) )
            {
                cat( "Some studentized intervals may be unstable\n" )
            }
        }
        if ( !is.null( perg ) ) {
            if ( ( perg[ 1 ] <= 1 ) || ( perg[ 2 ] >= R ) )
            {
                cat( "Warning: Percentile Intervals used Extreme Quantiles\n" )
            }
            if ( ( perg[ 1 ] <= 10 ) || ( perg[ 2 ] >= R-9 ) )
            {
                cat( "Some percentile intervals may be unstable\n" )
            }
        }
        if ( !is.null( bcarg ) ) {
            if ( ( bcarg[ 1 ] <= 1 ) || ( bcarg[ 2 ] >= R ) )
            {
                cat( "Warning: BCa Intervals used Extreme Quantiles\n" )
            }
            if ( ( bcarg[ 1 ] <= 10 ) || ( bcarg[ 2 ] >= R - 9 ) )
            {
                cat( "Some BCa intervals may be unstable\n" )
            }
        }
        invisible( ci.out )
    }





###########################################################################




# main function


    # set error message options
    # [commented out on instructions of B.Ripley]
    ##o <- options()
    ##options( warn = 1, error = recover )

    if ( trace )
    {
        cat( titan.version )
        cat( "\n" )
    }

    # load required libraries, installing if necessary
    require( "MASS" )
    require( "tcltk" )
    require( "splines" )
    require( "boot" )
    require( "lattice" )

    # get data and options from input
    titan.ip <- titan.input(
        data,
        trace,
        widget,
        dataFile,
        outFile,
        pdfFile,
        flagRaw,
        flagFitted,
        freqLo,
        freqHi,
        reg,
        term,
        sel,
        alpha,
        rx0,
        rx1,
        gene0,
        gene1,
        R,
        seed,
        ciConf,
        ciType
    )

    # record output
    if ( titan.ip$opt$outFile != "" )
    {
        txt <- file( titan.ip$opt$outFile, open = "w" )
        sink( txt )
        sink( type = "message" )
    }

    # run regression
    titan.cf <- titan.create.formulae( titan.ip )
    if ( length( titan.cf$gene.valid ) == 0 )
    {
        return( NULL )
    }
    titan.reg <- titan.regression( titan.ip, titan.cf, trace )
    if ( is.null( titan.reg ) )
    {
        return( NULL )
    }

    # manually flag data points
    while (
        widget &&
        interactive() &&
        (
            tclvalue
            (
                tkmessageBox(
                    message = "Review flagged spots?",
                    type = "yesno",
                    default = "yes",
                    icon = "question",
                    title = titan.version
                )
            ) == "yes"
        )
    )
    {
        titan.ip$data$dat$flag <- titan.manualFlag( titan.ip, titan.reg, titan.cf )
        titan.ip$data <- titan.process.data( titan.ip$data )
        titan.cf <- titan.create.formulae( titan.ip )
        titan.reg <- titan.regression( titan.ip, titan.cf, trace )
        if ( is.null( titan.reg ) )
        {
            break
            return( NULL )
        }
    }

    # interpolate x values at y == 0
    titan.interpolation <- titan.interpolate(
        titan.ip,
        titan.reg,
        titan.cf$gene.valid,
        titan.cf$generx.valid,
        trace
    )
    #
    # restrict test treatments & housekeeping genes to those with valid data
    titan.ip$opt$gene0 <- titan.interpolation$gene0
    titan.ip$opt$rx1 <- titan.interpolation$rx1

    # plot results
    titan.plot( titan.ip, titan.reg, titan.cf, trace = trace )

    # run bootstrap function
    titan.boot.res <- titan.bootstrap( titan.ip, titan.reg, titan.cf$gene.valid, trace )

    #restore options
    ##options( o )

    # return object of class 'titan'
    titan.res <- list(
        data = titan.ip$data$dat[ , c( "freq", "gene", "rx", "comp.conc", "flag" ) ],
        opt = titan.ip$opt,
        reg = titan.reg,
        cf = titan.cf,
        interpolation = titan.interpolation,
        boot = titan.boot.res$boot,
        bootci = titan.boot.res$bootci
    )
    class( titan.res ) <- "titan"
    if ( sink.number() )
    {
        sink()
    }

    return( invisible( titan.res ) )
}

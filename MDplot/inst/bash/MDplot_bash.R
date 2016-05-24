require( MDplot )

# get arguments and look for function call
VEC_inputArguments <- commandArgs( trailingOnly = TRUE )
if( length( VEC_inputArguments ) < 2 )
{
  stop( "Error due to missing arguments. You need to supply at least a plot selection and an input file." )
}
STRING_function <- VEC_inputArguments[ 1 ]
LIST_arguments <- parse_arguments( VEC_inputArguments[ -1 ] ) # -1: function name excluding
#########

# prepare vectors for checks
VEC_requiredForAll <- c(  )
VEC_allowedForAll <- c( VEC_requiredForAll, "files", "size", "outformat",
                        "outfile", "title", "subtitle",
                        "enableProtocol", "colours", "resolution",
                        "axisNames", "names", "printLegend",
                        "mdEngine", "help" )
VEC_allowedForAllDesc <- c( "<input file(s), separated by ','>", "<dimensions of the plot> (optional)", "['png'/'pdf'/'tiff'] (optional)",
                            "<outputfile> (optional)", "<plot main title> (optional)", "<plot subtitle> (optional)",
                            "<protocol steps in plot generation> (optional)", "<vector of colours used, separated by ','> (optional)", "<resolution> (optional)",
                            "<vector of names for the axes> (optional)", "<vector of names for the data sets> (optional)", "['true'/'false'] (optional)",
                            "<name of molecular dynamics engine used> (default: GROMOS)",
                            "<if set to 'true', all other options are ignored and help is printed> (optional)" )
#########

# set settings for all plots to be followed
BOOL_printLegend = TRUE
if( isKeySet( LIST_arguments, "printLegend" ) )
  if( getValue( LIST_arguments, "printLegend" ) == "FALSE" )
    BOOL_printLegend = FALSE
VEC_size <- c( 640, 640 )
if( isKeySet( LIST_arguments, "size" ) )
  VEC_size <- unlist( strsplit( getValue( LIST_arguments, "size" ),
                                ",",
                                fixed = TRUE ) )
VEC_dataNames <- NA
if( isKeySet( LIST_arguments, "names" ) )
  VEC_dataNames <- unlist( strsplit( getValue( LIST_arguments, "names" ),
                                ",",
                                fixed = TRUE ) )
VEC_axisNames <- NA
if( isKeySet( LIST_arguments, "axisNames" ) )
  VEC_axisNames <- unlist( strsplit( getValue( LIST_arguments, "axisNames" ),
                                     ",",
                                     fixed = TRUE ) )
STRING_outformat <- "png"
if( isKeySet( LIST_arguments, "outformat" ) )
  STRING_outformat <- getValue( LIST_arguments, "outformat" )
STRING_outfile <- "MDplot_out"
if( isKeySet( LIST_arguments, "outfile" ) )
  STRING_outfile <- getValue( LIST_arguments, "outfile" )
REAL_resolution <- 150
if( isKeySet( LIST_arguments, "resolution" ) )
  REAL_resolution <- as.numeric( getValue( LIST_arguments, "resolution" ) )
STRING_mdEngine <- "GROMOS"
if( isKeySet( LIST_arguments, "mdEngine" ) )
  STRING_mdEngine <- as.numeric( getValue( LIST_arguments, "mdEngine" ) )
#########

# define plot device and options
if( STRING_outformat == "pdf" )
{
  pdf( file = paste( STRING_outfile ),
       width = as.numeric( VEC_size[ 1 ] ) / 96,
       height = as.numeric( VEC_size[ 2 ] ) / 96 )
} else {
  if( STRING_outformat == "png" )
  {
    png( paste( STRING_outfile ),
         width = as.numeric( VEC_size[ 1 ] ),
         height = as.numeric( VEC_size[ 2 ] ),
         units = "px",
         res = REAL_resolution,
         type = "cairo" )
  } else {
    if( STRING_outformat == "tiff" )
    {
      tiff( filename = paste( STRING_outfile ),
            width = as.numeric( VEC_size[ 1 ] ),
            height = as.numeric( VEC_size[ 2 ] ),
            units = "px",
            compression = "none",
            res = REAL_resolution )
    } else {
      stop( paste( "Error, the specified output format '",
                   STRING_outformat,
                   "' is not known.",
                   sep = "" ) )
    }
  }
}
#########

# check, which plot has been selected
if( STRING_function == "dssp_summary" )
{
  # check, if input is sane for this plot and get input files
  VEC_dsspSumAll <- c( "plotType", "showResidues", "showValues" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_dsspSumAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_dsspSumAll,
                   VEC_allowedForAll ),
                c( "['dots'/'curves'/'bars'] (optional)",
                   "<range of residues to show, separated by ','> (optional)",
                   "<range of values to show, separated by ','> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  VEC_residues <- NA
  if( isKeySet( LIST_arguments, "showResidues" ) )
    VEC_residues <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "showResidues" ),
                                                  ",",
                                                  fixed = TRUE ) ) )
  VEC_values <- NA
  if( isKeySet( LIST_arguments, "showValues" ) )
    VEC_values <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "showValues" ),
                                                ",",
                                                fixed = TRUE ) ) )
  
  
  # plot
  MDplot::dssp_summary( MDplot::load_dssp_summary( VEC_files,
                                                   mdEngine = STRING_mdEngine ),
                        printLegend = BOOL_printLegend,
                        showResidues = VEC_residues,
                        showValues = VEC_values,
                        plotType = ifelse( isKeySet( LIST_arguments, "plotType" ),
                                           getValue( LIST_arguments, "plotType" ),
                                           "dots" ),
                        main = ifelse( isKeySet( LIST_arguments, "title" ),
                                       getValue( LIST_arguments, "title" ),
                                       NA ) )
}



if( STRING_function == "dssp_ts" )
{
  # check, if input is sane for this plot and get input files
  VEC_dsspTSAll <- c( "timeBoundaries", "residueBoundaries", "timeUnit", "snapshotsPerTimeInt" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_dsspTSAll,
                  VEC_allowedForAll ), LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_dsspTSAll,
                   VEC_allowedForAll ),
                c( "<range of time plotted, separated by ','> (optional)",
                   "<range of residues plotted, separated by ','> (optional)",
                   "<time unit, often: 'ns'> (optional)",
                   "<snapshots per time unit> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  VEC_timeBoundaries <- NA
  if( isKeySet( LIST_arguments, "timeBoundaries" ) )
    VEC_timeBoundaries <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "timeBoundaries" ),
                                                        ",",
                                                        fixed = TRUE )))
  VEC_residueBoundaries <- NA
  if( isKeySet( LIST_arguments, "residueBoundaries" ) )
    VEC_residueBoundaries <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "residueBoundaries" ),
                                                           ",",
                                                           fixed = TRUE )))
  
  # plot
  MDplot::dssp_ts( MDplot::load_dssp_ts( VEC_files ), 
                   timeBoundaries = VEC_timeBoundaries,
                   residueBoundaries = VEC_residueBoundaries,
                   timeUnit = ifelse( isKeySet( LIST_arguments, "timeUnit" ),
                                      getValue( LIST_arguments, "timeUnit" ),
                                      NA ),
                   snapshotsPerTimeInt = ifelse( isKeySet( LIST_arguments, "snapshotsPerTimeInt" ),
                                                 as.numeric( getValue( LIST_arguments, "snapshotsPerTimeInt" ) ),
                                                 NA ),
                   main = ifelse( isKeySet( LIST_arguments, "title" ),
                                  getValue( LIST_arguments, "title" ),
                                  NA ) )
}



if( STRING_function == "xrmsd" )
{
  # check, if input is sane for this plot and get input files
  VEC_xrmsdAll <- c( "xaxisRange", "yaxisRange", "factor" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_xrmsdAll,
                  VEC_allowedForAll ), LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_xrmsdAll,
                   VEC_allowedForAll ),
                c( "<range of x-axis data points plotted, separated by ','> (optional)",
                   "<range of y-axis data points plotted, separated by ','> (optional)",
                   "<factor by which the values should be divided> (default: 10000)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  VEC_xaxisRange <- NA
  if( isKeySet( LIST_arguments, "xaxisRange" ) )
    VEC_xaxisRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "xaxisRange" ),
                                                    ",",
                                                    fixed = TRUE ) ) )
  VEC_yaxisRange <- NA
  if( isKeySet( LIST_arguments, "yaxisRange" ) )
    VEC_yaxisRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "yaxisRange" ),
                                                    ",",
                                                    fixed = TRUE ) ) )
  
  # plot
  MDplot::xrmsd( MDplot::load_xrmsd( VEC_files, factor = ifelse( isKeySet( LIST_arguments, "factor" ),
                                                                 as.numeric( getValue( LIST_arguments, "factor" ) ),
                                                                 10000 ) ),
                 main = ifelse( isKeySet( LIST_arguments, "title" ),
                                getValue( LIST_arguments, "title" ),
                                NA ),
                 xaxisRange = VEC_xaxisRange,
                 yaxisRange = VEC_yaxisRange )
}



if( STRING_function == "rmsf" )
{
  # check, if input is sane for this plot and get input files
  VEC_rmsfAll <- c( "residuewise", "rmsfUnit", "range" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_rmsfAll,
                  VEC_allowedForAll ), LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_rmsfAll,
                   VEC_allowedForAll ),
                c( "<specifies, whether the protein is given in atoms or residues> (optional)",
                   "<abbreviation of rmsf unit> (default: nm)",
                   "<range of residues to be plotted> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  BOOL_resWise <- FALSE
  if( isKeySet( LIST_arguments, "residuewise" ) )
    BOOL_resWise <- ifelse( ( getValue( LIST_arguments, "residuewise" ) == "TRUE" ),
                            TRUE,
                            FALSE )
  VEC_range <- NA
  if( isKeySet( LIST_arguments, "range" ) )
    VEC_range <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "range" ),
                                               ",",
                                               fixed = TRUE ) ) )
  
  # plot
  MDplot::rmsf( MDplot::load_rmsf( VEC_files ),
                names = VEC_dataNames,
                residuewise = BOOL_resWise,
                printLegend = BOOL_printLegend,
                range = VEC_range,
                rmsfUnit = ifelse( isKeySet( LIST_arguments, "rmsfUnit" ),
                                   getValue( LIST_arguments, "rmsfUnit" ),
                                   "nm" ),
                main = ifelse( isKeySet( LIST_arguments, "title" ),
                               getValue( LIST_arguments, "title" ),
                               NA ) )
}



if( STRING_function == "rmsd" )
{
  # check, if input is sane for this plot and get input files
  VEC_rmsdAll <- c( "snapshotsPerTimeInt", "timeUnit", "rmsdUnit" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_rmsdAll,
                  VEC_allowedForAll ), LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_rmsdAll,
                   VEC_allowedForAll ),
                c( "<factor by which the values should be divided> (default: 1000)",
                   "<abbreviation of unit used for time-axis> (default: ns)",
                   "<abbreviation of unit used for rmsd-axis> (default: nm)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }

  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }    
  
  # plot
  MDplot::rmsd( MDplot::load_rmsd( VEC_files ),
                names = VEC_dataNames,
                snapshotsPerTimeInt = ifelse( isKeySet( LIST_arguments, "snapshotsPerTimeInt" ),
                                              as.numeric( getValue( LIST_arguments, "snapshotsPerTimeInt" ) ),
                                              1000 ),
                timeUnit = ifelse( isKeySet( LIST_arguments, "timeUnit" ),
                                   getValue( LIST_arguments, "timeUnit" ),
                                   "ns" ),
                rmsdUnit = ifelse( isKeySet( LIST_arguments, "rmsdUnit" ),
                                   getValue( LIST_arguments, "rmsdUnit" ),
                                   "nm" ),
                main = ifelse( isKeySet( LIST_arguments, "title" ),
                               getValue( LIST_arguments, "title" ),
                               NA ) )
}



if( STRING_function == "MDplot_RMSD_average" )
{
  # check, if input is sane for this plot and get input files
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( VEC_allowedForAll, LIST_arguments )

  # plot
  MDplot_RMSD_average( LIST_input = list( list( name = "one", 
                                                files = getValue( LIST_arguments, "files" ) ) ),
                       main = ifelse( isKeySet( LIST_arguments, "title" ),
                                      getValue( LIST_arguments, "title" ),
                                      NA ) )
}



if( STRING_function == "ramachandran" )
{
  # check, if input is sane for this plot and get input files
  VEC_ramaAll <- c( "bins",
                    "angleColumns",
                    "plotType",
                    "heatFun",
                    "heatUnits",
                    "plotContour",
                    "shiftAngles",
                    VEC_allowedForAll )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                VEC_ramaAll,
                c( "<number of bins to divide data in (x, y)> (optional)",
                   "<columns in file containing dihedrals>",
                   "['sparse'/'comic'/'fancy'] (optional)",
                   "<function to treat heat with, default 'log'> (optional)",
                   "<units, in which heat is given> (optional)",
                   "<plot contour as well, default 'false'> (optional)",
                   "<if angle interval is not -180 to 180, a shift can be specified> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( VEC_ramaAll, LIST_arguments )
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  VEC_bins <- c( 450, 450 )
  if( isKeySet( LIST_arguments, "bins" ) )
    VEC_bins <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "bins" ),
                                              ",",
                                              fixed = TRUE ) ) )
  VEC_angleColumns <- c( 1, 2 )
  if( isKeySet( LIST_arguments, "angleColumns" ) )
    VEC_angleColumns <- as.numeric( unlist( strsplit( getValue( LIST_arguments, "angleColumns" ),
                                                      ",",
                                                      fixed = TRUE ) ) )
  STRING_heatfunction <- "log"
  if( isKeySet( LIST_arguments, "heatFun" ) )
    STRING_heatfunction <- getValue( LIST_arguments, "heatFun" )
  
  STRING_plotType <- "comic"
  if( isKeySet( LIST_arguments, "plotType" ) )
    STRING_plotType <- getValue( LIST_arguments, "plotType" )
  
  STRING_heatUnits <- NA
  if( isKeySet( LIST_arguments, "heatUnits" ) )
    STRING_heatUnits <- paste( "[", getValue( LIST_arguments, "heatUnits" ), "]", sep = "" )

  # plot
  MDplot::ramachandran( MDplot::load_ramachandran( VEC_files,
                                                   angleColumns = VEC_angleColumns,
                                                   shiftAngles = ifelse( isKeySet( LIST_arguments, "shiftAngles" ),
                                                                         as.numeric( getValue( LIST_arguments, "shiftAngles" ) ),
                                                                         NA ) ),
                        xBins = VEC_bins[ 1 ],
                        yBins = VEC_bins[ 2 ],
                        plotType = STRING_plotType,
                        heatFun = STRING_heatfunction,
                        printLegend = ifelse( isKeySet( LIST_arguments, "printLegend" ),
                                              TRUE,
                                              FALSE ),
                        plotContour = ifelse( isKeySet( LIST_arguments, "plotContour" ),
                                              TRUE,
                                              FALSE ),
                        heatUnits = STRING_heatUnits,
                        main = ifelse( isKeySet( LIST_arguments, "title" ),
                                       getValue( LIST_arguments, "title" ),
                                       NA ) )
}



if( STRING_function == "TIcurve" )
{
  # check, if input is sane for this plot and get input files
  VEC_TIcAll <- c( "invertedBackwards", "printValues", "printErrors", "errorBarThreshold" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_TIcAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_TIcAll,
                   VEC_allowedForAll),
                c( "<if 'TRUE', the backward points are inverted> (optional)",
                   "<sets, whether the values are to be plotted> (default: TRUE)",
                   "<sets, whether the error bars are to be plotted> (default: TRUE)",
                   "<sets a threshold, below which the error bars are not plotted> (default: 0)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  BOOL_inverted = FALSE
  if( isKeySet( LIST_arguments, "invertedBackwards" ) )
    BOOL_inverted = ifelse( getValue( LIST_arguments, "invertedBackwards" ) == "TRUE",
                            TRUE,
                            FALSE )
  
  # plot
  MDplot::TIcurve( MDplot::load_TIcurve( VEC_files ),
                   invertedBackwards = BOOL_inverted,
                   printValues = ifelse( isKeySet( LIST_arguments, "printValues" ),
                                         FALSE,
                                         TRUE ),
                   printErrors = ifelse( isKeySet( LIST_arguments, "printErrors" ),
                                         FALSE,
                                         TRUE ),
                   errorBarThreshold = ifelse( isKeySet( LIST_arguments, "errorBarThreshold" ),
                                               getValue( LIST_arguments, "errorBarThreshold" ),
                                               0 ),
                   main = ifelse( isKeySet( LIST_arguments, "title" ),
                                  getValue( LIST_arguments, "title" ),
                                  NA ) )
}

if( STRING_function == "timeseries" )
{
  # check, if input is sane for this plot and get input files
  VEC_rmsdAll <- c( "snapshotsPerTimeInt", "timeUnit", "valueUnit", "valueName" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_rmsdAll,
                  VEC_allowedForAll ), LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_rmsdAll,
                   VEC_allowedForAll ),
                c( "<factor by which the values should be divided> (default: 1000)",
                   "<abbreviation of unit used for time-axis> (default: ns)",
                   "<abbreviation of unit used for response variable> (optional)",
                   "<name of response variable> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }    
  
  # plot
  MDplot::timeseries( MDplot::load_timeseries( VEC_files ),
                      names = VEC_dataNames,
                      snapshotsPerTimeInt = ifelse( isKeySet( LIST_arguments, "snapshotsPerTimeInt" ),
                                                    as.numeric( getValue( LIST_arguments, "snapshotsPerTimeInt" ) ),
                                                    1000 ),
                      timeUnit = ifelse( isKeySet( LIST_arguments, "timeUnit" ),
                                         getValue( LIST_arguments, "timeUnit" ),
                                         "ns" ),
                      valueUnit = ifelse( isKeySet( LIST_arguments, "valueUnit" ),
                                          getValue( LIST_arguments, "valueUnit" ),
                                          NA ),
                      main = ifelse( isKeySet( LIST_arguments, "title" ),
                                     getValue( LIST_arguments, "title" ),
                                     NA ) )
}



if( STRING_function == "clusters" )
{
  # check, if input is sane for this plot and get input files
  VEC_clustAll <- c( "clustersNumber" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_clustAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_clustAll,
                   VEC_allowedForAll ),
                c( "<number of clusters (sorted by population) to be shown in the plot> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  
  # load matrix
  MAT_input <- MDplot::load_clusters( VEC_files,
                                      names = VEC_dataNames )

  # plot
  MDplot::clusters( MAT_input,
                    clustersNumber = ifelse( isKeySet( LIST_arguments, "clustersNumber" ),
                                             getValue( LIST_arguments, "clustersNumber" ),
                                             NA ),
                    xlab = ifelse( is.na( VEC_axisNames ), "clusters", VEC_axisNames[ 1 ] ),
                    ylab = ifelse( is.na( VEC_axisNames ), "# configurations", VEC_axisNames[ 2 ] ),
                    main = ifelse( isKeySet( LIST_arguments, "title" ),
                                   getValue( LIST_arguments, "title" ),
                                   NA ) )
}



if( STRING_function == "clusters_ts" )
{
  # check, if input is sane for this plot and get input files
  VEC_clustTSAll <- c( "clustersNumber",
                       "selectTraj",
                       "selectTime",
                       "timeUnit",
                       "snapshotsPerTimeInt",
                       "lengths" )
  testRequired( c( VEC_requiredForAll ),
                LIST_arguments )
  testAllowed( c( VEC_clustTSAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_clustTSAll,
                   VEC_allowedForAll ),
                c( "<number of cluster printed> (optional)",
                   "<vector of selected trajectories, separated by ','> (optional)",
                   "<range of selected time, separated by ','> (optional)",
                   "<print time in units> (optional)",
                   "<number of snapshots per time unit (see above)> (optional)",
                   "<vector of trajectory-lengths, separated by ','>",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  
  # load matrix and set names
  LIST_input <- MDplot::load_clusters_ts( path = VEC_files,
                                          lengths = as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                                                  "lengths" ),
                                                                                  split = ",",
                                                                                  fixed = TRUE ) ) ),
                                          names = NA )
  if( isKeySet( LIST_arguments, "names" ) )
  {
    VEC_names <- unlist( strsplit( getValue( LIST_arguments, "names" ),
                                   ",",
                                   fixed = TRUE ) )
    if( length( VEC_names ) != length( LIST_input ) )
      stop( paste( "Error while assigning user specified trajectory ",
                   "names, since the numbers (",
                   length( VEC_names ),
                   ",",
                   length( LIST_input ),
                   ") do not match.",
                   sep = "" ) )
    for( i in 1:length( LIST_input ) )
      LIST_input[[ i ]][[ 1 ]] <- VEC_names[ i ]
  }
  
  # plot
  MDplot::clusters_ts( LIST_input,
                       clustersNumber = ifelse( isKeySet( LIST_arguments, "clustersNumber" ),
                                                as.numeric( getValue( LIST_arguments, "clustersNumber" ) ),
                                                NA ),
                       timeUnit = ifelse( isKeySet( LIST_arguments, "timeUnit" ),
                                          getValue( LIST_arguments, "timeUnit" ),
                                          NA ),
                       snapshotsPerTimeInt = ifelse( isKeySet( LIST_arguments, "snapshotsPerTimeInt" ),
                                                     as.numeric( getValue( LIST_arguments, "snapshotsPerTimeInt" ) ),
                                                     1000 ),
                       main = ifelse( isKeySet( LIST_arguments, "title" ),
                                      getValue( LIST_arguments, "title" ),
                                      NA ) )
}



if( STRING_function == "hbond" )
{
  # check, if input is sane for this plot and get input files
  VEC_hbAll <- c( "acceptorRange", "donorRange" )
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_hbAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_hbAll,
                   VEC_allowedForAll ),
                c( "<range of selected acceptors, separated by ','> (optional)",
                   "<range of selected donors, separated by ','> (optional)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  VEC_acceptorRange <- NA
  if( isKeySet( LIST_arguments, "acceptorRange" ) )
    VEC_acceptorRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                                 "acceptorRange" ),
                                             ",",
                                             fixed = TRUE ) ) )
  VEC_donorRange <- NA
  if( isKeySet( LIST_arguments, "donorRange" ) )
    VEC_donorRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                              "donorRange" ),
                                                    ",",
                                                    fixed = TRUE ) ) )
  
  # plot
  MDplot::hbond( MDplot::load_hbond( VEC_files ),
                 printLegend = BOOL_printLegend,
                 donorRange = VEC_donorRange,
                 acceptorRange = VEC_acceptorRange,
                 main = ifelse( isKeySet( LIST_arguments, "title" ),
                                getValue( LIST_arguments, "title" ),
                                NA ) )
}



if( STRING_function == "hbond_ts" )
{
  VEC_hbtsAll <- c( "acceptorRange", "donorRange", "plotOccurences",
                    "printNames", "namesToSingle", "timeUnit",
                    "snapshotsPerTimeInt", "hbondIndices", "printAtoms" )
  # check, if input is sane for this plot and get input files
  testRequired( VEC_requiredForAll, LIST_arguments )
  testAllowed( c( VEC_hbtsAll,
                  VEC_allowedForAll ),
               LIST_arguments )
  if( isKeySet( LIST_arguments, "help" )
      && getValue( LIST_arguments, "help" ) == "TRUE" )
  {
    print_help( STRING_function,
                c( VEC_hbtsAll,
                   VEC_allowedForAll ),
                c( "<range of acceptors, separated by ','> (optional)",
                   "<range of donors, separated by ','> (optional)",
                   "<boolean, setting the right subplot with occurences> (default: FALSE)",
                   "<boolean, setting whether hbonds should be plotted with names> (default: FALSE)",
                   "<boolean, setting wether the names are in single-letter code> (default: FALSE)",
                   "<abbreviation for time unit> (optional)",
                   "<number of snapshots per time unit (see above)> (default: 1000)",
                   "<range of hbonds defined by their indices> (optional)",
                   "<boolean, which sets the labels to contain the atom names if true> (default: FALSE)",
                   VEC_allowedForAllDesc ) )
    quit( save = "no", status = 0, runLast = TRUE )
  }
  VEC_files <- getFiles( getValue( LIST_arguments, "files" ) )
  for( i in 1:length( VEC_files ) )
  {
    if( !file.exists( VEC_files[ i ] ) )
      stop( paste( "Error in file checking: seemingly, file",
                   VEC_files[ i ], "does not exist." ) )
  }
  if( length( VEC_files ) != 2 )
    stop( paste( "Error in file checking: seemingly, the number of provided files is",
                 length( VEC_files ), "and not 2, as expected." ) )
  VEC_acceptorRange = NA
  if( isKeySet( LIST_arguments, "acceptorRange" ) )
    VEC_acceptorRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                                 "acceptorRange" ),
                                                       ",",
                                                       fixed = TRUE ) ) )
  VEC_donorRange <- NA
  if( isKeySet( LIST_arguments, "donorRange" ) )
    VEC_donorRange <- as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                              "donorRange" ),
                                                    ",",
                                                    fixed = TRUE ) ) )
  VEC_hbondIndices <- NA
  if( isKeySet( LIST_arguments, "hbondIndices" ) )
    VEC_hbondIndices <- as.numeric( unlist( strsplit( getValue( LIST_arguments,
                                                                "hbondIndices" ),
                                                      ",",
                                                      fixed = TRUE ) ) )
  
  # plot
  MDplot::hbond_ts( MDplot::load_hbond_ts( VEC_files[ 1 ] ),
                    MDplot::load_hbond( VEC_files[ 2 ] ),
                    printNames = ifelse( isKeySet( LIST_arguments, "printNames" ),
                                         TRUE,
                                         FALSE ),
                    plotOccurences = ifelse( isKeySet( LIST_arguments, "plotOccurences" ),
                                                       TRUE,
                                                       FALSE ),
                    acceptorRange = VEC_acceptorRange,
                    donorRange = VEC_donorRange,
                    namesToSingle = ifelse( isKeySet( LIST_arguments, "namesToSingle" ),
                                                      TRUE,
                                                      FALSE ),
                    printAtoms = ifelse( isKeySet( LIST_arguments, "printAtoms" ),
                                         TRUE,
                                         FALSE ),
                    timeUnit = ifelse( isKeySet( LIST_arguments, "timeUnit" ),
                                       getValue( LIST_arguments, "timeUnit" ),
                                       NA ),
                    snapshotsPerTimeInt = ifelse( isKeySet( LIST_arguments,
                                                            "snapshotsPerTimeInt" ),
                                                  as.numeric( getValue( LIST_arguments,
                                                                        "snapshotsPerTimeInt" ) ),
                                                  1000 ),
                    hbondIndices = list( VEC_hbondIndices ),
                    main = ifelse( isKeySet( LIST_arguments, "title" ),
                                   getValue( LIST_arguments, "title" ),
                                   NA ) )
}
#########

# end
dev.off()

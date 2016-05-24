#' Plot table level statistics, histograms, correlations and scatterplots in one go.
#' 
#' \code{showData} is the basic plotting function in the \code{toaster} package, designed to produce set of 
#' standard visualizations (see parameter \code{format}) in a single call. Depending on the \code{format} it 
#' is a wrapper to other functions or simple plotting function. It does all work in a single call by combining 
#' database round-trip (if necessary) and plotting functionality.
#' 
#' All formats support parameters \code{include} and \code{except} to include and exclude table columns respectively.
#' The \code{include} list guarantees that no columns outside of the list will be included in the results. 
#' The \code{excpet} list guarantees that its columns will not be included in the results.
#' 
#' Format \code{overview}: produce set of histograms - one for each statistic measure - across table columns. Thus,
#' it allows to compare averages, IQR, etc. across all or selected columns.
#' 
#' Format \code{boxplot}: produce boxplots for table columns. Boxplots can belong to the same plot or can be placed
#' inside facet each (see logical parameter \code{facet}).
#' 
#' Format \code{histogram}: produce histograms - one for each column - in a single plot or in facets (see logical 
#' parameter \code{facet}).
#' 
#' Format \code{corr}: produce correlation matrix of numeric columns.
#' 
#' Format \code{scatterplot}: produce scatterplots of sampled data.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param tableInfo pre-built summary of data to use (parameters \code{channel}, 
#'   \code{tableName}, \code{where} may not apply depending on \code{format}).
#'   See \code{\link{getTableSummary}}.
#' @param include a vector of column names to include. Output never contains attributes other than in the list.
#' @param except a vector of column names to exclude. Output never contains attributes from the list.
#' @param type what type of data to visualize: numerical (\code{"numeric"}), character (\code{"character"} or 
#'   date/time (\code{"temporal"})
#' @param format type of plot to use: \code{'overview'}, \code{'histogram'}, \code{'boxplot'}, \code{'corr'} for correlation 
#'   matrix or \code{'scatterplot'}
#' @param measures applies to format \code{'overview'} only. Use one or more of the following with \code{'numieric'} \code{type}:
#'   maximum,minimum,average,deviation,0%,10%,25%,50%,75%,90%,100%,IQR. Use one or more of the following with \code{'character'}
#'   \code{type}: distinct_count,not_null_count. By default all measures above are used per respeictive type.
#' @param title plot title
#' @param corrLabel column name to use to label correlation table: \code{'value'}, \code{'pair'}, or \code{'none'} (default)
#' @param digits number of digits to use in correlation table text (when displaying correlation coefficient value)
#' @param shape shape of correlation figure (default is 21)
#' @param shapeSizeRange correlation figure size range
#' @param facet Logical - if TRUE then divide plot into facets for each COLUMN (defualt is FALSE - no facets). 
#'   When set to TRUE and format is 'boxplot' scales defalut changes from 'fixed' to 'free'. Has no effect 
#'   when format is 'corr'.
#' @param numBins number of bins to use in histogram(s)
#' @param useIQR logical indicates use of IQR interval to compute cutoff lower and upper 
#'   bounds for values to be included in boxplot or histogram: \code{[Q1 - 1.5 * IQR, Q3 + 1.5 * IQR], IQR = Q3 - Q1}, 
#'   if FALSE then maximum and minimum are bounds (all values)
#' @param extraPoints vector contains names of extra points to add to boxplot lines. 
#' @param extraPointShape extra point shape (see 'Shape examples' in \link{aes_linetype_size_shape}).
#' @param sampleFraction sample fraction to use in the sampling of data for \code{'scatterplot'}
#' @param sampleSize if \code{sampleFraction} is not specified then size of sample must be specified 
#'   for \code{'scatterplot'}. 
#' @param pointColour name of column with values to colour points in \code{'scatterplot'}. 
#' @param facetName name(s) of the column(s) to use for faceting when \code{format} is \code{'scatterplot'}. 
#'   When single name then facet wrap kind of faceting is used. When two names then facet grid kind of 
#'   faceting is used. It overrides \code{facet} value in case of \code{'scatterplot'}. Must be part of 
#'   column list (e.g. \code{include}).
#' @param regressionLine logical if TRUE then adds regression line to scatterplot.
#' @param ncol Number of columns in facet wrap.
#' @param scales Are scales shared across all facets: \code{"fixed"} - all are the same, 
#'   \code{"free_x"} - vary across rows (x axis), \code{"free_y"} - vary across columns (Y axis) (default),
#'   \code{"free"} - both rows and columns (see in \code{facet_wrap} parameter \code{scales}. 
#'   Also see parameter \code{facet} for details on default values.)
#' @param coordFlip logical flipped cartesian coordinates so that horizontal becomes vertical, 
#'   and vertical, horizontal (see \link{coord_flip}).
#' @param paletteName palette name to use (run \code{display.brewer.all} to see available palettes).
#' @param baseSize base font size.
#' @param baseFamily base font family.
#' @param legendPosition legend position.
#' @param defaultTheme plot theme to use, default is \code{theme_bw}.
#' @param themeExtra any additional \code{ggplot2} theme attributes to add.
#' @param where SQL WHERE clause limiting data from the table (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param test logical: when applicable if TRUE show what would be done, only 
#'   (similar to parameter \code{test} in \link{RODBC} functions like \link{sqlQuery}
#'   and \link{sqlSave}). Doesn't apply when no sql expected to run, e.g. format
#'   is \code{'boxplot'}.
#' @return a ggplot object
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # get summaries to save time
#' pitchingInfo = getTableSummary(conn, 'pitching_enh')
#' battingInfo = getTableSummary(conn, 'batting_enh')
#' 
#' # Boxplots
#' # all numerical attributes
#' showData(conn, tableInfo=pitchingInfo, format='boxplot', 
#'          title='Boxplots of numeric columns')
#' # select certain attributes only
#' showData(conn, tableInfo=pitchingInfo, format='boxplot', 
#'          include=c('wp','whip', 'w', 'sv', 'sho', 'l', 'ktobb', 'ibb', 'hbp', 'fip', 
#'                    'era', 'cg', 'bk', 'baopp'), 
#'          useIQR=TRUE, title='Boxplots of Pitching Stats')
#' # exclude certain attributes
#' showData(conn, tableInfo=pitchingInfo, format='boxplot', 
#'          except=c('item_id','ingredient_item_id','facility_id','rownum','decadeid','yearid',
#'                   'bfp','ipouts'),
#'          useIQR=TRUE, title='Boxplots of Pitching Stats')
#' # flip coordinates
#' showData(conn, tableInfo=pitchingInfo, format='boxplot', 
#'          except=c('item_id','ingredient_item_id','facility_id','rownum','decadeid','yearid',
#'                   'bfp','ipouts'),
#'          useIQR=TRUE, coordFlip=TRUE, title='Boxplots of Pitching Stats')
#' 
#' # boxplot with facet (facet_wrap)
#' showData(conn, tableInfo=pitchingInfo, format='boxplot',
#'          include=c('bfp','er','h','ipouts','r','so'), facet=TRUE, scales='free',
#'          useIQR=TRUE, title='Boxplots Pitching Stats: bfp, er, h, ipouts, r, so')
#' 
#' # Correlation matrix
#' # on all numerical attributes
#' showData(conn, tableName='pitching_enh', tableInfo=pitchingInfo, 
#'          format='corr')
#' 
#' # correlation matrix on selected attributes
#' # with labeling by attribute pair name and
#' # controlling size of correlation bubbles
#' showData(conn, tableName='pitching', tableInfo=pitchingInfo, 
#'          include=c('era','h','hr','gs','g','sv'), 
#'          format='corr', corrLabel='pair', shapeSizeRange=c(5,25))
#'
#' # Histogram on all numeric attributes
#' showData(conn, tableName='pitching', tableInfo=pitchingInfo, include=c('hr'), 
#'          format='histogram')
#' 
#' # Overview is a histogram of statistical measures across attributes
#' showData(conn, tableName='pitching', tableInfo=pitchingInfo, 
#'          format='overview', type='numeric', scales="free_y")
#' 
#' # Scatterplots
#' # Scatterplot on pair of numerical attributes
#' # sample by size with 1d facet (see \code{\link{facet_wrap}})
#' showData(conn, 'pitching_enh', format='scatterplot', 
#'          include=c('so', 'er'), facetName="lgid", pointColour="lgid", 
#'          sampleSize=10000, regressionLine=TRUE,
#'          title="SO vs ER by League 1980-2000",
#'          where='yearid between 1980 and 2000')
#' 
#' # sample by fraction with 2d facet (see \code{\link{facet_grid}})
#' showData(conn, 'pitching_enh', format='scatterplot', 
#'          include=c('so','er'), facetName=c('lgid','decadeid'), pointColour="lgid",
#'          sampleFraction=0.1, regressionLine=TRUE,
#'          title="SO vs ER by League by Decade 1980 - 2012",
#'          where='yearid between 1980 and 2012')
#' }
showData <- function(channel = NULL, tableName = NULL, tableInfo = NULL, 
                     include = NULL, except = NULL, 
                     type = 'numeric', format = 'histgoram', measures = NULL,
                     title = paste("Table", toupper(tableName), format, "of", type, "columns"),
                     numBins = 30, 
                     useIQR = FALSE, extraPoints = NULL, extraPointShape = 15,
                     sampleFraction = NULL, sampleSize = NULL, pointColour = NULL,
                     facetName = NULL, regressionLine = FALSE,
                     corrLabel = 'none', digits = 2, 
                     shape = 21, shapeSizeRange = c(1,10),
                     facet = ifelse(format == 'overview', TRUE, FALSE), 
                     scales = ifelse(facet & format %in% c('boxplot','overview'),"free", "fixed"),
                     ncol = 4, coordFlip = FALSE, paletteName = "Set1", 
                     baseSize = 12, baseFamily = "sans",
                     legendPosition = "none",
                     defaultTheme=theme_tufte(base_size = baseSize, base_family = baseFamily), 
                     themeExtra = NULL, 
                     where = NULL, test = FALSE) {
  
  # match argument values
  type = match.arg(type, c('numeric','character','temporal'))
  format = match.arg(format, c('overview','histogram','boxplot','scatterplot', 'corr'))
  corrLabel = match.arg(corrLabel, c('none','value','pair'))
  
  # facetName needs to be included if missing
  if (!missing(include) & any(!(facetName %in% include))) {
    include = unique(append(include, facetName))
  }
  
  # colourPoint needs to be included if missing
  if (!missing(include) & !missing(pointColour)) { 
    if (!(pointColour %in% include)) {
      include = unique(append(include, pointColour))
    }
  }
  
  if (missing(tableName) && format %in% c('histogram','scatterplot','corr'))
    stop("Must provide table name when format one of ('histogram','scatterplot','corr')")
  
  if (missing(tableInfo) && test) {
    stop("Must provide tableInfo when test==TRUE.")
  }
  
  tableName = normalizeTableName(tableName)
      
  if (missing(tableInfo)) {
    summary = getTableSummary(channel, tableName, include=include, except=except, where=where)
  }else {
    summary = includeExcludeColumns(tableInfo, include, except)
  }
  
  # check that we have info for all required columns
  if (!is.null(include)) {
    if (!is.null(except)) {
      include = setdiff(include, except)
    }
    if (!all(include %in% summary$COLUMN_NAME)) {
      stop(paste("Not all specified columns are in the table summary (tableInfo):", include[!(include %in% summary$COLUMN_NAME)]))
    }
  }
  
  dataNum = getNumericColumns(summary, names.only=FALSE)
  dataChar = getCharacterColumns(summary, names.only=FALSE)
  dataTemp = getTemporalColumns(summary, names.only=FALSE)
  
  getPalette = getDiscretePaletteFactory(paletteName)
  
  if (format=='boxplot') {
    if (type != 'numeric') {
      warning("Automatically regressing to numerical types for format 'boxplot'")
    }
    
    # exclude columns with all NULLs
    dataNum = dataNum[dataNum$not_null_count > 0, ]
    
    checkNonEmpty(dataNum)
    
    # compute boxplot bounds based on IQR flag
    if (useIQR) {
      dataNum$lower_bound = apply(dataNum[, c('minimum', "25%", "IQR")], 1, FUN=function(x) max(x['minimum'], x['25%']-1.5*x['IQR']))
      dataNum$upper_bound = apply(dataNum[, c('maximum', "75%", "IQR")], 1, FUN=function(x) min(x['maximum'], x['75%']+1.5*x['IQR']))
    }else {
      dataNum$lower_bound = dataNum$minimum
      dataNum$upper_bound = dataNum$maximum
    }

    p = ggplot(dataNum, aes_string(x = 'COLUMN_NAME')) +
      geom_boxplot(aes_string(ymin = 'lower_bound', lower = '`25%`', middle = '`50%`', upper = '`75%`', ymax = 'upper_bound'),
                   stat="identity", position="dodge") +
      labs(title=title, x='Columns') 
    
    if (!missing(extraPoints))
      for (point in extraPoints)
        p = p + geom_point(aes_string(x = 'COLUMN_NAME', y = point), shape=extraPointShape)

  }
  else if (format=='corr') {
    if (type != 'numeric') {
      warning("Ignoring non-numeric types for format 'corr'")
    }
    
    checkNonEmpty(dataNum)
    
    corrmat = computeCorrelations(channel, tableName, tableInfo=tableInfo, 
                                  include=include, except=except, 
                                  where=where, test=test)
    if (test) {
      return (corrmat)
    } 
    
    corrmat = cbind(corrmat, valuePretty=prettyNum(corrmat[, 'value'], digits=digits) , none='')
    corrmat$value = abs(corrmat$value)
    corrLabelName = list('none', 'valuePretty', 'corr')[match(corrLabel, c('none','value','pair'))]
    p = createBubblechart(corrmat, "metric1", "metric2", "value", label=unlist(corrLabelName), fill="sign",
                          shape=shape, shapeSizeRange=shapeSizeRange, labelSize=5, labelVJust=0,
                          title=title, legendPosition=legendPosition,
                          defaultTheme=defaultTheme, 
                          themeExtra=theme(axis.title = element_blank()))
    
    # remove fill legend 
    if (legendPosition == "none") 
      p = p + guides(fill = FALSE)
    
    # force parameter facet to FALSE
    facet = FALSE
  }
  else if (format=='histogram') {
    all_hist = data.frame()
    
    if (type=='numeric') {
      for(columnName in dataNum$COLUMN_NAME) {
        hist = computeHistogram(channel, tableName, columnName, tableInfo=summary,
                                numbins=numBins, useIQR=useIQR, where=where, test=test)
        if (test) {
          return (hist)
        }
        
        if (nrow(hist) == 0) {
          warning(paste("Histogram for column '", columnName, "' is empty - skipping it." ))
        }else {
          hist = cbind(COLUMN_NAME=columnName, hist)
          all_hist= rbind(all_hist, hist)
        }
      }
      
      p = createHistogram(all_hist, x='bin', fill='COLUMN_NAME', facet='COLUMN_NAME', 
                          title=title, xlab="Bins", ylab="Frequency",
                          paletteValues = getPalette(length(unique(all_hist$COLUMN_NAME))),
                          legendPosition=legendPosition,
                          defaultTheme=defaultTheme, 
                          themeExtra=themeExtra)
      
      # remove fill legend 
      if (legendPosition == "none") 
        p = p + guides(fill = FALSE)
      
      # force parameter facet to FALSE
      facet = FALSE
      
    }
    else if (type=='character') {
      stop("Factor histograms are not supported by showData - use computeHistogram/createHistogram functions instead.")
    }
    else if (type=='temporal') {
      stop("Datetime histograms are not supported by showData - use computeHistogram/createHistogram functions instead.")
    }
  }
  else if (format=='scatterplot') {
    # Validations
    if (missing(include) | length(include)<2) {
      stop("Scatterplot format requires parameter 'include' define x and y coordiantes.")
    }
    if (!all(include[1:2] %in% dataNum$COLUMN_NAME)) {
      stop("Scatterplot format is valid for numerical data only.")
    }
    
    data = computeSample(channel, tableName, sampleFraction, sampleSize, include=summary$COLUMN_NAME,
                         where=where, stringsAsFactors=default.stringsAsFactors(),
                         test=test)
    if (test) {
      return (data)
    }
    
    # factor data for facet and colours
    if (!all(is.null(facetName), is.null(pointColour))) {
      data[union(facetName, pointColour)] = lapply(data[union(facetName, pointColour)], factor)
    }
    
    p = ggplot(data, aes_string(x=include[1], y=include[2])) +
      (if (missing(pointColour)) {
        geom_point()
      }else {
        geom_point(aes_string(colour=pointColour)) 
      }) +
      (if (!missing(pointColour)) {
        scale_colour_manual(values = getPalette(length(unique(data[,pointColour]))))
      }) +
      labs(title=title, x=include[1], y=include[2]) +
      theme(legend.position=legendPosition) 
    
    if (regressionLine) {
      p = p + geom_smooth(method="lm")
    }
    
    if (!missing(facetName) & length(facetName)>0) {
      if (length(facetName)==1) {
        p = p + facet_wrap(stats::as.formula(paste("~", facetName)), ncol=ncol, scales=scales)
      }else {
        p = p + facet_grid(stats::as.formula(paste(facetName[[1]],"~",facetName[[2]])), scales=scales)
      }
      
    }
    
    # remove colour legend 
    if (legendPosition == "none") 
      p = p + guides(colour = FALSE)
    else 
      p = p + guides(colour = guide_legend())
    
    facet = FALSE
  }
  else if (format=='overview') {
    if (type=='character' && missing(measures)) {
      measures = c('distinct_count', 'not_null_count', 'null_count')
    }
    else if (type=='numeric' && missing(measures)) {
      measures = c('distinct_count','not_null_count','null_count',
                   'maximum','minimum','average','deviation',
                   '25%','50%','75%','IQR')
    }
    
    if (type=='character') 
      data = dataChar
    else if (type=='numeric')
      data = dataNum
    
    overview = melt(data, id.vars='COLUMN_NAME', measure.vars=measures)
    
    p = ggplot(overview, aes_string(x='COLUMN_NAME')) +
      geom_bar(aes_string(y='value', fill='COLUMN_NAME'), stat="identity", position="dodge") +
      scale_fill_manual(values = getPalette(nrow(data))) +
      facet_wrap(~variable, ncol=1, scales=scales) +
      labs(title=title, x='Columns')
    
    # remove fill legend 
    if (legendPosition == "none") 
      p = p + guides(fill = FALSE)
    
    # force facet to FALSE
    facet = FALSE
  }
  
  if (facet) {
    p = p + facet_wrap(~COLUMN_NAME, ncol=ncol, scales=scales)
  }
  
  if (coordFlip) {
    p = p + coord_flip()
  }
  
  p = p +
    defaultTheme + themeExtra
  
  return(p)
}


# Checks if vector or data set contains at least one element or column.
# Note that for data frame it will NOT check if it has 1 or more rows.
# 
checkNonEmpty <- function(data) {
  if (length(data) == 0L || nrow(data) == 0L) stop("Nothing to show: check lists of columns, include/except, type and plot compatibilities.")
}



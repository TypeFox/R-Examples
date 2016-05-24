#' Get spatial map of the input dataset.
#' 
#' @param dataset A list containing different information, should be the result of reading netcdf file using
#' \code{loadNcdf}.
#' @param method A string showing different calculating method for the map. More information please refer to
#' details.
#' @param member A number showing which member is selected to get, if the dataset has a "member" dimension. Default
#' is NULL, if no member assigned, and there is a "member" in dimensions, the mean value of the members will be
#' taken.
#' @param ...  several arguments including x, y, title, catchment, point, output, name, info, scale, color, 
#' type in \code{?getSpatialMap_mat} for details.
#' @return A matrix representing the raster map is returned, and the map is plotted.
#' @details
#' There are following methods to be selected, 
#' "meanAnnual": annual rainfall of each year is plotted.  
#' "winter", "spring", "autumn", "summer": MEAN seasonal rainfall of each year is plotted.
#' Month(number 1 to 12): MEAN month rainfall of each year is plotted, e.g. MEAN march rainfall of each year.
#' "mean", "max", "min": mean daily, maximum daily, minimum daily precipitation.
#' @examples
#' 
#' 
#' \dontrun{
#' #gridData provided in the package is the result of \code {loadNcdf}
#' data(tgridData)
#' getSpatialMap(tgridData, method = 'meanAnnual')
#' getSpatialMap(tgridData, method = 'winter')
#' 
#' 
#' getSpatialMap(tgridData, method = 'winter', catchment = testCat)
#' 
#' file <- system.file("extdata", "point.txt", package = "hyfo")
#' point <- read.table(file, header = TRUE, sep = ',' )
#' getSpatialMap(tgridData, method = 'winter', catchment = testCat, point = point)
#' }
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
getSpatialMap <- function(dataset, method = NULL, member = 'mean', ...) {

  #check input dataset
  checkHyfo(dataset)
  
  #range of the dataset just loaded 
  lon <- dataset$xyCoords$x
  lat <- dataset$xyCoords$y
  startTime <- as.POSIXlt(dataset$Dates$start, tz = 'GMT')
  yearIndex <- startTime$year + 1900
  monthIndex <-startTime$mon + 1
  data <- dataset$Data
  
  # Dimension needs to be arranged. Make sure first and second dimension is lat and lon.
  data <- adjustDim(data, ref = c('lon', 'lat', 'time'))
  
  # Because in the following part, only 3 dimensions are allowed, so data has to be processed.
  if (member == 'mean' & any(attributes(data)$dimensions == 'member')) {
    dimIndex3 <- which(attributes(data)$dimensions != 'member')
    data <- apply(data, MARGIN = dimIndex3, FUN = mean, na.rm = TRUE)
    message('Mean value of the members are returned.')
    
  } else if (member != 'mean' & any(attributes(data)$dimensions == 'member')) {
    dimIndex3 <- which(attributes(data)$dimensions == 'member')
    data <- chooseDim(data, dimIndex3, member, drop = TRUE)
    
  } else if (member != 'mean' & !any(attributes(data)$dimensions == 'member')){
    stop('There is no member part in the dataset, but you choose one, check the input
         dataset or change your arguments.')
  }
  
  
  
  
  if (is.null(method)) {
    
    warning('You should shoose a method, unless input is a matrix directly to be plotted.')
    #in case the dataset is ready to plot and no need to calculate
    
  } else if (method == 'meanAnnual') { 
    #mean value of the annual precipitation over the period of the data 
    #time <- proc.time()
    if (length(unique(monthIndex)) < 12) {
      warning ('There are less than 12 months in a year, the results may be inaccurate.')
    }
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex,  method = 'annual')
    #newTime <- proc.time() - time
    title_d  <- 'Mean Annual Precipitation (mm / year)'
    
  } else if (method == 'winter') {
    #mean value of the seasonal precipitation, in this case, winter 
    #time <- proc.time()
    wm <- match(c(12, 1, 2), unique(monthIndex))
    if (length(which(!is.na(wm))) < 3) {
      stop ('Winter has less than 3 months, check data and try to calculate every month
  seperately or choose another season.')
    }
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex, monthIndex = monthIndex, 
                      method = 'winter')
    #newTime <- proc.time() - time
    title_d <- 'Mean Winter Precipitation (mm / winter)'
    
  } else if (method == 'spring') {
    wm <- match(c(3, 4, 5), unique(monthIndex))
    if (length(which(!is.na(wm))) < 3) {
      stop ('Spring has less than 3 months, check data and try to calculate every month
  seperately or choose another season.')
    }
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex, monthIndex = monthIndex, 
                      method = 'spring')    
    title_d <- 'Mean Spring Precipitation (mm / spring)'
    
  } else if (method == 'summer') {
    wm <- match(c(6, 7, 8), unique(monthIndex))
    if (length(which(!is.na(wm))) < 3) {
      stop ('Summer has less than 3 months, check data and try to calculate every month
  seperately or choose another season.')
    }
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex, monthIndex = monthIndex, 
                      method = 'summer')    
    title_d <- 'Mean Summer Precipitation (mm / summer)'
    
  } else if (method == 'autumn') {
    
    wm <- match(c(9, 10, 11), unique(monthIndex))
    if (length(which(!is.na(wm))) < 3) {
      stop ('Autumn has less than 3 months, check data and try to calculate every month
  seperately or choose another season.')
    }
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex, monthIndex = monthIndex, 
                      method = 'autumn')    
    title_d <- 'Mean Autumn Precipitation (mm / autumn)'
    
  } else if (method == 'mean') {
    
    #sum value of the dataset, this procedure is to get the mean value
    data_new <- apply(data, MARGIN = c(2, 1), FUN = mean, na.rm = TRUE)
    title_d <- 'Mean Daily Precipitation (mm / day)'
    
  } else if (method == 'max') {
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = suppressWarnings(max), na.rm = TRUE)
    title_d <- 'Max Daily Precipitation (mm / day)'
    
  } else if (method == 'min') {
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = suppressWarnings(min), na.rm = TRUE)
    title_d <- 'Min Daily Precipitation (mm / day)'
    
  } else if (is.numeric(method)) {
    
    data_new <- apply(data, MARGIN = c(2, 1), FUN = getMeanPreci, yearIndex = yearIndex, monthIndex = monthIndex, 
                      method = method)    
    title_d <- paste(month.abb[method], 'Precipitation (mm / month)', sep = ' ')
    
  } else {
    wrongMethod <- method
    stop(paste('no method called', wrongMethod))
  }
  # This is to give attributes to the matrix and better be melted in ggplot.
  colnames(data_new) <- round(lon, 2)
  rownames(data_new) <- round(lat, 2)
  
  # If ... also has a title argument, this will cause conflicts. so title has to be renamed as title_d
  # This has to be paid a lot of attention when use ... to pass arguments.
  output <- getSpatialMap_mat(matrix = data_new, title_d = title_d, ...)
  return(output)
}





#' Replot raster matrix
#' 
#' replot the matrix output from \code{getSpatialMap}, when \code{output = 'data'} or output is default
#' value.
#' 
#' @param matrix A matrix raster, should be the result of \code{getSpatialMap()}, output should be default
#' or 'data'
#' @param title_d A string showing the title of the plot, defaut is NULL.
#' @param catchment A catchment file geting from \code{shp2cat()} in the package, if a catchment is available for background.
#' @param point A dataframe, showing other information, e.g., location of the gauging stations. The 
#' the data.frame should be with columes "name, lon, lat, z, value".
#' @param output A string showing the type of the output, if \code{output = 'ggplot'}, the returned 
#' data can be used in ggplot and \code{getSpatialMap_comb()}; if \code{output = 'plot'}, the returned data is the plot containing all 
#' layers' information, and can be plot directly or used in grid.arrange; if not set, the raster matrix data
#' will be returned.
#' @param name If \code{output = 'ggplot'}, name has to be assigned to your output, in order to differentiate
#' different outputs in the later multiplot using \code{getSpatialMap_comb}.
#' @param info A boolean showing whether the information of the map, e.g., max, mean ..., default is FALSE.
#' @param scale A string showing the plot scale, 'identity' or 'sqrt'.
#' @param color Most of time you don't have to set this, but if you are not satisfied with the 
#' default color, you can set your own palette here. e.g., \code{color = c('red', 'blue')}, then
#' the value from lowest to highest, will have the color from red to blue. More info about color,
#' please check ?palette().
#' @param ... \code{title, x, y} showing the title and x and y axis of the plot. e.g. \code{title = 'aaa'}
#'default is about precipitation.
#' @return A matrix representing the raster map is returned, and the map is plotted.
#' @examples
#' 
#' \dontrun{
#' data(tgridData)# the result of \code{loadNcdf}
#' #the output type of has to be default or 'data'.
#' a1 <- getSpatialMap(tgridData, method = 'mean')
#' a2 <- getSpatialMap(tgridData, method = 'max')
#' a3 <- getSpatialMap(tgridData, method = 'winter')
#' a4 <- getSpatialMap(tgridData, method = 'summer')
#' #For example, if we want to investigate the difference between mean value and max.
#' 
#' a5 <- a2 - a1
#' getSpatialMap_mat(a4)
#' 
#' #Or to investigate the difference between winter value and summer value.
#' a6 <- a3 - a4
#' getSpatialMap_mat(a6)
#' 
#' }
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @import ggplot2 plyr maps maptools rgeos
#' @importFrom stats median
#' @importFrom reshape2 melt
#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' 
#' \item Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software,
#' 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
#' 
#' \item Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical
#' Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
#' 
#' \item Original S code by Richard A. Becker and Allan R. Wilks. R version by Ray Brownrigg. Enhancements
#' by Thomas P Minka <tpminka at media.mit.edu> (2015). maps: Draw Geographical Maps. R package version
#' 2.3-11. http://CRAN.R-project.org/package=maps
#' 
#' \item Roger Bivand and Nicholas Lewin-Koh (2015). maptools: Tools for Reading and Handling Spatial
#' Objects. R package version 0.8-36. http://CRAN.R-project.org/package=maptools
#' 
#' \item Roger Bivand and Colin Rundel (2015). rgeos: Interface to Geometry Engine - Open Source (GEOS). R
#' package version 0.3-11. http://CRAN.R-project.org/package=rgeos
#' 
#' }
#' 
#' 
#' 
#' 
#' 
getSpatialMap_mat <- function(matrix, title_d = NULL, catchment = NULL, point = NULL, output = 'data', 
                              name = NULL, info = FALSE, scale = 'identity', color = NULL, ...) {
  #check input
  checkWord <- c('lon', 'lat', 'z', 'value')
  if (is.null(attributes(matrix)$dimnames)) {
    stop('Input matrix is incorrect, check help to know how to get the matrix.')
  } else if (!is.null(catchment) & class(catchment) != "SpatialPolygonsDataFrame") {
    stop('Catchment format is incorrect, check help to get more details. ')
  } else if (!is.null(point) & any(is.na(match(checkWord, attributes(point)$names)))) {
    stop('point should be a dataframe with colnames "lon, lat, z, value".')
  }
  
  #ggplot
  #for the aes option in ggplot, it's independent from any other command through all ggplot, and aes() function
  #get data from the main dataset, in this case, data_ggplot. for other functions in ggplot, if it wants to use 
  #data from the main dataset as parameters, it has to use aes() function. if not, it has to use data available 
  #in the environment.
  #in other words, all the parameters in aes(), they have to come from the main dataset. Otherwise, just put them
  #outside aes() as normal parameters.
  
  if (info == TRUE) {  
    plotMax <- round(max(matrix, na.rm = TRUE), 2)
    plotMin <- round(min(matrix, na.rm = TRUE), 2)
    plotMean <- round(mean(matrix, na.rm = TRUE), 2)
    plotMedian <- round(median(matrix, na.rm = TRUE), 2)
    word <- paste('\n\n', paste('Max', '=', plotMax), ',', paste('Min', '=', plotMin), ',',
                  paste('Mean', '=', plotMean), ',', paste('Median', '=', plotMedian))
  } else {
    word <- NULL
  }
  
  x_word <- paste('Longitude', word)
  world_map <- map_data('world')
  
  # For some cases, matrix has to be reshaped, because it's too fat or too slim, to make
  # it shown on the map, the ratio is x : y is 4 : 3.
  matrix <- reshapeMatrix(matrix)
  
  
  # cannot remove NA, or the matrix shape will be changed.
  data_ggplot <- melt(matrix, na.rm = FALSE) 
  
  colnames(data_ggplot) <- c('lat', 'lon', 'value')
  theme_set(theme_bw())
  
  if (is.null(color)) color <- c('yellow', 'orange', 'red')
  # if (is.null(color)) color <- rev(rainbow(n = 20, end = 0.7))
  
  mainLayer <- with(data_ggplot, {
    
    ggplot(data = data_ggplot) +
    geom_tile(aes(x = lon, y = lat, fill = value)) +
    #scale_fill_discrete()+
    scale_fill_gradientn(colours = color, na.value = 'transparent') +#usually scale = 'sqrt'
                        #guide = guide_colorbar, colorbar and legend are not the same.
    guides(fill = guide_colourbar(title='Rainfall (mm)', barheight = rel(9), trans = scale)) +#usually scale = 'sqrt'
    geom_map(data = world_map, map = world_map, aes(map_id = region), fill = 'transparent', 
             color='black') +
    #    guides(fill = guide_colorbar(title='Rainfall (mm)', barheight = 15))+
    xlab(x_word) +
    ylab('Latitude') +
    ggtitle(title_d) +
    labs(empty = NULL, ...) +#in order to pass "...", arguments shouldn't be empty.
    theme(plot.title = element_text(size = rel(1.8), face = 'bold'),
          axis.title.x = element_text(size = rel(1.7)),
          axis.title.y = element_text(size = rel(1.7)),
          axis.text.x = element_text(size = rel(1.9)),
          axis.text.y = element_text(size = rel(1.9)),
          legend.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3)))
#    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim)
    
#   geom_rect(xmin=min(lon)+0.72*(max(lon)-min(lon)),
#             xmax=min(lon)+0.99*(max(lon)-min(lon)),
#             ymin=min(lat)+0.02*(max(lat)-min(lat)),
#             ymax=min(lat)+0.28*(max(lat)-min(lat)),
#             fill='white',colour='black')+
#   annotate('text', x = min(lon), y = min(lat), label=word, hjust = 0, vjust = -1)
  
  })
  
  printLayer <- mainLayer
  
  #catchment conversion
  if (is.null(catchment) == FALSE) {
    a <- catchment
    a@data$id <- rownames(a@data)
    b <- fortify(a, region = 'id')
    c <- join(b, a@data, by = 'id')
    catchmentLayer <- with(c, {
      geom_polygon(data = c, aes(long, lat, group = group), color = 'black', 
                                   fill = 'transparent')
    })
      
    
    printLayer <- printLayer + catchmentLayer
  }
  #plot point
  if (is.null(point) == FALSE) {
    pointLayer <- with(point, {
      geom_point(data = point, aes(x = lon, y = lat, size = value, colour = z),
                 guide = guide_legend(barheight = rel(3)))
        
        
    })
    
    printLayer <- printLayer + pointLayer
  }
  
  print(printLayer)
  
  if (output == 'ggplot') {
    if (is.null(name)) stop('"name" argument not found, 
                            If you choose "ggplot" as output, please assign a name.')
    data_ggplot$Name <- rep(name, dim(data_ggplot)[1])
    return (data_ggplot)
  } else if (output == 'plot') {
    return(printLayer)
  } else {
    return(matrix)
  }
}


#' Combine maps together
#' @param ... different maps generated by \code{getSpatialMap(, output = 'ggplot')}, see details.
#' @param nrow A number showing the number of rows.
#' @param list If input is a list containing different ggplot data, use \code{list = inputlist}.
#' @param x A string of x axis name.
#' @param y A string of y axis name.
#' @param title A string of the title.
#' @param output A boolean, if chosen TRUE, the output will be given.
#' @return A combined map.
#' @examples
#' 
#' 
#' \dontrun{
#' data(tgridData)# the result of \code{\link{loadNcdf}}
#' #The output should be 'ggplot'
#' a1 <- getSpatialMap(tgridData, method = 'summer', output = 'ggplot', name = 'a1')
#' a2 <- getSpatialMap(tgridData, method = 'winter', output = 'ggplot', name = 'a2')
#' a3 <- getSpatialMap(tgridData, method = 'mean', output = 'ggplot', name = 'a3')
#' a4 <- getSpatialMap(tgridData, method = 'max', output = 'ggplot', name = 'a4')
#' getSpatialMap_comb(a1, a2)
#' 
#' # or you can put them into a list.
#' getSpatialMap_comb(list = list(a1, a2), nrow = 2)
#' }
#' 
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @details
#' For \code{getSpatialMap_comb}, the maps to be compared should be with same size and resolution, 
#' in other words, they should be fully overlapped by each other.
#' 
#' If they have different resolutions, use \code{interpGridData{ecomsUDG.Raccess}} to interpolate.
#' 
#' @export
#' @import ggplot2 maps
#' @references 
#' 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' }
getSpatialMap_comb <- function(..., list = NULL, nrow = 1, x = '', y = '', title = '', 
                               output = FALSE) {
  
  
  if (!is.null(list)) {
    data_ggplot <- do.call('rbind', list)
  } else {
    maps <- list(...)
    checkBind(maps, 'rbind')
    data_ggplot <- do.call('rbind', maps)
  }
  
  if (!class(data_ggplot) == 'data.frame') {
    warning('Your input is probably a list, but you forget to add "list = " before it.
            Try again, or check help for more information.')
  } else if (is.null(data_ggplot$Name)) {
    stop('No "Name" column in the input data, check the arguments in getSpatialMap(), if 
         output = "ggplot" is assigned, more info please check ?getSpatialMap().')
  }
  
  data_ggplot$Name <- factor(data_ggplot$Name, levels = unique(data_ggplot$Name), ordered = TRUE)
  
#   lim <- getLim(data_ggplot$lon, data_ggplot$lat)
#   xlim <- lim[[1]]                                  
#   ylim <- lim[[2]]
  
  world_map <- map_data('world')
  theme_set(theme_bw())
  mainLayer <- with(data_ggplot, {
    ggplot(data = data_ggplot) + 
    geom_tile(aes(x = lon, y = lat, fill = value)) +
    #scale_fill_gradient(high = 'red', low = 'yellow')+
    scale_fill_gradientn(colours = c('yellow', 'orange', 'red'), na.value = 'transparent') +#usually scale = 'sqrt'
    geom_map(data = world_map, map = world_map, aes(map_id = region), fill = 'transparent', color = 'black') +
#    guides(fill = guide_colourbar(title='Rainfall (mm)', barheight = rel(9), trans = scale)) +#
    facet_wrap(~ Name, nrow = nrow) +
    theme(plot.title = element_text(size = rel(1.8), face = 'bold'),
          axis.title.x = element_text(size = rel(1.7)),
          axis.title.y = element_text(size = rel(1.7)),
          axis.text.x = element_text(size = rel(1.9)),
          axis.text.y = element_text(size = rel(1.9)),
          legend.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3))) +
    # no solultion for some very fat or very slim, in facet ggplot2, so, it's not buitiful.
    #coord_equal() +
    labs(x = x, y = y, title = title)
  })
  
  
  suppressWarnings(print(mainLayer))
  
  if (output == TRUE) return(data_ggplot)
}



reshapeMatrix <- function(matrix) {
  # This is for the map plot to keep the ratio x : y == 4:3
  # mainly used in map plot in ggplot2.
  
  
  # So the input matrix should be reshaped, add in some NA values, 
  # in order to be shown appropriately in ggplot.
  
  lon <- as.numeric(colnames(matrix))
  lat <- as.numeric(rownames(matrix))
  
  dx <- mean(diff(lon))
  dy <- mean(diff(lat))
  
  lx <- max(lon) - min(lon)
  ly <- max(lat) - min(lat)
  
  
  if (0.75 * lx < ly) {
    # In this case, x needs to be made longer
    
    xhalf <- 0.67 * ly
    xadd <- xhalf - lx / 2
    # calculate how many columns needs to be added.
    nxadd <- abs(round(xadd / dx))
    
    madd1 <- matrix(data = NA, nrow = length(lat), ncol = nxadd)
    madd2 <- madd1
    colnames(madd1) <- seq(to = min(lon) - dx, length = nxadd, by = dx)
    colnames(madd2) <- seq(from = max(lon) + dx, length = nxadd, by = dx)
    
    matrix_new <- cbind(madd1, matrix, madd2)  
    
    
  } else if (0.75 * lx > ly) {
    
    yhalf <- 0.38 * lx
    yadd <- yhalf - ly / 2
    nyadd <- abs(round(yadd / dy))
    
    madd1 <- matrix(data = NA, nrow = nyadd, ncol = length(lon))
    madd2 <- madd1  
      
    rownames(madd1) <- seq(to = max(lat) + dy, length = nyadd, by = -dy)
    rownames(madd2) <- seq(from = min(lat) - dx, length = nyadd, by = -dy)
    
    matrix_new <- rbind(madd1, matrix, madd2)
    
  } else {
    matrix_new <- matrix
  }
  
  return(matrix_new)
}

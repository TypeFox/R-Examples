#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, Februrary 2015                       #
#######################################################################################################
#' Gets met data from the redimentioned CFSR data set and outputs SWAT IO format meteorological input files
#' @param centroids - data.frame object or location of csv file with two columns, the first being lattitude and second being longidute in decimal degrees. These should be the centroids of the swat model subbasins and are the locations of the ouput met data.
#' @param outDir - Directory where ouput files will be saved.
#' @return returns cfsr met data in swat IO format
#' @examples
#' \dontrun{
#' centroids = data.frame(lat = 38, lon = 79)
#' outDir = "test"
#' getSWATcfsr(centroids=centroids, outDir=outDir)
#' }
#' @importFrom EcoHydRology get_cfsr_latlon
#' @importFrom utils read.csv write.table
#' @export


getSWATcfsr <- function(centroids, outDir = getwd()){

  tag1 <- '\n... Requesting Data from cfsr.bse.vt.edu/swat-cfsr-v02.pl ...\n'

  cat(gettext(tag1),"\n", sep <- "", file = stderr())


  #read in csv file of centorids
  setwd(outDir)


  if(class(centroids) == 'character'){
    cent <- read.csv(centroids)
  } else if(class(centroids) == 'data.frame'){
    cent <- centroids
  }

  # row number of centroids
  ln <- nrow(cent)

  ####################################format for tmax
  #make headers with lat lon and evelatiuon la

  cfsr <- list()
  tmpls <- list()
  pcpls <- list()
  hmdls <- list()
  slrls <- list()
  wndls <- list()

  for(i in 1:ln){

    cfsr[[i]] <- get_cfsr_latlon(cent[i,1], cent[i,2])

  }

  for(a in 1:ln){


    #tmax.agg
    date.full <- getSWATdates(cfsr[[a]]$DATE)

    ###Recombine for matrix

    tmax.format <- data.frame(date.full, cfsr[[a]]$TMX)

    colnames(tmax.format) <- c('date', 'tmax')


    tmax.format$tmax <- formatC(tmax.format$tmax, digits = 1, width = 5, flag = '0', format = 'f')



    ############################################################################################################
    ############################################################################################################
    ##format for tmin

    tmin.format <- data.frame(tmax.format$date, cfsr[[a]]$TMN)

    colnames(tmin.format) <- c('date', 'tmin')

    tmin.format$tmin <- formatC(tmin.format$tmin, digits = 1, width = 5, flag = '0', format ='f')

    colnames(tmin.format) <- c('date', 'tmin')

    ############################################################################################################

    ############################################################################################################
    ############################################################################################################
    ###Combine max and min into one temp file
    temp.format <- data.frame(tmax.format, tmin.format$tmin)

    colnames(temp.format) <- c('date', 'tmax', 'tmin')


    ############################################################################################################
    ############################################################################################################
    ##for rhum

    rhum.format <- data.frame(tmax.format$date, cfsr[[a]]$AVGRH)

    colnames(rhum.format) <- c('date', 'rhum')

    rhum.format$rhum <- formatC(rhum.format$rhum, digits = 3, width = 8, flag = '0', format ='f')

    colnames(rhum.format) <- c('date', 'rhum')

    ############################################################################################################
    ############################################################################################################
    ##for pcp
    pcp.format <- data.frame(tmax.format$date, cfsr[[a]]$PRECIP)

    colnames(pcp.format) <- c('date', 'pcp')

    pcp.format$pcp <- formatC(pcp.format$pcp, digits = 1, width = 5, flag = '0', format ='f')

    colnames(pcp.format) <- c('date', 'pcp')

    ############################################################################################################
    ############################################################################################################
    ##for wnd
    wnd.format <- data.frame(tmax.format$date, cfsr[[a]]$WIND)

    colnames(wnd.format) <- c('date', 'wnd')


    wnd.format$wnd<- formatC(wnd.format$wnd, digits = 3, width = 8, flag = '0', format = 'f')

    colnames(wnd.format) <- c('date', 'wnd')


    ############################################################################################################
    ############################################################################################################
    ##for slr
    slr.format <- data.frame(tmax.format$date, cfsr[[a]]$SOLAR)

    colnames(slr.format) <- c('date', 'slr')


    slr.format$slr <- formatC(slr.format$slr, digits = 3, width = 8, flag = '0', format ='f')

    colnames(slr.format) <- c('date', 'slr')


    ############################################################################################################
    ############################################################################################################
    #everthing into there own big ass lists
    tmpls[[a]] <- temp.format
    pcpls[[a]] <- pcp.format
    hmdls[[a]] <- rhum.format
    slrls[[a]] <- slr.format
    wndls[[a]] <- wnd.format

    #temp is different than the rest, but get the lists ready for c bind
    tmpls[[a]][,1] <- NULL
    tmpls[[a]]$full <- paste(tmpls[[a]][,1],tmpls[[a]][,2], sep = '')
    tmpls[[a]]$tmax <- NULL
    tmpls[[a]]$tmin <- NULL
    pcpls[[a]][,1] <- NULL
    hmdls[[a]][,1] <- NULL
    slrls[[a]][,1] <- NULL
    wndls[[a]][,1] <- NULL
    ############################################################################################################
    ############################################################################################################


  }


  #bind all the stations together
  alltmps <- do.call('cbind', tmpls)
  allpcps <- do.call('cbind', pcpls)
  allhmds <- do.call('cbind', hmdls)
  allslrs <- do.call('cbind', slrls)
  allwnds <- do.call('cbind', wndls)

  # get the full thing for printing
  tmpprint <- data.frame(tmax.format$date, alltmps)
  pcpprint <- data.frame(tmax.format$date, allpcps)
  hmdprint <- data.frame(tmax.format$date, allhmds)
  slrprint <- data.frame(tmax.format$date, allslrs)
  wndprint <- data.frame(tmax.format$date, allwnds)

  ############################################################################################################
  ############################################################################################################
  #make the headers for tmp and pcp

  Lati <- centroids[,1]

  Long <- centroids[,2]

  headtmp <- c('Station  tmp	Source cfsr', paste('Lati', paste(Lati, collapse = '   '), sep = '   '), paste('Long', paste(Long, collapse = '   '), sep = '   '), 'Elev')

  headpcp <- c('Station  pcp	Source cfsr', paste('Lati', paste(Lati, collapse = '   '), sep = '   '), paste('Long', paste(Long, collapse = '   '), sep = '   '), 'Elev')

  headhmd <- 'Relative Humidity % 	Source cfsr'

  headwnd <- 'Wind Speed m/s 	Source cfsr'

  headslr <- 'Solar Radiation MJ/m^2 	Source cfsr'


  ############################################################################################################
  ############################################################################################################
  #write tmp and pcp header files and append the data

  writeLines(headtmp, 'tmp1.tmp')

  writeLines(headpcp, 'pcp1.pcp')

  writeLines(headhmd, 'hmd.hmd')

  writeLines(headwnd, 'slr.slr')

  writeLines(headslr, 'wnd.wnd')

  write.table(tmpprint, 'tmp1.tmp', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  write.table(pcpprint, 'pcp1.pcp', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  write.table(hmdprint, 'hmd.hmd', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  write.table(slrprint, 'slr.slr', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  write.table(wndprint, 'wnd.wnd', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  ############################################################################################################
  ############################################################################################################

}

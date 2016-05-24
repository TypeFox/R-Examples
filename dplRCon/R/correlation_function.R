#' Performs correlations functions 
#' 
#' Uses the imported climate variables and tree ring data to produce seasonal correlation functions. Also, uses the bootstrapped chronologies to produce confidence intervals.
#' @param climate.anom.season.data climate data anomalies for seasons 
#		(including prior seasons), matrix: rows=year, col=season
#' @param site.chron.data site chronologies, matrix: rows=year, col=subset.1, subset.2,fullest
#' @param site.boot.data bootstrapped site chronologies, list: matrices for each subset	with: row=year, col=bootstrapped series
#' @param period.RF the period used to calculate response functions, vector: (start,end)
#' @param col.names.season  col.names.season<- list("SON_2", "DJF_2", "MAM_2", "JJA_2", "SON_1", "DJF_1", "MAM_1", "JJA_1", "SON", "DJF", "MAM", "JJA")
#' @param Climate.name name of the climate variable for which correlation functions are being calculated
#' @param Subset.name names given to each of the subsets. 
#' @note site.chron.data, and site.boot.data must be in the same order and confidence intervals plotted are for the 1st two subsets. 
#' @return \item{corr.site.1}{The correlations between the climate variable and the site chronology for the 1st subset.}
#' \item{corr.site.2}{The correlations between the climate variable and the site chronology for the 12st subset.}  
#' \item{percentile.ci.1}{The percentile confidence intervals for the 1st subset.} 
#' \item{percentile.ci.2}{The percentile confidence intervals for the 2st subset.} 
#' @return Other summary varibles \item{summary.ci.1}{}\item{summary.ci.2}{}\item{t.mean}{Test for correlation equal zero}\item{t.meanequal}{Test for correlations from the two subsets are equal}\item{percentile.ci}{}
#' @examples
#' \dontrun{ 
#' period.RF<-c(1900,1990)
#' col.names.season <- list("SON_2", "DJF_2", "MAM_2", "JJA_2", "SON_1", "DJF_1", "MAM_1","JJA_1", 
#'          "SON", "DJF", "MAM", "JJA")
#' ##  Full dataset
#' site.full <- site.chron(spline200.sub0.2000.n$sub.series.stand, aver.by.tree=F)
#' site.chron.data <- cbind(site.full$aver.site, site.full$aver.site)
#' site.boot.full <- ts(boot.full$boot.series.mean, start=tsp(site.full$aver.site)[1] )
#' site.boot.data<-list(site.boot.full, site.boot.full) 
#' 	
#' corr.SOI.full<-correlation.function(SOI.anom.season.data, site.chron.data,site.boot.data, 
#'       period.RF, col.names.season,
#'       Climate.name="SOI", Subset.name=c("0-20cm","20-200cm" ) )
#' corr.prec.full<-correlation.function(prec.anom.season.data, site.chron.data,site.boot.data, 
#'      period.RF, col.names.season, 
#'      Climate.name="SOI", Subset.name=c("0-20cm","20-200cm" ) )
#' corr.temp.full<-correlation.function(temp.anom.season.data, site.chron.data,site.boot.data, 
#'      period.RF, col.names.season, 
#'      Climate.name="SOI", Subset.name=c("0-20cm","20-200cm" ) )
#' 
#' ##	Near vs Far
#' site.0.20  <- site.chron(spline200.sub0.20.n$sub.series.stand, aver.by.tree=F)
#' site.20.200 <- site.chron(spline200.sub20.2000.n$sub.series.stand, aver.by.tree=F)
#' site.chron.data <- cbind(site.0.20$aver.site, site.20.200$aver.site)
#'
#'site.boot.0.20 <- ts(boot.0.20$boot.series.mean, start=tsp(site.0.20$aver.site)[1] )
#'site.boot.20.200 <- ts(boot.20.2000$boot.series.mean, start=tsp(site.20.200$aver.site)[1] )
#'site.boot.data<-list(site.boot.0.20, site.boot.20.200) 
#'
#'corr.SOI<-correlation.function(SOI.anom.season.data, site.chron.data, site.boot.data, 
#'    period.RF, col.names.season, 
#'    Climate.name="SOI",Subset.name=c("0-20cm","20-200cm" ) )
#'corr.prec<-correlation.function(prec.anom.season.data, site.chron.data,          site.boot.data, 
#'    period.RF, col.names.season, 
#'    Climate.name="SOI", Subset.name=c("0-20cm","20-200cm" ) )
#'corr.temp<-correlation.function(temp.anom.season.data, site.chron.data, site.boot.data, 
#'    period.RF, col.names.season, 
#'    Climate.name="SOI", Subset.name=c("0-20cm","20-200cm" ) )
#'}
#' @export

correlation.function<-function(climate.anom.season.data, site.chron.data, site.boot.data, period.RF,col.names.season, Climate.name,Subset.name){
  
  climate.season.wind <- as.matrix(window(climate.anom.season.data, 
                                          start=period.RF[1], end=period.RF[2]))
  
  site.chron.wind <- window(site.chron.data, 
                            start=period.RF[1], end=period.RF[2])
  
  
  ########## sub-functions: calculate correlations
  fun.corr.site.1 <- function(x){
    cor.test(as.numeric(x), site.chron.wind[, 1])$estimate
  }
  fun.corr.site.2 <- function(x){
    cor.test(as.numeric(x), site.chron.wind[, 2])$estimate
  }
  ##########
  
  corr.site.1<- apply(climate.season.wind, 2, fun.corr.site.1)
  
  corr.site.2 <- apply(climate.season.wind, 2, fun.corr.site.2)
  
  
  ########## sub-function:bootstrapped site chron correlations (confidence intervals)
  ### sub-sub-function:corrlations for bootstrapped series
  fun.corr.site.boot <- function(x, y){
    cor.test(x, y)$estimate
  }
  ###
  
  fun.boot.CI<-function(site.boot.data,  climate.season.wind, period.RF, i){
    site.boot.i<-site.boot.data[[i]]
    site.boot.wind.i<- window(site.boot.i[, 2:dim(site.boot.i)[2]],
                              start=period.RF[1], end=period.RF[2])
    
    corr.boot.i <- matrix(NA, dim(climate.season.wind)[2], 
                          dim(site.boot.wind.i)[2])
    
    for( r in 1:dim(site.boot.wind.i)[2]){
      for ( m in 1:dim(climate.season.wind)[2]){
        corr.boot.i[m, r] <- cor.test(as.numeric(site.boot.wind.i[, r]), 
                                      as.numeric(climate.season.wind[, m]))$estimate
      }
    }
    
    sort.corr.boot.i <- apply(corr.boot.i, 1, sort)
    
    percentile.ci.i <- cbind( sort.corr.boot.i [ 999*0.025, ], 
                              sort.corr.boot.i [ 999*0.975, ])
    summary.ci <- summary(sort.corr.boot.i)
    sd.ci <- sqrt(apply(sort.corr.boot.i, 2, var))
    return(list(percentile.ci.i=percentile.ci.i, summary.ci = summary.ci, sd.ci = sd.ci))
  }
  ###########
  percentile.ci<-NULL
  
  for( i in 1:2 ){
    percentile.ci [[ i]]<-fun.boot.CI(site.boot.data, climate.season.wind, 
                                      period.RF, i)
  }
  
  ############
  # t-test and p-values
  # ho: mean = 0
  t.mean <- rbind(corr.site.1/ percentile.ci [[ 1]]$sd.ci, corr.site.2/ percentile.ci [[ 2]]$sd.ci)
  # ho mean1=mean2
  t.meanequal <- (corr.site.1- corr.site.2)/(sqrt(percentile.ci [[ 1]]$sd.ci^2 + percentile.ci [[ 2]]$sd.ci^2 ))
  
  
  ############ plotting 
  plot(NULL, xlim=c(1, 12), ylim=c(-0.3, 0.3), xaxt = "n", xlab="", ylab="")
  axis(1,  at = 1:12,  las=2, labels =unlist( col.names.season))
  mtext("Pearson Correlation coeffients", side=2, line=3)
  abline(h=0)
  lines(corr.site.1,  col='blue')
  lines(corr.site.2,  col='red')
  arrows(1:12, percentile.ci[[1]]$percentile.ci.i[,1], 1:12,  percentile.ci[[1]]$percentile.ci.i[, 2], code = 3, length = 0.05, col = 4 , angle = 90)
  arrows(1:12, percentile.ci[[2]]$percentile.ci.i[, 1],1:12,  percentile.ci[[2]]$percentile.ci.i[, 2], code = 3, length = 0.05, col = 'red' , angle = 90)
  
  ############ outputs
  output<-list(corr.site.1=corr.site.1, corr.site.2=corr.site.2,  percentile.ci.1=percentile.ci[[1]]$percentile.ci.i, 
               percentile.ci.2=percentile.ci[[2]]$percentile.ci.i, summary.ci.1 = percentile.ci[[1]]$summary.ci, 
               summary.ci.2 = percentile.ci[[2]]$summary.ci, 
               t.mean = t.mean, t.meanequal = t.meanequal, percentile.ci[[1]]$sd.ci)
  return(output)
  
}

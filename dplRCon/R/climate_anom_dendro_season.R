#' Seasonal Data
#' 
#' Changes monthly climate data to seasonal climate data. 
#' @param climate.data This is a matrix, normally loaded into R in CSV form with the columns: year, jan, feb, mar,..,dec
#' @param year1 is the start of the period used to calculate climate anomalies.
#' @param year2 is the end of the period used to calculate climate anomalies. 
#' @param year.start is the start date of the climate data,
#' @param is.anomaly  accepts true or false. True will calculate the anomalies; departures from the mean from the period specified in year1 and year2
#' @details This is required to change the monthly data, such as the data used in this package, to the seasons; Sept-Oct-Nov, Dec-Jan-Feb, Mar-Apr-May, Jun-Jul-Aug. 
#' @return Return a matrix of time series with col names season for the seasons including lagged seasons; "SON_2","DJF_2","MAM_2","JJA_2","SON_1","DJF_1","MAM_1", "JJA_1","SON","DJF","MAM","JJA"; 
#' @author Maryann Pirie
#' @examples 
#' data(SOI)  # this is the Southern Oscillation Index data loaded with this package.
#' data(temperature) # this is the temperature data loaded with this package.
#' data(precipitation) # this is the precipitation data loaded with this package.
#' SOI.anom.season.data  <- climate.anom.dendro.season( SOI, 1933, 1992, 1876,
#'     is.anomaly="TRUE")
#' temp.anom.season.data	<- climate.anom.dendro.season( temperature, 1933, 1992, 1876, 
#'    is.anomaly="TRUE")
#' prec.anom.season.data	<- climate.anom.dendro.season( precipitation, 1933, 1992, 1876, 
#'    is.anomaly="TRUE")
#' @export
climate.anom.dendro.season<-function(climate.data,year1,year2,year.start,is.anomaly) {
  
  ############
  # calculate climate annomalies, bases on period c(year1,year2)
  colnames(climate.data)[1] <- "YEAR"
  x<-climate.data
  aver.month.year1.year2<-function(x,year1,year2){
    xx<-subset(x=x,subset=("YEAR">=year1))
    xxx<-subset(x=xx,subset=("YEAR"<=year2))
    colMeans(xxx)
  }
  year<-as.vector(climate.data$YEAR,"character")
  
  
  if (is.anomaly == "TRUE"){
    xxx<-climate.data[2:13] }else{
      aver.month<-aver.month.year1.year2(climate.data,year1,year2)
      month<-dimnames(x)[[2]][2:13]	
      xxx<-matrix("na",dim(x)[1],dim(x)[2]-1)
      rownames(xxx)<-year
      colnames(xxx)<-month	
      for( j in 1:12){
        xxx[,month[j]]<-apply(x,1,function(y){y[month[j]]-aver.month[month[j]]})	
      }
      
    }
  
  
  ########
  # anomialies dendro year Sept-Aug
  # ring 1901 is SON D 1901, and JF MAM JJA 1902
  
  Jan.Aug.dendro<-rbind(xxx[2:length(year),1:8],NA)
  rownames(Jan.Aug.dendro)<-year[1:(length(year))]
  Sep.Dec.dendro<-xxx[,9:12]
  dendro.anom<-cbind(Jan.Aug.dendro,Sep.Dec.dendro)
  
  
  ###########
  # Calculate seasonal means dendro year
  
  JJA<-matrix("na",length(year),1)
  rownames(JJA)<-year
  for( k in 1:length(year)){
    JJA[k]<-if(is.na(dendro.anom[year[k],"Jun"])==FALSE &
                 is.na(dendro.anom[year[k],"Jul"])==FALSE &
                 is.na(dendro.anom[year[k],"Aug"])==FALSE){
      mean(as.numeric(c(dendro.anom[year[k],"Jun"],
                        dendro.anom[year[k],"Jul"],
                        dendro.anom[year[k],"Aug"])))}else{"NA"}
  }
  
  MAM<-matrix("na",length(year),1)
  rownames(MAM)<-year
  for( k in 1:length(year)){
    MAM[k]<-if(is.na(dendro.anom[year[k],"Mar"])==FALSE &
                 is.na(dendro.anom[year[k],"Apr"])==FALSE &
                 is.na(dendro.anom[year[k],"May"])==FALSE){
      mean(as.numeric(c(dendro.anom[year[k],"Mar"],
                        dendro.anom[year[k],"Apr"],
                        dendro.anom[year[k],"May"])))}else{"NA"}
  }
  
  DJF<-matrix("na",length(year),1)
  rownames(DJF)<-year
  for( k in 1:length(year)){
    DJF[k]<-if(is.na(dendro.anom[year[k],"Dec"])==FALSE &
                 is.na(dendro.anom[year[k],"Jan"])==FALSE &
                 is.na(dendro.anom[year[k],"Feb"])==FALSE){
      mean(as.numeric(c(dendro.anom[year[k],"Dec"],
                        dendro.anom[year[k],"Jan"],
                        dendro.anom[year[k],"Feb"])))}else{"NA"}
  }
  
  SON<-matrix("na",length(year),1)
  rownames(SON)<-year
  for( k in 1:length(year)){
    SON[k]<-if(is.na(dendro.anom[year[k],"Sep"])==FALSE &
                 is.na(dendro.anom[year[k],"Oct"])==FALSE &
                 is.na(dendro.anom[year[k],"Nov"])==FALSE){
      mean(as.numeric(c(dendro.anom[year[k],"Sep"],
                        dendro.anom[year[k],"Oct"],
                        dendro.anom[year[k],"Nov"])))}else{"NA"}
  }
  
  ##########
  # season means lags
  
  JJA_1<-matrix(c("NA",JJA[1:(length(JJA)-1)]))
  rownames(JJA_1)<-year[1:(length(year))]
  MAM_1<-matrix(c("NA",MAM[1:(length(MAM)-1)]))
  rownames(MAM_1)<-year[1:(length(year))]
  DJF_1<-matrix(c("NA",DJF[1:(length(DJF)-1)]))
  rownames(DJF_1)<-year[1:(length(year))]
  SON_1<-matrix(c("NA",SON[1:(length(SON)-1)]))
  rownames(SON_1)<-year[1:(length(year))]
  
  JJA_2<-matrix(c("NA",JJA_1[1:(length(JJA_1)-1)]))
  rownames(JJA_2)<-year[1:(length(year))]
  MAM_2<-matrix(c("NA",MAM_1[1:(length(MAM_1)-1)]))
  rownames(MAM_2)<-year[1:(length(year))]
  DJF_2<-matrix(c("NA",DJF_1[1:(length(DJF_1)-1)]))
  rownames(DJF_2)<-year[1:(length(year))]
  SON_2<-matrix(c("NA",SON_1[1:(length(SON_1)-1)]))
  rownames(SON_2)<-year[1:(length(year))]
  
  climate.season<-cbind(SON_2,DJF_2,MAM_2,JJA_2,SON_1,DJF_1,MAM_1,JJA_1,SON,DJF,MAM,JJA)
  col.names.season<-list("SON_2","DJF_2","MAM_2","JJA_2","SON_1","DJF_1","MAM_1",
                         "JJA_1","SON","DJF","MAM","JJA")
  colnames(climate.season)<-(col.names.season)
  climate.season<- ts ( climate.season , start=year.start ) 
  climate.season
} 


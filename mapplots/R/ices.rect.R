ices.rect <-
function(rectangle){
  #this function converts a vector of ices rectangles into a dataframe with the lat
  #and long of the mid-points of the rectangles. It is can deal with rectangles 
  #like 48E2 that have been turned into numbers (48000 or 4.8e+3)
  if(is.factor(rectangle)) rectangle <- as.character(rectangle)
  lat <- function(r){
    split <- unlist(strsplit(as.character(r),''))
    if(split[1]%in%0:9 & split[2]%in%0:9 & split[3]%in%LETTERS & split[4]%in%0:9)
      ry <- as.numeric(substr(r,1,2)) else
        ry <- as.numeric(r)/10^floor(log10(as.numeric(r))-1)
    return((ry+71.5)/2)
    }
  lon <- function(r){
    split <- unlist(strsplit(as.character(r),''))
    if(split[1]%in%0:9 & split[2]%in%0:9 & split[3]%in%LETTERS & split[4]%in%0:9)
      rx <- as.numeric(paste(match(split[3],LETTERS),split[4],sep='')) else
        rx <- floor(log10(as.numeric(r))-1)+50
    return(rx-59.5)
    }
  data.frame(lon=unlist(lapply(rectangle,lon)),lat=unlist(lapply(rectangle,lat)))
}

ices.rect2 <-
function(lon,lat){
  x <- floor(lon+60)+1000
  y <- floor(lat*2)-71+100
  num1<- substr(y,2,3)
  lett <- LETTERS[as.numeric(substr(x,2,3))]
  num2 <- substr(x,4,4)
  paste(num1,lett,num2,sep='')
  }


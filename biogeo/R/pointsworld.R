pointsworld <-
function(world,dat,x="x",y="y",ext=c(-180,180,-90,90)){
  # pointsworld(fold,dat,x,y,ext=c(-180,180,-90,90))
  # plots points in object dat on a map of the world
  # tests for accidental transpose of latitude and longitude
  # colour points red that are in the sea 
  # x and y must be coordinates in decimal degrees, NA values are allowed
  # set extent see
  
  stopifnot(nrow(dat)>0) # nrow of dat must be >0
  nn<-names(world@data) # countries field name must be "NAME"
  fc<-match("NAME",nn) # which column
  if (is.na(fc)) stop("The field name for countries in world must be: 'NAME'")
  fieldsmissing(dat,fields=c("ID"))
  if(length(dat$ID)!=length(unique(dat$ID))) stop("ID contains duplicate numbers")
  nr<-nrow(dat)  # number of rows in dat
  cn<-names(dat) # names of columns
  f1<-match(x,cn) # which column is the x-coordinate
  f2<-match(y,cn) # which column is the y-coordinate
  xn<-coord2numeric(dat[,f1])
  yn<-coord2numeric(dat[,f2])
  xy<-cbind(xn,yn) # combine
  
  ex<-getextent(xn,yn,ext) # beyond from this function can be used below
  
  # impossible coordinate values
  x1<-(abs(xn)>180)*1
  y1<-(abs(yn)>90)*1
  
  if(any(x1+y1>0)){
    ff<-which(x1+y1>0)
    er<-str_c(dat$ID[ff],sep="",collapse=";")
    msg<-paste("coordinate values impossible see ID:",as.character(er))
    stop(msg)
  } 
    
  g<-SpatialPoints(xy) # convert points into SpatialPoints format
  plot(world,border="gray",xlim=ex$xlm,ylim=ex$ylm) # plot map of world and set limits to bbox of points
  world@proj4string=CRS(as.character(NA)) #
  s1<-over(g,world,returnList=FALSE) # overlay the points on the world map
  country_ext<-s1$NAME # country names 
  f<-which(is.na(country_ext)) # which points are in the sea (points that do not have a country name)
  f2<-which(!is.na(country_ext))
  dat1<-data.frame(dat,country_ext) # join country to dat
  points(xn[f2],yn[f2],pch=20,col="blue") # plot the points on the map
  points(xn[f],yn[f],pch=20,col="red") # points in sea in red
  #points(xn[f1],yn[f1],pch=22,col="black")
  return(dat1)
}

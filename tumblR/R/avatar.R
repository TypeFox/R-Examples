avatar <-
function(base_hostname=NA,size=64){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
   
  size_range<-c(16,24,30,40,48,64,96,128,512)
    
  if(!(size %in% size_range))
    stop("Avaliable values for size are: 16,24,30,40,48,64,96,128,512")
  
  url<-GET(paste("http://api.tumblr.com/v2/blog/",base_hostname,"/avatar/",as.character(size),sep=""))
  
  image_url<-url[[1]]
  
  return(image_url)
}

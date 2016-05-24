### TO CREATE THE mixsep.env ON START UP

.onLoad <- function(...){
  if(!exists("mixsep.env"))
    mixsep.env <<- new.env(hash=TRUE,size=NA)
}

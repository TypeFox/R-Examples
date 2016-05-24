default_memes <- 
function(...){
    out <- imgurGET('memegen/defaults', ...)
    lapply(out, `class<-`, 'imgur_image')
}
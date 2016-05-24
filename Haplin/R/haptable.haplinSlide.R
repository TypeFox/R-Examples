haptable.haplinSlide <- function(object){

.ut <- f.haptable.list(object)
names(.ut)[names(.ut) == "element"] <- "window"

return(.ut)
}

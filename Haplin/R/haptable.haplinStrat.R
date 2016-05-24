haptable.haplinStrat <- function(object){

.ut <- f.haptable.list(object)
names(.ut)[names(.ut) == "element"] <- "stratum"

return(.ut)
}

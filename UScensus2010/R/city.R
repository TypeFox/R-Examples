city<-function (name = NULL, state = NULL, statefips = FALSE, sp.object = NULL, proj = NULL) 
{
    city.aux <- function(name = NULL, state = NULL, statefips = FALSE, 
        sp.object = NULL, proj = NULL) {
        state <- check.state(state, statefips)
        if (is.null(state)) {
            stop("Not a State! \n")
        }
        if (is.null(sp.object) == FALSE) {
            if (class(sp.object)[1] != "SpatialPolygonsDataFrame") {
                stop("Not a SpatialPolygonsDataFrame object! \n")
            }
            assign(paste(state, ".cdp10", sep = ""), sp.object)
        }
        else {
            data(list = paste(state, ".cdp10", sep = ""), envir = parent.frame())
        }
        temp.cdp <- get(paste(state, ".cdp10", sep = ""))
        
        ####Need to fix this permantly 
        if(statefips){
        temp <- which(tolower(temp.cdp$place) %in%name == TRUE)
        }else{
        temp <- which(tolower(temp.cdp$name) %in% tolower(name) == TRUE)
        }
        
        if (length(temp) == 0) {
            stop(paste(name, "is not in this SpatialPolygonsDataFrame object \n or this city does not exist in this state!"))
        }
        out <- temp.cdp[temp, ]
        if (is.null(proj) == FALSE) {
            require(rgdal)
            out <- spTransform(out, proj)
        }
        out
    }
    out <- city.aux(name, state, statefips, sp.object, proj)
}
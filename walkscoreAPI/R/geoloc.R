geoloc <- function(address,apikey){
    place <- address
    place <- gsub(" ","+",place)

    URL <- paste("http://maps.google.com/maps/geo?q=",place,"&output=json&oe=utf8\
            &sensor=true_or_false&key=",apikey,sep="")

    X <- character(0)
    X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))

    coord <- grep("\"coordinates\":",X)
    coord <- strsplit(X[coord],"\"coordinates\"")
    coord <- coord[[1]][2]
    coord <- gsub(": [ ","",coord,fixed=TRUE)
    coord <- strsplit(coord[[1]],",")

    long <- as.numeric(coord[[1]][1])
    lat <- as.numeric(coord[[1]][2])

    acc <- grep("Accuracy",X)
    acc <- strsplit(X[acc]," : ")
    acc <- gsub(",","",acc[[1]][2])
    acc <- as.numeric(acc)

    loc <- grep("LocalityName", X)
    loc <- strsplit(X[loc], " : ")
    loc <- gsub("\"","",loc[[1]][2],fixed=TRUE)
    loc <- gsub(",","",loc)
     
    admin <- grep("AdministrativeAreaName",X)
    admin <- strsplit(X[admin], " : ")
    admin <- gsub("\"","",admin[[1]][2],fixed=TRUE)
    admin <- gsub(",","",admin)

    country <- grep("CountryName",X)
    country <- strsplit(X[country], " : ")
    country <- gsub("\"","",country[[1]][2],fixed=TRUE)
    country <- gsub(",","",country)

    object <- list()
    class(object) <- "GoogleGeoloc"
    object$coordinates <- c(long,lat)
    object$accuracy <- acc
    object$city <- loc
    object$state <- admin
    object$country <- country
    return(object)

}

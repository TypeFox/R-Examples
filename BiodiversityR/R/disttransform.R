`disttransform` <-
function(x, method="hellinger") {
    x <- as.matrix(x)
    METHODS <- c("hellinger", "chord", "profiles", "chi.square", "log", "square", "pa",
        "Braun.Blanquet", "Domin", "Hult", "Hill", "fix", "coverscale.log")
    method <- match.arg(method,METHODS)
    switch(method, hellinger = {
        x <- decostand(x,"hellinger")
    }, profiles = {
        x <- decostand(x,"total")
    }, chord = {
        x2 <- x^2
        rowtot <- apply(x2,1,sum)
        for (i in 1:length(rowtot)) {if (rowtot[i]==0) {rowtot[i] <- 1}}
        rowtot <- rowtot^0.5
        x <- x/rowtot
    }, chi.square = {
        x <- decostand(x, "chi.square")
    }, log = {
        x <- log(x+1)
    }, square = {
        x <- x^0.5
    }, pa = {
        x <- decostand(x, "pa")
#
    }, Braun.Blanquet = {
        x <- coverscale(x, "Braun.Blanquet")
    }, Domin = {
        x <- coverscale(x, "Domin")
    }, Hult = {
        x <- coverscale(x, "Hult")
    }, Hill = {
        x <- coverscale(x, "Hill")
    }, fix = {
        x <- coverscale(x, "fix")
    }, coverscale.log = {
        x <- coverscale(x, "log")
    })
#
    for (i in 1:ncol(x)) {x[,i] <- as.numeric(x[,i])}
    x <- data.frame(x)
    return(x)
}


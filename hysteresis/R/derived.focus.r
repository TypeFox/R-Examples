derived.focus <- function (semi.major,semi.minor,rote.rad) {
focus.distance <- sqrt(semi.major^2-semi.minor^2)
focus.x <- cos(rote.rad)*focus.distance
focus.y <- sin(rote.rad)*focus.distance
eccentricity <- sqrt((semi.major^2-semi.minor^2)/semi.major^2)
z <- c(focus.x,focus.y, eccentricity)
z
}

### Hard-coded Google+ API URLs

base.url <- "https://www.googleapis.com/plus/v1/"

start.people <- "people/"
close.page1 <- "/activities/public?maxResults=100"
close.page2 <- "&key="

close.people <- "?key="

gp <- new.env(parent=emptyenv())
assign("apikey", NULL, envir=gp)

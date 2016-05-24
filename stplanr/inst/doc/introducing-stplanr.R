## ---- include=FALSE------------------------------------------------------
library(stplanr)
library(sp) # needed for geographical objects
library(leaflet) # needed for plotting data

## ---- message=FALSE------------------------------------------------------
library(stplanr)

## ------------------------------------------------------------------------
data("flow", package = "stplanr")
head(flow[c(1:3, 12)])

## ------------------------------------------------------------------------
data("cents", package = "stplanr")
as.data.frame(cents[1:3,-c(3,4)])

## ------------------------------------------------------------------------
l <- od2line(flow = flow, zones = cents)

# remove lines with no length
l <- l[!l$Area.of.residence == l$Area.of.workplace,]

## ------------------------------------------------------------------------
plot(l, lwd = l$All / 10)

## ------------------------------------------------------------------------
# test the connection
x = try(download.file("http://www.google.com", "test.html"), silent = T) 
offline <- grepl("Error", x)
if(!offline) file.remove("test.html")

data("routes_fast")
data("routes_slow")

plot(l)
lines(routes_fast, col = "red")

## ---- eval=FALSE---------------------------------------------------------
#  # ggmap dependency may be down
#  if(!offline & !Sys.getenv("GRAPHHOPPER") == ""){
#    ny2oaxaca1 <- route_graphhopper("New York", "Oaxaca", vehicle = "bike")
#    ny2oaxaca2 <- route_graphhopper("New York", "Oaxaca", vehicle = "car")
#    plot(ny2oaxaca1)
#    plot(ny2oaxaca2, add = T, col = "red")
#  
#    ny2oaxaca1@data
#    ny2oaxaca2@data
#  }

## ------------------------------------------------------------------------
lgb <- spTransform(l, CRSobj = CRS("+init=epsg:27700"))
l$d_euclidean <- rgeos::gLength(lgb, byid = T)
l$d_fastroute <- routes_fast@data$length
plot(l$d_euclidean, l$d_fastroute,
  xlab = "Euclidean distance", ylab = "Route distance")
abline(a = 0, b = 1)
abline(a = 0, b = 1.2, col = "green")
abline(a = 0, b = 1.5, col = "red")

## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
routes_fast$All <- l$All
rnet <- overline(routes_fast, "All", fun = sum)
w <- rnet$All / mean(rnet$All) * 3
if(offline){
  print("Figure not shown: requires internet connection")
} else {
  n <- 30
  bb <- bbox(routes_fast[n,])
  leaflet() %>% addTiles() %>%
    addPolylines(data = rnet, color = "red", weight = w) 
}

## ------------------------------------------------------------------------
l$pwalk <- l$On.foot / l$All
plot(l$d_euclidean, l$pwalk, cex = l$All / 50,
  xlab = "Euclidean distance (m)", ylab = "Proportion of trips by foot")

## ------------------------------------------------------------------------
lm1 <- lm(pwalk ~ d_euclidean, data = l@data, weights = All)
lm2 <- lm(pwalk ~ d_fastroute, data = l@data, weights = All)
lm3 <- glm(pwalk ~ d_fastroute + I(d_fastroute^0.5), data = l@data, weights = All, family = quasipoisson(link = "log"))

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  summary(lm1)
#  summary(lm2)
#  summary(lm3)

## ------------------------------------------------------------------------
plot(l$d_euclidean, l$pwalk, cex = l$All / 50,
  xlab = "Euclidean distance (m)", ylab = "Proportion of trips by foot")
l2 <- data.frame(d_euclidean = 1:5000, d_fastroute = 1:5000)
lm1p <- predict(lm1, l2)
lm2p <- predict(lm2, l2)
lm3p <- predict(lm3, l2)
lines(l2$d_euclidean, lm1p)
lines(l2$d_euclidean, exp(lm2p), col = "green")
lines(l2$d_euclidean, exp(lm3p), col = "red")


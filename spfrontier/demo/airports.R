library(spdep)
library(pastecs)
library(car)
library(maptools)
library(maps)
library(ggmap)
library(corpcor)

rm(list=ls())
data(airports, envir = environment())


airports$Country<-as.factor(airports$Country)
airports$Routes <- airports$RoutesDeparture+airports$RoutesArrival
vars <- c("ICAO", "Country", "AirportName", "longitude", "latitude","Year", "PAX", "ATM", "Cargo", "Population100km", "Population200km", "Island", 
          "GDPpc","RunwayCount","CheckinCount","GateCount","ParkingSpaces", "RoutesDeparture", "RoutesArrival", "Routes")
airportsEU <- airports[vars]

airports2011 <- subset(airportsEU, Year==2011)
airports2011.complete <- with(airports2011, airports2011[complete.cases(PAX,ATM,Cargo,Population100km,Population200km,Island, GDPpc, RoutesDeparture,RoutesArrival),])
nrow(airports2011.complete)
stat.desc(airports2011.complete)

data <- airports2011.complete
nrow(data)
W <- constructW(cbind(data$longitude, data$latitude),data$ICAO)

formula <- log(PAX) ~ log(Population100km) + log(Routes) + log(GDPpc)


model000 <- spfrontier( formula , data=data)
summary(model000 )

model100 <- spfrontier(formula , data=data, W_y=W)
summary(model100 )

model010 <- spfrontier(formula , data=data, W_v=W, logging="debug")
summary(model010)

model001 <- spfrontier(formula , data=data, W_u=W, logging="debug")
summary(model001)
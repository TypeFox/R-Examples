context("TimeSeries")

seed <- 1.0

set.seed(seed)

offset <- runif(1) #cache braking offset. used in all tests

print(cat("Random seed is ",seed,"; Cache braking offset is",offset));

test_that("fcTimeSeriesYearly timeseries correct length for single point", {
        a <- fcTimeSeriesYearly(variable="airt",latitude=76.0+offset, longitude=57.7,firstYear=1950,lastYear=2000)
        expect_equal(ncol(a$values), 51,label=paste("wrong time series length. expected 51 but got",ncol(a$values)))
})

test_that("fcTimeSeriesDaily timeseries correct length for single point", {
        a<-fcTimeSeriesDaily(variable="airt",latitude=76.1+offset, longitude=57.7,firstYear=1950,lastYear=2000)
        expect_equal(ncol(a$values), 365,label=paste("wrong time series length. expected 365 but got",ncol(a$values)))
})

test_that("fcTimeSeriesHourly timeseries correct length for single point.", {
        a<-fcTimeSeriesHourly(variable="airt",latitude=75.1+offset, longitude=57.7,firstYear=1950,lastYear=2000,startHour=0,stopHour=24)
        expect_equal(ncol(a$values), 24,label=paste("wrong time series length. expected 24 but got",ncol(a$values)))
})

test_that("fcTimeSeriesHourly timeseries correct length for single point. timestamp set", {
        a<-fcTimeSeriesHourly(variable="airt",latitude=75.1+offset, longitude=57.7,firstYear=1950,lastYear=2000,startHour=0,stopHour=24,reproduceFor="2015-10-01")
        expect_equal(ncol(a$values), 24,label=paste("wrong time series length. expected 24 but got",ncol(a$values)))
})

test_that("fcTimeSeriesDaily timeseries correct length for a point set. timestamp set", {
        data(quakes) #the only built-in dataset with locations which I've found. The Fiji earthquakes
        a<-fcTimeSeriesDaily( #fetching day-to-day temperature variations at earthquake locations
                "airt", 
                quakes$lat+offset, quakes$long,
                firstYear=1981, lastYear=2000 #averaging across 20 years
        )
        expect_equal(ncol(a$values), 365,label=paste("wrong time series length. expected 365 but got",ncol(a$values)))
        expect_equal(nrow(a$values), length(quakes$lat),label=paste("wrong time series count expected",length(quakes$lat),"but got",ncol(a$values)))
})

test_that("fcTimeSeriesHourly with overriden url", {
        ts <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0+offset, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fc.dgrechka.net/',
                dataSets="CRU CL 2.0")
})

context("Grid")

test_that("fcGrid succeeds", {
        a<- fcGrid("airt",40+offset,80+offset,10,10,200,10,firstYear=1950,lastYear=2000,firstDay=1,lastDay=31)
})

test_that("fcGrid with long timestamp succeeds", {
        a<- fcGrid("airt",40+offset,80+offset,10,10,200,10,firstYear=1950,lastYear=2000,firstDay=1,lastDay=31,reproduceFor="2015-10-01 05:00:00")
})

test_that("fcGrid with short timestamp succeeds", {
        a<- fcGrid("airt",40+offset,80+offset,10,10,200,10,firstYear=1950,lastYear=2000,firstDay=1,lastDay=31,reproduceFor="2015-10-01")
})

test_that("fcGrid succeeds with dataset specified", {
        a<- fcGrid("airt",40+offset,80+offset,10,10,200,10,firstYear=1950,lastYear=2000,firstDay=1,lastDay=31,dataSets=c("NCEP/NCAR Reanalysis 1 (regular grid)"))
})

test_that("fcGrid succeeds with two datasets specified", {
        a<- fcGrid("airt",40+offset,80+offset,10,10,200,10,firstYear=1950,lastYear=2000,firstDay=1,lastDay=31,dataSets=c("NCEP/NCAR Reanalysis 1 (regular grid)","CRU CL 2.0"))
})


context("Fc paper figure related fetches")

test_that("original Fc Paper figure 1 fetch succeeds", {
        ts <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                reproduceFor='2015-05-27')
});


test_that("original Fc Paper figure 2 fetch succeeds", {
        africaJulyTemp <- fcGrid(variable="airt",
                                 latitudeFrom=-35, latitudeTo=35, latitudeBy=1,
                                 longitudeFrom=-20,longitudeTo=60,longitudeBy=1,
                                 firstDay=182,lastDay=212, #July
                                 firstYear=1950,lastYear=2000,
                                 url='http://fetchclim.cloudapp.net/',
                                 reproduceFor='2015-05-27') 
});


test_that("original Fc Paper figure 3a fetch succeeds", {
        ts <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet="GHCNv2",
                reproduceFor='2015-05-27')
});

test_that("original Fc Paper figure 3b fetch succeeds", {
        ts2 <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet ="NCEP/NCAR Reanalysis 1 (regular grid)",
                reproduceFor='2015-05-27')    
});
test_that("original Fc Paper figure 3c fetch succeeds", {
        ts3 <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet ="CESM1-BGC airt",
                reproduceFor='2015-05-27')
});

test_that("Fc Paper figure 1 fetch succeeds (random shifted)", {
        ts <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0+offset, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                reproduceFor='2015-05-27')
});


test_that("Fc Paper figure 2 fetch succeeds (random shifted)", {
        africaJulyTemp <- fcGrid(variable="airt",
                                 latitudeFrom=-35+offset, latitudeTo=35+offset, latitudeBy=1,
                                 longitudeFrom=-20,longitudeTo=60,longitudeBy=1,
                                 firstDay=182,lastDay=212, #July
                                 firstYear=1950,lastYear=2000,
                                 url='http://fetchclim.cloudapp.net/',
                                 reproduceFor='2015-05-27') 
});


test_that("Fc Paper figure 3a fetch succeeds (random shifted)", {
        ts <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0+offset, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet="GHCNv2",
                reproduceFor='2015-05-27')
});

test_that("Fc Paper figure 3b fetch succeeds (random shifted)", {
        ts2 <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0+offset, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet ="NCEP/NCAR Reanalysis 1 (regular grid)",
                reproduceFor='2015-05-27')    
});
test_that("Fc Paper figure 3c fetch succeeds (random shifted)", {
        ts3 <- fcTimeSeriesYearly(
                variable="airt",
                latitude=8.0+offset, longitude=10.0,
                firstDay=152,lastDay=243,
                firstYear=1950,lastYear=2050,
                url='http://fetchclim.cloudapp.net/',
                dataSet ="CESM1-BGC airt",
                reproduceFor='2015-05-27')
});

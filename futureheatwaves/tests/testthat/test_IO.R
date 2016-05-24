library(futureheatwaves)
library(dplyr)
context("Input / output of climate projection files and model info")

test_cities <- data.frame(city = c("balt", "ny", "nwk"),
                          lat = c(39.30080, 40.66980, 40.72410),
                          long = c(283.3894, 286.0562, 285.8268))

test_cities_2 <- test_cities %>%
        rename(latitude = lat,
               longitude = long) %>%
        mutate(extra = rep(1)) %>%
        select(extra, longitude, latitude, city)

expect_equal(process_cities_file(test_cities, c("lat", "long")),
             process_cities_file(test_cities_2, c("latitude", "longitude")))

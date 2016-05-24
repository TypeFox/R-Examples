pop.convert<- function(){
# data(pop)
### script written by Beth Ebert to convert data into binary obs.
### to convert the pop text into binary observations ###
### Note: cat0 = rain <= 0.2 mm 
###       cat1 = 0.2 < rain <= 4.4 mm
###       cat2 = 4.4 mm < rain
### Make observations into logical variables
d <- verification::pop
d$obs_norain <- d$obs <= 0.2 
d$obs_light  <- d$obs > 0.2 & d$obs <= 4.4    # light rain only 
d$obs_heavy  <- d$obs > 4.4   # heavy rain only
d$obs_rain   <- d$obs_light | d$obs_heavy


### Rename probability variables and compute probabilities for all
### rain

d$p24_norain <- d$p24_cat0 
d$p24_light  <- d$p24_cat1 
d$p24_heavy <- d$p24_cat2

d$p48_norain <- d$p48_cat0 
d$p48_light  <- d$p48_cat1 
d$p48_heavy  <- d$p48_cat2

d$p24_rain <- d$p24_light + d$p24_heavy 
d$p48_rain <- d$p48_light + d$p48_heavy 
# assign("d",d, envir = .GlobalEnv)
return(invisible(d) )
}

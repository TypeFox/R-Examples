hmdcountry <- function(Country, sex, username, password)
{
 file = c("Births.txt","Deaths_1x1.txt","Population.txt","Exposures_1x1.txt","Mx_1x1.txt","E0per.txt")
 yname = c("Birth count", "Death count", "Population", "Exposure", "Mortality rate", "Life expectancy")
 biglist <- list()
 for(i in 1:6){
     biglist[[i]] <- read.hmd(country = Country, sex = sex, file = file[i], username = username, password = password, yname = yname[i])
 }
 data = list(Birthcount = biglist[[1]], Deathcount = biglist[[2]], Population = biglist[[3]], Exposure = biglist[[4]], Mortalityrate = biglist[[5]], Lifeexpectancy = biglist[[6]])
 return(data)
}



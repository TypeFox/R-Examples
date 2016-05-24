hmdstatistic <-
function(sex, type = c("birth count","death count","population","exposure","mortality rate","life expectancy"), username, password)
{
   type = match.arg(type)
   Country=c("AUS","AUT","BLR","BEL","BGR","CAN","CHL","CZE","DNK",
             "DEUTNP","DEUTFRG", "DEUTGDR", "EST","FIN","FRATNP",
             "FRACNP","HUN","ISL","IRL","ISR","ITA","JPN","LVA",
             "LTU","LUX","NLD","NZL_NP","NZL_MA","NZL_NM","NOR",
             "POL","PRT","RUS","SVK","SVN","ESP","SWE","CHE","TWN",
             "GBR_NP","GBRTENW","GBRCENW","GBR_SCO","GBR_NIR","USA",
             "UKR")
   file = switch(type, "birth count" = "Births.txt", "death count" = "Deaths_1x1.txt", "population" = "Population.txt", "exposure" = "Exposures_1x1.txt", "mortality rate" = "Mx_1x1.txt", "life expectancy" = "E0per.txt")  
   yname = switch(type, "birth count" = "Birth count", "death count" = "Death count", "population" = "Population", "exposure" = "Exposure", "mortality rate" = "Mortality rate", "life expectancy" = "Life expectancy")
   biglist = list()
   for(i in 1:46){
       biglist[[i]] = read.hmd(country = Country[i], sex = sex, file = file, username = username, password = password, yname = yname)
   }
   data = list(Australia = biglist[[1]], Austria = biglist[[2]], Belarus = biglist[[3]], Belgium = biglist[[4]], 
               Bulgaria = biglist[[5]], Canada = biglist[[6]], Chile = biglist[[7]], CzechRepublic = biglist[[8]], 
               Denmark = biglist[[9]], Germanytotal = biglist[[10]], Germanywest = biglist[[11]], 
               Germanyeast = biglist[[12]], Estonia = biglist[[13]], Finland = biglist[[14]],
               Francetotal = biglist[[15]], Francecivilian = biglist[[16]], Hungary = biglist[[17]], Iceland = biglist[[18]],
               Ireland = biglist[[19]], Israel = biglist[[20]], Italy = biglist[[21]], Japan = biglist[[22]], 
               Latvia = biglist[[23]], Lithuania = biglist[[24]], Luxembourg = biglist[[25]], Netherlands = biglist[[26]], 
               NewZealandtotal = biglist[[27]], NewZealandmaori = biglist[[28]], NewZealandnonmaori = biglist[[29]], 
               Norway = biglist[[30]], Poland = biglist[[31]], Portugal = biglist[[32]], Russia = biglist[[33]], 
               Slovakia = biglist[[34]], Slovenia = biglist[[35]], Spain = biglist[[36]],
               Sweden = biglist[[37]], Switzerland = biglist[[38]], Taiwan = biglist[[39]], UKtotal = biglist[[40]], 
               EnglandWalestotal = biglist[[41]], EnglandWalescivilian = biglist[[42]], Scotland = biglist[[43]],
               NorthIreland = biglist[[44]], USA = biglist[[45]], Ukraine = biglist[[46]]) 
   return(data)
}



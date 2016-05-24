# TODO: Add comment
# 
# Author: ecor
###############################################################################


##############################################################
library(RSMET)
library(ggmap)

data(meteofrance)


dates <- as.Date(meteofrance$timestamp)

data=meteofrance[dates==dates[1],]


####
map <- get_map(location ="France", zoom = 6)

size <- 3

gsnow <- ggmap(map) +
		geom_point(data = data,aes(x = longitude, y = latitude),size=size,  alpha
						=1, color="blue",show.legend  = F)

## Uncomment if you want to save in PDF format the otput of gsnow
## ggsave("test-map.pdf", gsnow,width=10,height=10)

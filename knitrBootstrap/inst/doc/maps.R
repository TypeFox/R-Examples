## ----maps_arrests, dev='png', fig.show='hold', warning=FALSE-------------
library(ggplot2)
crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
library(reshape2) # for melt
crimesm <- melt(crimes, id = 1)
library(maps)
states_map <- map_data("state")
ggplot(crimes, aes(map_id = state)) + geom_map(aes(fill = Murder), map = states_map) + expand_limits(x = states_map$long, y = states_map$lat)
last_plot() + coord_map()
ggplot(crimesm, aes(map_id = state)) + geom_map(aes(fill = value), map = states_map) + expand_limits(x = states_map$long, y = states_map$lat) + facet_wrap( ~ variable)
great_lakes_states = c('michigan', 'illinois', 'ohio', 'wisconsin', 'indiana')
great_lakes_map = subset(states_map, region %in% great_lakes_states)
ggplot(subset(crimesm, state %in% great_lakes_states), aes(map_id = state)) + geom_map(aes(fill = value), map = great_lakes_map) + expand_limits(x=great_lakes_map$long, y=great_lakes_map$lat) + facet_wrap( ~ variable)


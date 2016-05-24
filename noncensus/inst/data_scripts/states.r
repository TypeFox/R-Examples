library(dplyr)
library(datasets)

# For info about census regions, see:
# http://en.wikipedia.org/wiki/List_of_regions_of_the_United_States#Census_Bureau-designated_regions_and_divisions

# Region - Midwest

# Division - East North Central
east_north_central <- data.frame(
    state=c("IL", "IN", "MI", "OH", "WI"),
    region="Midwest",
    division="East North Central"
)

# Division - West North Central
west_north_central <- data.frame(
    state=c("IA", "KS", "MN", "MO", "ND", "NE", "SD"),
    region="Midwest",
    division="West North Central"
)

# Region - Northeast

# Division - Mid-Atlantic
mid_atlantic <- data.frame(
    state=c("NJ", "NY", "PA"),
    region="Northeast",
    division="Mid-Atlantic"
)

# Division - New England
new_england <- data.frame(
    state=c("CT", "MA", "ME", "NH", "RI", "VT"),
    region="Northeast",
    division="New England"
)

# Region - South

# Division - East South Central
east_south_central <- data.frame(
    state=c("AL", "KY", "MS", "TN"),
    region="South",
    division="East South Central"
)

# Division - South Atlantic
south_atlantic <- data.frame(
    state=c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA", "WV"),
    region="South",
    division="South Atlantic"
)

# Division - West South Central
west_south_central <- data.frame(
    state=c("AR", "LA", "OK", "TX"),
    region="South",
    division="West South Central"
)

# Region - West

# Division - Mountain
mountain <- data.frame(
    state=c("AZ", "CO", "ID", "MT", "NM", "NV", "UT", "WY"),
    region="West",
    division="Mountain"
)

# Division - Pacific
pacific <- data.frame(
    state=c("AK", "CA", "HI", "OR", "WA"),
    region="West",
    division="Pacific"
)

states <- rbind(east_north_central,
                west_north_central,
                mid_atlantic,
                new_england,
                east_south_central,
                south_atlantic,
                west_south_central,
                mountain,
                pacific)
states$state <- as.character(states$state)
states <- states[order(states$state), ]


state_capitals <- c('Montgomery', 'Juneau', 'Phoenix', 'Little Rock',
                    'Sacramento', 'Denver', 'Hartford', 'Dover', 'Tallahassee',
                    'Atlanta', 'Honolulu', 'Boise', 'Springfield',
                    'Indianapolis', 'Des Moines', 'Topeka', 'Frankfort',
                    'Baton Rouge', 'Augusta', 'Annapolis', 'Boston', 'Lansing',
                    'Saint Paul', 'Jackson', 'Jefferson City', 'Helena',
                    'Lincoln', 'Carson City', 'Concord', 'Trenton', 'Santa Fe',
                    'Albany', 'Raleigh', 'Bismarck', 'Columbus', 'Oklahoma City',
                    'Salem', 'Harrisburg', 'Providence', 'Columbia', 'Pierre',
                    'Nashville', 'Austin', 'Salt Lake City', 'Montpelier',
                    'Richmond', 'Olympia', 'Charleston', 'Madison', 'Cheyenne')

state_names <- data.frame(state=state.abb,
                          name=state.name,
                          area=state.area,
                          capital=state_capitals,
                          stringsAsFactors=FALSE)

states <- left_join(states, state_names)
states <- within(states, {
                 name[state == "DC"] <- "District of Columbia"
                 area[state == "DC"] <- 68.3
                 })


# Linked from: http://www.census.gov/popest/data/state/totals/2013/index.html
census_data <- read.csv("http://www.census.gov/popest/data/national/totals/2013/files/NST_EST2013_ALLDATA.csv",
                        header=T)

states$population <- as.integer(merge(states,
                                      census_data,
                                      by.x="name",
                                      by.y="NAME")$CENSUS2010POP)

# Reorder columns
cols_states <- c('state', 'name', 'region', 'division', 'capital', 'area',
                 'population')
states <- states[cols_states]

# Additional info from:
# http://en.wikipedia.org/wiki/List_of_capitals_in_the_United_States#Insular_area_capitals
states <- rbind(states,
                c("AS", "American Samoa", NA, NA, "Pago Pago", 76.8, 55519),
                c("GU", "Guam", NA, NA, "Hagåtña", 212, 159358),
                c("MP", "Northern Mariana Islands", NA, NA, "Saipan", 179.01, 53833),
                c("PR", "Puerto Rico", NA, NA, "San Juan", 3515, 3725789),
                c("VI", "U.S. Virgin Islands", NA, NA, "Charlotte Amalie", 133.73, 106405))

states$state <- factor(states$state)
states$region <- factor(states$region)
states$division <- factor(states$division)


save(states, file="../../data/states.RData")

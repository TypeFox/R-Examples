setwd("C:/Mtech_Sem3")
require(RGoogleAnalyticsPremium)

# Authorize the Google Analytics account
# This need not be executed in every session once the token object is created
# and saved
token <- Auth(client.id,client.secret)

# Save the token object for future sessions
save(token,file="./token_files")
load("token_files")

ValidateToken(token)

# Initialize the query parameters
query.list <- Init(end.date = "2015-10-01",
                 metrics = "ga:users,ga:bounces,ga:avgTimeOnSite",
                 start.date = "2014-01-01",
                 title = "Behaviour_visits_and_OS",
                 dimensions = "ga:userType,ga:operatingSystem,ga:landingpagepath",
                 segment = "dynamic::ga:operatingSystem==Android")



# Create the Query Builder object so that the query parameters are validated
ga.query <- QueryBuilder(query.list)


# Specify the google analytics premium account id, webproperty id and view id
# from which you want to query the data
accountid <- "530XXX5"
webpropertyid <- "UA-530XXX5-1"
profileid <- "18XXXXX8"

# Extracts the unsampled data and stores it in R object (data)

path <- GetFile(ga.query, token, accountid, webpropertyid, profileid)
data <- read.csv(path,header = FALSE,sep = ",",comment.char = "#")
View(data)


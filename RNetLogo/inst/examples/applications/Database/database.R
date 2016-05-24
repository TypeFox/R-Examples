# load required libraries
library(RNetLogo)

# Note: there are different packages providing database access
#       for example: RMySQL, RPostgreSQL, ROracle, RJDBC and RODBC
# In this example, we will use RSQLite, which provides a connection to SQLite databases.
# These are small and easy-to-use databases and enough for a short example.
library(RSQLite)

# load database driver
m <- dbDriver("SQLite")

#----
# TODO: adapt these paths
#----
nl.path = "C:/Program Files/NetLogo 5.3/app"
database.path = "C:/Users/jthiele/Documents/test_netlogo.db"

# create connection to the database 
# (if the database does not exist alredy, this will create a file test_netlogo.db)
con <- dbConnect(m, dbname = database.path)
   
# start NetLogo session
NLStart(nl.path, gui=FALSE)
# load a NetLogo Model
NLLoadModel(paste(c(nl.path,"models/Sample Models/Earth Science/Fire.nlogo"),collapse="/"))
# setup the model
NLCommand("setup")

# run the model for 10 time steps and save the results (ticks and burned-trees) in table "Burned1" of database
dbWriteTable(con, "Burned1", 
             NLDoReport(10,"go",c("ticks","burned-trees"),
                        as.data.frame=T,df.col.names=c("time","burned")), 
             row.names=F, append=F)

# our first query: how many lines has the new table? 
dbGetQuery(con, "select count(*) from Burned1")[[1]]

# our second query: select all row from table Burned10 where time is more than 5
rs <- dbSendQuery(con, "select * from Burned1 where time > 5")
# get the result of the query
data <- fetch(rs, n = -1)
# show the result
print(data)

# clear query
dbClearResult(rs)

# append further results to existing table
dbWriteTable(con, "Burned1", 
             NLDoReport(10,"go",c("ticks","burned-trees"),
                        as.data.frame=T,df.col.names=c("time","burned")), 
             row.names=F, append=T)

# show content of our table
dbGetQuery(con, "select * from Burned1")

#---------------------------------------------------


# create a second table and save the result of 10 repeated simulation of 20 simulation steps each
for (x in 1:10)
{
NLCommand("setup")
dbWriteTable(con, "Burned2", 
             NLDoReport(20,"go",c("ticks","burned-trees"),
                        as.data.frame=T,df.col.names=c("time","burned")), 
             row.names=F, append=T)
}
            
# calculate the mean of burned trees (out of the 10 repetitions) for each time
rs <- dbSendQuery(con, "select avg(burned) as mean_burned from Burned2 group by time")
# get the result of the query
data <- fetch(rs, n = -1)
# show the result
print(data)

# clear query
dbClearResult(rs)

# clean up
dbDisconnect(con)
NLQuit()
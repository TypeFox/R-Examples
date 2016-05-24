require(sqlutils)
require(RSQLite)

sqlfile <- paste(system.file(package='sqlutils'), '/db/students.db', sep='')
m <- dbDriver("SQLite")
conn <- dbConnect(m, dbname=sqlfile)

#This will return the path(s) where query files will be loaded from
sqlPaths()

#List of available queries
getQueries()

#Return documentation of the queries
sqldoc('StudentSummary')
sqldoc('StudentsInRange')

#Execute the query
q1 <- execQuery('StudentSummary', connection=conn)
head(q1)
#Can always get the SQL statement to examine
getSQL('StudentSummary')

q2 <- execQuery('StudentsInRange', connection=conn)
head(q2)
#This query that has parameters will have their values replaced.
getSQL('StudentsInRange')

#Cache query
fn <- tempfile(fileext='.rda')
q3 <- cacheQuery('StudentSummary', filename=fn, connection=conn)
names(q3); nrow(q3)

#Since this will read from the cache, we don't need to specify the connection.
q4 <- cacheQuery('StudentSummary', filename=fn) 
names(q4); nrow(q4)

#Clean-up our session
dbDisconnect(conn)

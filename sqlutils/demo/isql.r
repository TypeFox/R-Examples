require(sqlutils)
require(RSQLite)

sqlfile <- paste(system.file(package='sqlutils'), '/db/students.db', sep='')
m <- dbDriver("SQLite")
conn <- dbConnect(m, dbname=sqlfile)

hist <- isql(conn=conn, sql=getSQL('StudentSummary'))
names(hist)
hist[['commands']]
hist[['sql']]

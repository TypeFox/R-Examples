# Importing counter table from a database
# 
# Function \code{get_counter_table} imports counter table from a database.
# 
# @usage get_counter_table(dbname = 'sejmrp', user = 'reader', password = 'qux94874', 
#    host = 'services.mini.pw.edu.pl', sorted_by_id = TRUE)
# 
# @param dbname name of database; default: 'sejmrp'
# @param user name of user; default: 'reader'
# @param password password of database; default: 'qux94874'
# @param host name of host; default: 'services.mini.pw.edu.pl'
# @param sorted_by_id information if table should be sorted by id; default: TRUE
# 
# @return data frame
# 
# @examples
# \dontrun{
# counter <- get_counter_table()
# names(counter)
# # [1] 'id' 'what' 'date'}
#
# @note
# Default parameters use privilages of 'reader'. It can only SELECT data from database.
# 
# All information is stored in PostgreSQL database.
# 
# @author Piotr Smuda


get_counter_table <- function(dbname = "sejmrp", user = "reader", password = "qux94874", host = "services.mini.pw.edu.pl",
                              sorted_by_id = TRUE) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), is.logical(sorted_by_id))
    
    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    
    # reading table
    if (sorted_by_id) {
        counter <- dbGetQuery(database_diet, "SELECT * FROM counter ORDER BY id")
    } else {
        counter <- dbGetQuery(database_diet, "SELECT * FROM counter")
    }
    
    suppressWarnings(dbDisconnect(database_diet))
    return(counter)
} 

####################################################################################################################
undigest <- function(db, str.md5) {
  #using a given db, the function search for a md5 match nd return the correspondent plain string
  ans <- db$db.ref[grep(str.md5, db$db.ref[,2], fixed = TRUE), 1]
  ans
}
####################################################################################################################
db.check <- function(db) {
  #check if the variable parsed as a parameter isn't null
  if (is.null(db)) {
    #if null, it triggers an error (check AQSys.err.R for details)
    AQSys.err("4")
  } else {
    # check if db is a list with length = 3 (db.cas, db.sys and db.ref)
    if ((is.list(db)) & (length(db) == 3)) {
      #Every element in db must be a dataframe
      for (i in 1:length(db)) {
        #if not, it triggers an error (check AQSys.err.R for details)
        if (is.data.frame(db[[i]]) == FALSE)
          AQSys.err("3", k = names(db)[i])
      }
      # if db isn't a list with length = 3, it triggers an error (check AQSys.err.R for details)
    } else {
      AQSys.err("7")
    }
  }
}
####################################################################################################################

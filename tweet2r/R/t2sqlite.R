#' @export
#' @import streamR rgdal RSQLite

t2sqlite<-function(fileprefix, import=TRUE){

  #list all the files of the folder and get the total number of files
  files <- list.files(pattern = ".json")
  num_files<-length(files)
  
  #opena database connection
  con <- dbConnect(SQLite(), fileprefix)
  
  
  #--------parse tweets using a loop to append next files to database table--------------#
  for(i in 1:num_files) { 
    
    # ------1) parse the JSON files
    #get the files names
    filename=paste(fileprefix,i-1,".json", sep="")
    
    #parse tweets and create a data frame
    tweets <- parseTweets(filename, simplify = FALSE, verbose = TRUE)
    
    #remove bad character before import it to database
    tweets$text<-gsub("[^[:alnum:]///' ]","", tweets$text)
    
    #------2) populate table with tweets
    dbWriteTable(con, fileprefix ,tweets,overwrite=FALSE, append=TRUE)
    
  }

  
  
  
  #transform to sqlite timestamp format created_at column
  dbSendQuery(con, 
              paste("ALTER TABLE", fileprefix,"ADD t_trans TIMESTAMP WITH TIME ZONE;
                UPDATE", fileprefix, "SET t_trans=created_at::TIMESTAMP WITH TIME ZONE;
                ALTER TABLE", fileprefix, "drop created_at;
                ALTER TABLE ", fileprefix, " rename t_trans TO created_at;"))
  
  
  #create view with only with tweets with only geolocated tweets
  dbSendQuery(con,
              paste("CREATE TABLE geo",fileprefix," AS
                  SELECT * 
                  FROM ", fileprefix, " 
                  WHERE lat is not null", sep=""))
  
  if (import==TRUE){
    tweets=dbReadTable(con, fileprefix)
  }else{
    remove(tweets)
  }
  
  
  #close connection with sqlite
  dbDisconnect(con)
  
  message(paste("Database created in ", getwd(),  "and tweets imported as a data frame.", sep=" ")) 
    
  return(tweets)
}


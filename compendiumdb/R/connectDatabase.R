connectDatabase <-
function(user, password, host = "localhost", dbname = "compendium", port = 3306)
{
  Sys.setenv(CYGWIN="nodosfilewarning")
  connect <- dbConnect(MySQL(), user=user, password=password, host=host, port=port)

  rs <- dbSendQuery(connect,paste("use",dbname))
  return(list(connect=connect,user=user,password=password,host=host,port=port,dbname=dbname))
}


.rimpalaEnv <- new.env()
rimpala.init <- function(impala_home=NULL, libs ="/usr/lib/impala/lib") {
  
  if(file.exists(libs)==TRUE)
    
  {
    rimpala.CP <- c(list.files(libs, full.names=TRUE, pattern="jar$", recursive=FALSE)
                    , list.files(paste(system.file(package="RImpala"),"java", sep=.Platform$file.sep ),pattern="jar$", full.names=TRUE))
    
    assign("classpath", rimpala.CP, envir =.rimpalaEnv )
    .jinit(classpath=rimpala.CP)
    return ("Classpath added successfully")
  }
  else
  {
    print("Error: Classpath does not exist. Please check the path\n")
  }
}


rimpala.query <-function (Q="show tables",isDDL="false",fetchSize="10000") {
  .jinit(parameters="-Xmx10240m")
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  rs = impalaObj$query(Q,isDDL,fetchSize)
  
  if(is.jnull(rs))
  {
    stop("SQL error")
  }
  
  
  arr = rs$toArray();#made list to array
  li = .jevalArray(arr) #made the result a list
  rw = lapply(li, .jevalArray) #now they are rows
  
  result = Reduce ("rbind", rw, init=NULL) #bind the rows together
  colNames = result[1,]
  
  if(nrow(result)<3){
    if(nrow(result)==1 && result[[1]] == "true")
    {
      print("Query Executed successfully")
    }
    else{
      print("Query Returned <0 Rows>")
    }
  }else{
    onlyData = data.frame(result[3:nrow(result),],stringsAsFactors=FALSE,row.names=NULL)
    colnames(onlyData)=colNames
    colTypes = result[2,]
    colNum.int = grep("int", colTypes,ignore.case=TRUE)
    if(length(colNum.int>0))
    {
      for(i in colNum.int)
      {
        onlyData[[i]] = as.integer(onlyData[[i]])
      }
    }
    
    colNum.double = grep("double", colTypes, ignore.case=TRUE)
    colNum.float = grep("float", colTypes, ignore.case=TRUE)
    colNum.double=c(colNum.double, colNum.float)
    if(length(colNum.double>0))
    {
      for(i in colNum.double)
      {
        onlyData[[i]] = as.double(onlyData[[i]])
      }
    }
    row.names(onlyData)=NULL
    return (onlyData)
  }
}

rimpala.connect <- function(IP="localhost",port="21050",principal="noSasl"){
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  
  #building the connection string
  #concat auth= or principal= depending on the user input to argument principal
  if(principal=="noSasl")
  {
    principal = paste("auth=",principal,sep="");
  } else  {
    principal = paste("principal=",principal,sep="");
  }
  
  return(impalaObj$connect(IP,port,principal))
  
}

rimpala.close <- function()
{
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  return(impalaObj$close())
}

rimpala.usedatabase <-function(db)
{
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  return(impalaObj$usedatabase(db))
}

rimpala.invalidate <- function(table=" "){
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  return(impalaObj$invalidate(table))
}

rimpala.refresh <- function(table="table_name"){
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  return(impalaObj$refresh(table))
}


rimpala.showtables <- function()
{
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  rs = impalaObj$showtables()
  if(is.jnull(rs))
  {
    stop("SQL error")
  }
  arr = rs$toArray();#made list to array
  li = .jevalArray(arr) #made the result a list
  rw = lapply(li, .jevalArray) #now they are rows
  
  result = Reduce ("rbind", rw) #bind the rows together
  colNames = result[1,]
  if(nrow(result)==1)
  {
    print("No data in table")
  }
  else
  {
    onlyData = data.frame(result[2:nrow(result),],stringsAsFactors=FALSE)
    colnames(onlyData)=colNames
  }
  return(onlyData)
}

rimpala.showdatabases <- function()
{
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  rs = impalaObj$showdatabases()
  if(is.jnull(rs))
  {
    stop("SQL error")
  }
  arr = rs$toArray();#made list to array
  li = .jevalArray(arr) #made the result a list
  rw = lapply(li, .jevalArray) #now they are rows
  
  result = Reduce ("rbind", rw) #bind the rows together
  colNames = result[1,]
  if(nrow(result)==1)
  {
    print("No data in table")
  }
  else
  {
    onlyData = data.frame(result[2:nrow(result),],stringsAsFactors=FALSE)
    colnames(onlyData)=colNames
  }
  return(onlyData)
}

rimpala.describe <-function(table)
{
  impalaObj = .jnew("com.musigma.ird.bigdata.RImpala")
  rs = impalaObj$describe(table)
  if(is.jnull(rs))
  {
    stop("SQL error")
  }
  arr = rs$toArray();#made list to array
  li = .jevalArray(arr) #made the result a list
  rw = lapply(li, .jevalArray) #now they are rows
  
  result = Reduce ("rbind", rw) #bind the rows together
  colNames = result[1,]
  if(nrow(result)==1)
  {
    print("No data in table")
  }
  else
  {
    onlyData = data.frame(result[2:nrow(result),],stringsAsFactors=FALSE)
    colnames(onlyData)=colNames
  }
  return(onlyData)
}




#draw a graph: histogram
xhist <- function(v="")
{
  tryCatch(
{
  #create a data frame
  conf = fromJSON(paste(.libPaths()[1], "x.ent/www/config/ini.json", sep='/'))
  lst_f <- xfile(sep=":")
  value <- c()
  date <- c()
  for(i in 1:length(lst_f))
  {
    value <- c(value,0)
    date <- c(date,0)
  }
  d <- data.frame(file=lst_f,date=date,value_date=value,visible=value)
  data <- readLines(conf$result$file)
  #check entity
  reg = ":\\$:"
  for(i in 1:length(data))
  {
    #get name of file
    f <- unlist(strsplit(data[i],":"))[1]
    #find all on data
    if(v == "")
    {
      #update date
      #format dd.mm.yyyyy
      reg_date1 <- ":([[:digit:]]{2}).([[:digit:]]{2}).([[:digit:]]{4}):" 
      if(grepl(pattern = reg_date1, x = data[i]))
      {
        date <- str_extract(data[i],reg_date1)
        year <- sub(pattern = reg_date1, replacement = "\\3", x = date)
        month <- sub(pattern = reg_date1, replacement = "\\2", x = date)
        d[d$file == f,2] <- paste(month,year,sep=".")
        d[d$file == f,3] <-  as.numeric(paste(year,month,sep=""))
      }
      #format mm.yyyyy
      reg_date1 <- ":([[:digit:]]{2}).([[:digit:]]{4}):" 
      if(grepl(pattern = reg_date1, x = data[i]))
      {
        date <- str_extract(data[i],reg_date1)
        year <- sub(pattern = reg_date1, replacement = "\\2", x = date)
        month <- sub(pattern = reg_date1, replacement = "\\1", x = date)
        d[d$file == f,2] <- paste(month,year,sep=".")
        d[d$file == f,3] <-  as.numeric(paste(year,month,sep=""))
      }
      if(!is.na(f))
      {
        d[d$file == f,4] <- 1  
      }
    }
    else
    {
      #format dd.mm.yyyyy
      reg_date1 <- ":([[:digit:]]{2}).([[:digit:]]{2}).([[:digit:]]{4}):" 
      if(grepl(pattern = reg_date1, x = data[i]))
      {
        date <- str_extract(data[i],reg_date1)
        year <- sub(pattern = reg_date1, replacement = "\\3", x = date)
        month <- sub(pattern = reg_date1, replacement = "\\2", x = date)
        d[d$file == f,2] <- paste(month,year,sep=".")
        d[d$file == f,3] <-  as.numeric(paste(year,month,sep=""))
      }
      #format mm.yyyyy
      reg_date1 <- ":([[:digit:]]{2}).([[:digit:]]{4}):" 
      if(grepl(pattern = reg_date1, x = data[i]))
      {
        date <- str_extract(data[i],reg_date1)
        year <- sub(pattern = reg_date1, replacement = "\\2", x = date)
        month <- sub(pattern = reg_date1, replacement = "\\1", x = date)
        d[d$file == f,2] <- paste(month,year,sep=".")
        d[d$file == f,3] <-  as.numeric(paste(year,month,sep=""))
      } 
      
      #lst <- c(lst,str_count(x = tolower(data[i]), pattern = reg1 , sep="\n"))
      #check 0 or 1 at the end of relation
      if(grepl(pattern = ":1$", x = v) || grepl(pattern = ":0$", x = v))
      {
        reg_entity = paste(":",v,sep="")
      }
      else
      {
        reg_entity = paste(":",v,":",sep="")
      }
      
      count <- str_count(x = tolower(data[i]), pattern = tolower(reg_entity) , sep="\n")
      if(count > 0)
      {
        #update value in the 
        d[d$file == f,4] <-  count
      } 
    }
  }
  if(v=="")
  {
    ylabel = "All documents" 
  }
  else
  {
    ylabel = paste("Documents contain the key:",v,sep = " ")
  }
  if(length(d[!is.na(d$value_date) & (d$value_date > 0) & (d$visible > 0),3])>0)
  {
    par(las = 0)#load default
    par(mfrow=c(1,1))#1 row, 1 column
    h1 = hist(d[!is.na(d$value_date) & (d$value_date > 0) & (d$visible > 0),3],plot = FALSE)
    hist(d[!is.na(d$value_date) & (d$value_date > 0) & (d$visible > 0),3], ylim=c(0,max(h1$count)+2) ,breaks = 12, col="blue",labels=TRUE,xaxt="n",border = "pink",main="Histogram of bulletin: date",xlab="Date",ylab= ylabel)    
    axis(side=1, at = d$value_date, labels= d$date)  
  }
  else
  {
    print("No data available")
  }
  #return a data frame for users check
  return(d)
},
error=function(cond) {
  message("Error: Parameters are incorrect or there are problems in the paths, please check your parameters!")
},
warning=function(cond) {
  message("Warning: Parameters are incorrect or there are problems in the paths, please check your parameters!")
},
finally={
  rm(list=ls())
})
}
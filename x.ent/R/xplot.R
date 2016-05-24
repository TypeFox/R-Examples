#draw graph: plot
xplot <- function(v1="",v2="",t="")
{
  tryCatch(
{
  #create a data frame
  conf = fromJSON(paste(.libPaths()[1], "x.ent/www/config/ini.json", sep='/'))
  lst_f <- xfile(sep=":")
  init <- c()
  init_pos <- c()
  for(i in 1:length(lst_f))
  {
    init <- c(init,0)
    init_pos <- c(init_pos,1)
  }
  #create a table for stocking data
  d <- data.frame(file=lst_f,date=init,value_date=init,visible=init_pos)
  data <- readLines(conf$result$file)
  #check entity e1
  #check entity
  reg = ":\\$:"
  reg_req = "";
  if((v1 == "") && (v2 == ""))
  {
    stop("Entity v1 or v2 must have a vulue")
  }
  if(length(v1) > 1)
  {
    stop("Entity v1 has only 0 or 1 value")
  }
  if((v1 != "") && (v2 == ""))
  {
    #find only a field v1
    reg_req = paste(":",v1,":",sep="")
    d[,v1] <- init
  }
  if((v1 != "") && (v2 != ""))
  {
    for(j in 1:length(v2))
    {
      d[,paste(v1,"-",v2[j],sep="")] <- init
    }
  }
  if((v1 == "") && (v2 != ""))
  {
    for(j in 1:length(v2))
    {
      d[,v2[j]] <- init
    }
  }
  for(i in 1:length(data))
  {
    #get name of file
    f <- unlist(strsplit(data[i],":"))[1]
    #fill data to data frame
    #add value of year
    #format dd.mm.yyyy
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
    #add value of entity
    if((v1 != "") && (v2 == "")) 
    {
      #find only a field v1
      if(grepl(pattern=reg,x=data[i]))
      {
        count <- str_count(x = tolower(data[i]), pattern = tolower(reg_req) , sep="\n")
        if(count > 0)
        {
          #update value in the 
          d[d$file == f,5] <-   count
        } 
      }
    }
    else if((v1 == "") && (v2 != ""))
    {
      if(grepl(pattern=reg,x=data[i]))
      {
        for(j in 1:length(v2))
        {
          reg_req = paste(":",v2[j],":",sep="")
          count <- str_count(x = tolower(data[i]), pattern = tolower(reg_req) , sep="\n")
          if(count > 0)
          {
            #update value in the 
            d[d$file == f,j+4] <-   count
          }
        }
      }
    }
    else
    {
      #more than deux entities
      for(j in 1:length(v2))
      {
        reg_req = paste(":",v1,":",v2[j],":",sep="")
        count <- str_count(x = tolower(data[i]), pattern = tolower(reg_req) , sep="\n")
        if(count > 0)
        {
          #update value in the 
          d[d$file == f,j+4] <-  count #d[d$file == f,j+4] + count
        } 
      }
    }
  }
  #filte the data following the time
  #if time
  if(any(t != ""))
  {
    reg_date = "([[:digit:]]{2}).([[:digit:]]{4})"
    if(length(t) == 1) #a precise day 
    { 
      if(grepl(pattern = reg_date, x = t))
      {
        month1 <- sub(pattern = reg_date, replacement = "\\1", x = t[1])
        year1 <- sub(pattern = reg_date, replacement = "\\2", x = t[1])
        date1 = as.numeric(paste(year1,month1,sep=""))
        if(!is.na(date1))
        {
          d[d$value_date != date1,4] <-  0 
        }
      }
      else
      {
        stop("Format of year (mm.yyyy) isn't valid, please check again!")
      }
    }
    if(length(t) == 2)
    {
      if((grepl(pattern=reg_date, x=t[1],perl=FALSE)) && (grepl(pattern=reg_date, x=t[2],perl=FALSE)))
      {
        month1 <- sub(pattern = reg_date, replacement = "\\1", x = t[1])
        year1 <- sub(pattern = reg_date, replacement = "\\2", x = t[1])
        date1 = as.numeric(paste(year1,month1,sep=""))
        month2 <- sub(pattern = reg_date, replacement = "\\1", x = t[2])
        year2 <- sub(pattern = reg_date, replacement = "\\2", x = t[2])
        date2 = as.numeric(paste(year2,month2,sep=""))
        if(!is.na(date1) & !is.na(date2))
        {
          d[d$value_date < date1 | d$value_date > date2 ,4] <-  0
        }
      }
      else
      {
        stop("Format of year (mm.yyyy) isn't valid, please check again!")
      }
    }
    if(length(t) > 2)#interval of date
    {
      stop("value of date isn't valid, there are two choise: a date (ex 02.2010) or interval of date (02.2010,02.2011) , please check again!")
    }
  }
  #draw graphe
  par(las = 0)#load default
  par(mfcol = c(ncol(d)-3,1),mar = c(0.5, 4.0, 0.5, 0.5), oma=c(1, 1, 4, 2))
  test = d[d$value_date > 0 & d$visible == 1,]
  test <- test[order(test$value_date),]
  label_h = 0
  for(i in 5:ncol(d))
  {
    if(nrow(test) > 0)
    {
      plot(test[,i], axes = TRUE, col = ifelse(test[test$value_date > 0,i] > 0, "red","purple"), xaxt="n", ylim=c(0,1), xlab = "",cex=2.0, ylab = colnames(test)[i] ,pch=15 , lty=5)
      if(label_h == 0)
      {
        axis(3, at= 1:length(test[test$value_date>0,3]),labels=test$date,col = "violet", las = 2,col.axis = "blue",cex.lab=0.7,cex=0.7,cex.axis=0.7)
        label_h <- 1
      }
    }
    else
    {
      print("No data available")
    }
  } 
  if(nrow(test) > 0)
  {
    axis(1, at= 1:length(test[test$value_date>0,1]),labels=test$file, col = "violet", las = 2,col.axis = "blue",cex.lab=0.7,cex=0.7,cex.axis=0.7)
    #title("Comparison of every entity in documents")
    #legend("bottomleft", inset=.05, title="Visible",
    #       c("0","1"), fill=terrain.colors(3), horiz=TRUE)
  }
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
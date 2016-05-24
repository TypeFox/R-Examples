#draw a graph: proposition
xprop <- function(v1,v2,type=1)
{
  tryCatch(
{
  #this is a dataframe that stores the results extracted
  if(is.null(v1) || is.null(v2))
  {
    print("Two arguments are required")
  }
  else
  {
    dta <- xdata()
    #list file in the corpus
    lst_f <- xfile(sep=":")
    df <- NULL
    #this is a variable that stores for xprop
    for(i in 1:length(v1))
    {
      #add columns
      df1 <- data.frame("val" <- 0)
      df1[,v2] <- 0
      df1[,"cat"] <- 0
      for(j in 1:length(v2))
      {
        str <- paste(v1[i],":",v2[j],":",sep="")
        for(k in 1:ncol(dta))
        {
          rows <- grep(str,dta[,k])
          #if column val is egal 0, update value this column, else add new row
          if(length(rows) > 0)
          {
            for(l in 1:length(rows))
            {
              df1[1,"val"] <- v1[i]
              df1[1,v2[j]] <- 1
              df1[1,"cat"] <- v2[j]
              if(is.null(df))
              {
                df= df1
              }
              else
              {
                df <- rbind(df,df1) 
              }
              #reset
              df1[1,"val"] <- 0
              df1[1,v2] <- 0
              df1[1,"cat"] <- 0
            } 
          }
        }
      }
    }
    df = df[rowSums(df[v2]) == 1, ]
    if(type == 1)
    {
      bp = ggplot(df,aes(x =val ,fill = cat)) + geom_bar(position = "fill")
      bp + coord_flip()
    }
    else
    {
      bp = ggplot(df,aes(x =  cat ,fill = val)) +     geom_bar(position = "fill")
      bp + coord_flip()
    }
  }
},
error=function(cond) {
  message("Parameters are incorrect or there are problems in the paths, please check your parameters!")
},
warning=function(cond) {
  message("Parameters are incorrect or there are problems in the paths, please check your parameters!")
},
finally={
  rm(list=ls())
})
}
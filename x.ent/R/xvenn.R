xvenn <- function(v,e=NULL)
{
  tryCatch(
{
  if(length(v) <= 1)
  {
    print("This vector must be greater than a element!")
  }
  else
  {
    dta <- xdata()
    df <- NULL # create a data frame NULL
    if(is.null(e))
    {
      #don't calcule relation with others entities
      #get the list file of corpus
      lst_f <- xfile(sep=":")
      df <- data.frame(val=lst_f)
      for(i in 1:length(v))
      {
        df[,v[i]] <- 0
      }
      for(i in 1:length(v))
      {
        for(j in 1: nrow(dta))
        {
          cols <- grep(v[i],dta[j,])
          if(length(cols) > 0)
          {
            #update the data frame
            df[df$val == dta[j,1],i+1] <- df[df$val == dta[j,1],i+1] + 1
          }
        }
      }
    }
    else
    {
      #get values of every entities 
      for(i in 1:length(e))
      {
        v1 = xdata_value(e[i])[["value"]]
        df1 <- data.frame(val=v1)
        if(length(v1) > 0)
        {
          if(is.null(df))
          {
            df = df1
          }
          else
          {
            df = rbind(df,df1)
          }
        }
      }
      if(!is.null(df) > 0)
      {
        for(i in 1:length(v))
        {
          df[,v[i]] <- 0
        } 
      }
      #finish creation a data frame
      for(i in 1:length(e))
      {
        #list value of entity
        v1 = xdata_value(e[i])[["value"]]
        if(length(v1) > 0)
        {
          for(j in 1:length(v1))
          {
            for(k in 1:length(v))
            {
              reg <- paste(v[k],":",v1[j],":",sep="")
              #find in the results and update the data frame
              for(h in 1: nrow(dta))
              {
                cols <- grep(reg,dta[h,])
                if(length(cols) > 0)
                {
                  #update the data frame
                  #df[df$val == as.character(v1[j]),k+1] <- df[df$val == as.character(v1[j]),k+1] + 1
                  df[df$val == as.character(v1[j]),k+1] <- 1
                }
              }
            }
          }
        }
      }
      #finish calcule
    }
    #Data Synthesis
    if(!is.null(df))
    {
      x <- as.matrix(df[,2:ncol(df)])
      nprobes <- nrow(x)
      ncontrasts <- ncol(x)
      names <- colnames(x)
      labels <- c()
      if(is.null(names)) names <- paste("Group",1:ncontrasts)
      noutcomes <- 2^ncontrasts
      outcomes <- matrix(0,noutcomes,ncontrasts)
      colnames(outcomes) <- names
      for (j in 1:ncontrasts)
        outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
      for(j in 1:ncontrasts)
        outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
      for(i in 1:noutcomes)
      {
        name <- ""
        for(j in 1:ncontrasts)
        {
          if(outcomes[i,j] == 1)
          {
            if(name == "")
              name <- names[j]
            else
              name <- paste(name,"&",names[j],sep="")
          }
        }
        if(name != "")
          labels <- c(labels,name)
      } 
      xlist <- list()
      for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
      counts <- as.vector(table(xlist))
      counts <- counts[-1]
      names(counts) <- labels
      #create venn diagram
      vd <- venneuler(counts)
      par(las = 0)#load default
      par(mfrow=c(1,1))#1 row, 1 column
      plot(vd,main="Compare values of entities or relations appearing simultaneously in bulletins")
      return(counts) 
    }
    else
    {
      print("No data available!")
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
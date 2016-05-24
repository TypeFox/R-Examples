xdata_value <- function(v, sort = "a")
{
  tryCatch(
  {
    if(is.null(v))
    {
      print("You must enter an entity represented in the results, please!")
    }
    else
    {
      conf = fromJSON(paste(.libPaths()[1], "x.ent/www/config/ini.json", sep='/'))
      path = conf$result$file;
      lst_f <- xfile(sep=":")
      #get all tags in the file config
      lst_tag <- xentity();
      lst_tag <-add_unique(lst_tag,xrelation())
      dta <- data.frame()
      reg_ent = ":\\$:"
      reg_rel = ":\\$\\$:"
      if(v %in% lst_tag)#only a column
      {
        lst <- c()
        #data frame, local parameter, use as hash
        dt <- data.frame(value=character(0),freq=integer(0),stringsAsFactors=FALSE)#store values of column and frequency found
        text <- readLines(path)
        for(i in 1:length(text))
        {
          #get name of file
          f_name <- unlist(strsplit(text[i],":"))[1]
          #creat a regular expression ex=> filename:p:$:
          reg = paste(f_name,":",v,reg_ent,sep="")
          if(grepl(pattern=reg,x=text[i]))#entity
          {
            #get value
            eles <- unlist(strsplit(text[i],":"))
            for(j in 4:length(eles))
            {
              #lst <- add_unique(lst,eles[j])
              if(length(dt$value[dt$value %in% eles[j]]) > 0)
              {
                dt[dt$value == eles[j],2] = dt[dt$value == eles[j],2] + 1;
              }
              else#add new value
              {
                dt <- rbind(dt,data.frame(value=eles[j],freq=1))
              }
            } 
            next
          }
          #creat a regular expression ex=> filename:p:$$:
          reg = paste(f_name,":",v,reg_rel,sep="")
          if(grepl(pattern=reg,x=text[i]))#relation
          {              
            #replace le result file:p:s:$$:e1:e2:1
            result = gsub(reg, "",text[i], perl=TRUE)
            #result = gsub(":$", "",result, perl=TRUE)
            #lst <- add_unique(lst,result)
            if(length(dt$value[dt$value %in% result]) > 0)
            {
              dt[dt$value == result,2] = dt[dt$value == result,2] + 1;
            }
            else#add new value
            {
              dt <- rbind(dt,data.frame(value=result,freq=1))
            }
          }
        }
        #convert to vector
        value <- as.vector(dt$value)
        freq <- as.vector(dt$freq)
        dt <- data.frame(value=value,freq=freq)
        if(sort == "a")
        {
          dt <- dt[order(dt$value,-dt$freq),]
        }
        else
        {
          dt <- dt[order(-dt$freq,dt$value),]
        }
        lst <- dt$value
        frq <- dt$freq
        dta <- data.frame(value=lst,freq=frq)#global parameter
        if(nrow(dta)>0)
        {
          dta$freq = formatC(dta$freq,digits=0,format="f") 
        }
      }
      else
      {
        print("The tag isn't valid!")
        return(NULL)
      }
      return(dta)
    }
  },
  error=function(cond) {
    message("Parameters are incorrect or there are problems in paths, please check your parameters!")
  },
  warning=function(cond) {
    message("Parameters are incorrect or there are problems in paths, please check your parameters!")
  },
  finally={
    rm(list=ls())
  })
}
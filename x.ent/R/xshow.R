xshow <- function(e=NULL,sort="a")
{
  tryCatch(
  {
    dta <- data.frame()
    if(is.null(e))#argument default, get all data
    {
      #all columns
      dta = xdata();
    }
    else#get a column from user
    {
      if(length(e) == 1)#if 1 argument or argument default
      {
        dta = xdata_value(e,sort)
      }
      else
      {
        dta = xdata(e)
      }
    }
    path = paste(.libPaths()[1],"x.ent/www/output.html",sep="/")
    html= print(xtable(dta),"html",file=path)
    html = gsub(";", "<br/>",html, perl=TRUE)
    html = paste("<meta http-equiv=Content-Type content=text/html; charset=utf-8>",html,sep="")
    write(html, file=path)
    browseURL(path)
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
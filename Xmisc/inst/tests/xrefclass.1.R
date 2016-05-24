
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 11 23:30:45 EDT 2014 -0400 (Week 32)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


test00.xRefClass <-
  function()
{
  
  obj <- xRefClass$new()
  logme(obj)
  ##logme(obj$trace())

  try(obk <- obj$copy())
  ## Error in (structure(function ()  : unused argument (quote("xrefclass")) ##

  obk2 <- obj$copy2()
  logme(obk2)
  
}

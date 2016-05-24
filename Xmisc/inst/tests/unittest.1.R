
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Tue Aug 12 09:56:07 EDT 2014 -0400 (Week 32)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


test00.UnitTest <-
  function()
{
  
  pkg <- 'Xmisc'
  test.obj <- UnitTest$new(pkg=pkg)
  test.obj$defineme()
  logme(test.obj)
  
}

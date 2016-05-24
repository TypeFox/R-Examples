# Hinda Haned, December 2008
# haned@biomserv.univ-lyon1.fr

# Accesor $

#simugeno
setMethod("$","simugeno",function(x,name) {
    return(slot(x,name))
})


#simumix
setMethod("$","simumix",function(x,name) {
    return(slot(x,name))
})



#tabfreq
setMethod("$","tabfreq",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","simugeno",function(x,name,value) {
   slot(x,name,check=TRUE) <- value
  return(x)
})
setMethod("$<-","simumix",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})
setMethod("$<-","tabfreq",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})



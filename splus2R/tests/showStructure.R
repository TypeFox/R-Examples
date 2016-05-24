if(is.R()) library(splus2R)

a <- c(m=1, n=2)
b <- diag(1:3)
cc <- cbind(a=1:5, b=2:6, c=letters[1:5])
d <- data.frame(cc)
attr(d, "dup.row.names") <- TRUE
e <- ts(1:10, frequency = 4, start = c(1959, 2))
f <- list(a,b=b)
setClass("track", representation(x="numeric", y="numeric"))
g <- new("track", x=1:5, y=1:5)

showStructure(a)
showStructure(b)
showStructure(cc)
showStructure(d)
showStructure(e)
showStructure(f)
showStructure(g)  # prints with @ rather than $
showStructure(list(a=a, b=b))
showStructure(list(cc=cc, d, list(a,e)))
showStructure(list(a=a, g=g))

if(!is.R()){
  g2 <- signalSeries(data.frame(x = 11:20, y = 21:30), 1:10,
                    units = c("mm", "cm"))
  showStructure(g2)
}

if(is.R()){
  g3 <- g
  attr(g3, "fish") <- "bonito"
  showStructure(g3)
}

if(is.R()){
  # More objects, from help(typeof)
  typeof(unclass)
  # "builtin"
  showStructure(unclass)
  # function[ length 1]  class: function 
  typeof(zapsmall)
  # "closure"
  showStructure(zapsmall)
  # function[ length 1]  class: function 
  typeof(get("while"))
  # "special"
  showStructure(get("while"))
  # function[ length 1]  class: function 

  # I am not trying the various objects that help(typeof)
  # says are unlikely to be seen at the user level.
}

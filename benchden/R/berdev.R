lben01 <- function()  {
 return(
  list(
   name = nberdev(1)  ,
   peaks = 0.5,
   support = matrix(c(0,1), nrow = 1),
   breaks = bberdev(1)
  )                                                               
 )
}
lben02 <- function()  {
 return(
  list(
   name = nberdev(2)  ,
   peaks = 0,
   support = matrix(c(0,Inf), nrow = 1),
   breaks = bberdev(2)
  )
 )
}
lben03 <- function()  {
 return(
  list(
   name = nberdev(3)  ,
   peaks = 1,
   support = matrix(c(0,Inf), nrow = 1),
   breaks = bberdev(3)
  )
 )
}
lben04 <- function()  {
 return(
  list(
   name = nberdev(4)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(4)
  )
 )
}
lben05 <- function()  {
 return(
  list(
   name = nberdev(5)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(5)
  )
 )
}
lben06 <- function()  {
 return(
  list(
   name = nberdev(6)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(6)
  )
 )
}
lben07 <- function()  {
 return(
  list(
   name = nberdev(7)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(7)
  )
 )
}
lben08 <- function()  {
 return(
  list(
   name = nberdev(8)  ,
   peaks = 0,
   support = matrix(c(0,1), nrow = 1),
   breaks = bberdev(8)
  )
 )
}
lben09 <- function()  {
 return(
  list(
   name = nberdev(9)  ,
   peaks = 1,
   support = matrix(c(1,Inf), nrow = 1),
   breaks = bberdev(9)
  )
 )
}
lben10 <- function()  {
 return(
  list(
   name = nberdev(10)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(10)
  )
 )
}
lben11 <- function()  {
 return(
  list(
   name = nberdev(11)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(11)
  )
 )
}
lben12 <- function()  {
 return(
  list(
   name = nberdev(12)  ,
   peaks = exp(-1),                 
   support = matrix(c(0,Inf), nrow = 1),
   breaks = bberdev(12)
  )
 )
}
lben13 <- function()  {
 return(
  list(
   name = nberdev(13)  ,
   peaks = 0,
   support = matrix(c(-5,5), nrow = 1),
   breaks = bberdev(13)
  )
 )
}
lben14 <- function()  {
 return(
  list(
   name = nberdev(14)  ,
   peaks = 0,
   support = matrix(c(-1/exp(1)^2,1/exp(1)^2), nrow = 1),
   breaks = bberdev(14)
  )
 )
}
lben15 <- function()  {
 return(
  list(
   name = nberdev(15)  ,
   peaks = 0,
   support = matrix(c(0,1), nrow = 1),
   breaks = bberdev(15)
  )
 )
}
lben16 <- function()  {
 return(
  list(
   name = nberdev(16)  ,
   peaks = 0,
   support = matrix(c(-1,1), nrow = 1),
   breaks = bberdev(16)
  )
 )
}
lben17 <- function()  {
 return(
  list(
   name = nberdev(17)  ,
   peaks = 0.5,
   support = matrix(c(0,1), nrow = 1),
   breaks = bberdev(17)
  )
 )
}
lben18 <- function()  {
 return(
  list(
   name = nberdev(18)  ,
   peaks =  0,
   support = matrix(c(0,Inf), nrow = 1),
   breaks = bberdev(18)
  )
 )
}
lben19 <- function()  {
 return(
  list(
   name = nberdev(19)  ,
   peaks = 0,
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(19)
  )
 )
}
lben20 <- function()  {
 return(
  list(
   name = nberdev(20)  ,
   peaks = 1/9,              
   support = matrix(c(0,Inf), nrow = 1),
   breaks = bberdev(20)
  )
 )
}
lben21 <- function()  {     
 return(
  list(
   name = nberdev(21)  ,
   peaks = c(-20,0),
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(21)
  )
 )
}
lben22 <- function()  {
 return(
  list(
   name = nberdev(22)  ,
   peaks = c(0.0005446896, 1.442525),              
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(22)
  )
 )
}
lben23 <- function()  {
 return(
  list(
   name = nberdev(23)  ,
   peaks = c(-0.9969638, -0.4978001, 0, 0.4978001, 0.9969638),        
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(23)
  )
 )
}
lben24 <- function()  {
 return(
  list(
   name = nberdev(24)  ,
   peaks = c(-1.476191, 0.8095009, 1.952370, 2.523804, 2.809521, 2.952380),
   support = matrix(c(-Inf,Inf), nrow = 1),
   breaks = bberdev(24)
  )
 )
}
lben25 <- function()  {
 return(
  list(
   name = nberdev(25)  ,
   peaks = c(-0.1, 0.1),
   support = matrix(bberdev(25), nrow =2, byrow = TRUE),
   breaks = bberdev(25)
  )
 )
}
lben26 <- function()  {
 return(
  list(
   name = nberdev(26)  ,
   peaks = c(-20.05, 0, 20.05),
   support = matrix(bberdev(26), nrow = 3, byrow = TRUE),
   breaks = bberdev(26)
  )
 )
}
lben27 <- function()  {
 return(
  list(
   name = nberdev(27)  ,
   peaks = c(-9,-7,-5,-3,-1,1,3,5,7,9),
   support = matrix(c(-10,10), nrow = 1),
   breaks = bberdev(27)
  )
 )
}
lben28 <- function()  {
 return(
  list(
   name = nberdev(28)  ,
   peaks = c(0,1),
   support = matrix(c(0,1), nrow = 1),
   breaks = bberdev(28)
  )
 )
}


`berdev` <- function(dnum=1) {
    if (is.nan(dnum) || ! dnum %in% 1:28)
        stop("dnum must be between 1 and 28")
    return (
            eval(
                parse(text = sprintf("lben%02d()", dnum)) # evaluate "ben[dnum]()"-string
            )
        )
}
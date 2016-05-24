# Copyright 2000-2001 (c) Nicholas Lewin-Koh 
# modifications 2001-2003 (c) Roger Bivand
# modifications 201506 Nick Bearman


# Returns a vector of HSV values
# modifications 2003 (c) Renaud Lancelot
color.ramp <- function (nclass, color = "red", nvec = NULL, type = "q"){
  eq.split <- function(ncl){
    mult <- rep((1 / ncl), (ncl - 1))
    mult * seq(1, (ncl - 1))
    }
  color.list <- list(cname = c("blue", "green", "yellow", "red"), hsvcol = c(0.7, 0.375, 0.19, 0))
  cind <- match(color, color.list$cname)
### change from "if(nvec)" to "if(!is.null(nvec))"
  if(!is.null(nvec)){
    if(type == "q"){
      pr <- eq.split(nclass)
### changes in min, quantile and max
      brks <- c(min(nvec, na.rm = TRUE),
                quantile(nvec, pr, names = FALSE, na.rm = TRUE),
                max(nvec, na.rm = TRUE))
      brks <- unique(brks)
      classvec <- cut(nvec, brks, labels = FALSE, include.lowest = TRUE)
      ramp <- hsv(rep(color.list$hsvcol[cind], nclass), c(pr, 1))
      return(list(ramp = ramp, col.class = classvec, breaks=brks))
      }
  else
    if(type == "e"){
      pr <- eq.split(nclass)
### changes in min, range and max
      brks <- c(min(nvec, na.rm = TRUE),
                pr * diff(range(nvec, na.rm = TRUE)),
                max(nvec, na.rm = TRUE))
      brks <- unique(brks)
      classvec <- cut(nvec, brks, labels = FALSE, include.lowest = TRUE)
      ramp <- hsv(rep(color.list$hsvcol[cind], nclass), c(pr, 1))
      return(list(ramp = ramp, col.class = classvec, breaks=brks))
      }
    }
  return(NULL)
}

#updated leglabs function to allow order of entries to be reversed

leglabs <- function(vec, under="under", over="over", between="-", reverse = FALSE) {
  x <- vec
  lx <- length(x)
  if(lx < 3) stop("vector too short")
  if(reverse) { #if user asks order to be reversed
    #reorder variable
#    x<-reverseorder(x)
    x<-rev(x)
    #set new under and over values
    under <- "over"
    over <- "under"
  }  
  res <- character(lx-1)
  res[1] <- paste(under, x[2])
  for (i in 2:(lx-2)) res[i] <- paste(x[i], between, x[i+1])
  res[lx-1] <- paste(over, x[lx-1])
  res
}

#reverseorder <- function(old.var){
  #function to reverseorder any set of values, designed for use with legends in maptools library
  #copy variable to be reordered
#  new.var <- old.var
  #loop going through each value in turn
#  for (i in 1:length(old.var)) {
    #save it in the appropiate place in the new variable
#    new.var[length(old.var)+1-i] <- old.var[i]
#  } 
  #return new variable
#  return(new.var)
#}

#x 

#x
#The set of classification methods is large (Dent p. 145), but there are a few to remember:

#        * Equal Intervals ("Constant Interval"): each class has same difference in value
#        * Quantile (N-tile): each class has same number of units
#        * Natural Breaks: visual examination; manual determination 

#    Then a lot of ones that you might need to use once in a while
  # Arithmetic progression: constant increase (decrease) in "width" of class
  # Geometric progression: constant multiplier used to derive width of class
  # Jenk's Iterative ("optimal") minimize within class standard deviations (variance) [ESRI calls this "natural breaks"]
#    (see Dent 147-149 on use of F-ratio and weighting)
  # Arbitrary breaks: given externally (laws, regulations, natural process)
  # Standard deviations: statistical distribution
  # Nested Means works by successive halving at the mean (2,4,8,16, ...)

#(Chrisman)




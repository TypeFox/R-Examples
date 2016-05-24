#
# qtl.comparison.plot.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2013
# first written Nov, 2012
# Contains: qtl.comparison.plot
#

# qtl.comparison.plot
#
# DESCRIPTION:
#  function to create new map and save it in cross object
# PARAMETERS:
#   - cross1 - object of class cross
#   - cross2 - object of class cross
#   - chr - which chromosome whoudl be plotted
#   - ... - sent to scanone
# OUTPUT:
#  An object of class cross
#
qtl.comparison.plot <- function(cross1, cross2, chr, ...){
  if(missing(cross1))   stop("Provide two objects of class cross!\n")
  if(missing(cross2))  stop("Provide two objects of class cross!\n")
  if(!(any(class(cross1)=="cross")))   stop("Provide two objects of class cross!\n")
  if(!(any(class(cross2)=="cross")))  stop("Provide two objects of class cross!\n")
  res1 <- scanone(cross1,...)
  res2 <- scanone(cross2,...)
  if((missing(chr))){
    chr <- 1
  }
  if(length(chr)>1){
    chr <- chr[1]
    warning("Function can plot only one chromosome at the time, plotting first of selected chromosomes!")
  }
  if(!(chr%in%1:min(nchr(cross1),nchr(cross2)))){
    stop("Chromosome asked is not present in the cross object\n")
  }
  plot(res1,res2,col=c("orange","black"),chr=chr,main="Figure 2 - QTL profiles comparison.")
}
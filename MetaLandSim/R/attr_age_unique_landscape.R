attr_age_unique_landscape <- function(rland, span_age, position)
{
  
  l1 <- rland$nodes.characteristics
  l2 <- span_age[[position]]
  l3 <- cbind(l1,l2$age)
  if(class(rland)=="metapopulation")names(l3)[10] <- "age"
  if(class(rland)=="landscape")names(l3)[9] <- "age"
  rland$nodes.characteristics <- l3
  return(rland)
}
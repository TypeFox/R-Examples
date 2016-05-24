f <- function # title
### description
(x, ##<< arg x
 y
### arg y
 ){
  ##value<< a list with elements
  list(x=x, ##<< original x value
       y=y, ##<< original y value
       sum=x+y) ##<< their sum
  ##end<<
}

.result <- 
 list(f = list(definition = "f <- function # title\n### description\n(x, ##<< arg x\n y\n### arg y\n ){\n  ##value<< a list with elements\n  list(x=x, ##<< original x value\n       y=y, ##<< original y value\n       sum=x+y) ##<< their sum\n  ##end<<\n}",  
     description = "description", `item{y}` = "arg y", `item{x}` = "arg x",  
     value = "a list with elements\n\\item{x}{original x value}\n\\item{y}{original y value}\n\\item{sum}{their sum}",  
     title = "title", format = "")) 

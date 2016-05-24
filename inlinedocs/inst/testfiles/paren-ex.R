myFun <- structure(

function # My function
  ### Fun. description
(
  x ##<< an argument
)
{
  x * 2
},

ex=function()
{
  x <- 5
  myFun(x)
})

.result <- 
 list(myFun = list(definition = "myFun <- structure(\n\nfunction # My function\n  ### Fun. description\n(\n  x ##<< an argument\n)\n{\n  x * 2\n},\n\nex=function()\n{\n  x <- 5\n  myFun(x)\n})",  
     description = "Fun. description", `item{x}` = "an argument",  
     title = "My function", format = "", examples = "\nx <- 5\nmyFun(x)\n")) 

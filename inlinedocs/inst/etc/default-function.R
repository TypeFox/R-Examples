## bt_jannis@yahoo.de



## I have a problem when I use functions as default values for argumnents
## in other functions. When I use curly brackets { here, I can not create a
## package with inlinedocs. It will give me the error when using
## package.skeleton() in my package structure:

## Error in parse(text = utxt) : <text>:4:0: unexpected end of input


## For example:

dummyfunction = function(filters = function(x) {b = 0; x > b} )
{
   # rest of code here
   return(filters)
}


## This seems to me as a legal function declaration but creates the above
## mentioned error. Is this an error of inlinedocs or do I misunderstand
## the R language? Or is there another way of using functions in such a way
## as arguments? In this case I could easily define this filters argument
## inside the function for cases when it is not supplied as an argument but
## I have some more complex functions where I really need to define
## something sequential as an argument like:

## dummyfunction = function(filters = {a = 1; b > a; b}) {print('test')}

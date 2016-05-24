genlibDLLC <- function(x)
{
	## TODO: Modify the codes below
	## Example of .C() with two arguments
	## See genlibDLL.cxx for the C/C++ (VC++) implementation of genlibDLLC() 

	len = length(x)
	y = .Call("genlibDLLC", as.double(x),	y=double(len),	len )$y
	y
}

genlibDLLCall <- function(x)
{
	## TODO: Modify the codes below
	## Example of .Call() with one argument
	## See genlibDLL.cxx for the C/C++ (VC++) implementation of genlibDLLCall() 

	.Call("genlibDLLCall", as.double(x))
}

genlibDLLCall2 <- function(x)
{
	## TODO: Modify the codes below
	## Example of .Call() with one argument
	## See genlibDLL.cxx for the C/C++ (VC++) implementation of genlibDLLCallSPL() 

	.Call("genlibDLLCall2", x)
}

genlibDLLTest <- function()
{
	v  = rnorm(10);
	vv = v*v
	vvC = genlibDLLC(v)
	vvCall = genlibDLLCall(v)
	vvCallSPL = genlibDLLCall2(v)
}

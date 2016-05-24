dircor.default<- function(veg,x.axis,y.axis,step,...)  {
	o.dircor<- dircor2(veg,x.axis,y.axis,step)
	o.dircor$call<- match.call()
	cat("Call:\n") 
	class(o.dircor) <- "dircor"
	print(o.dircor$call)
	o.dircor
}

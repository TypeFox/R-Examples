e1 = new.env()
assign("f", function(x)"f1", envir = e1)
e3 = new.env()
assign("g", function(x)f(x), envir = e3)
attach(e1)
attach(e3)
## Our intention is that g() calls f() in e1
g(1) #OK
## But now a new version of f() appears
e2 = new.env()
assign("f", function(x)"f2", envir = e2)
attach(e2)
search()
g(1)
## Because of search list position, the wrong function was called
## Furthermore, even having f() in the same environment as g() doesn't help:
assign("f", function(x)"f3", envir = e3)
objects(3)
detach(3)
attach(e3, pos=3)
g(1)
## And the same danger exists if the call to f() is from within the environment
assign("h", function(x)g(x), envir = e3)
detach(3)
attach(e3, pos=3)
h(1)

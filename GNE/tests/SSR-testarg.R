if(!require("GNE"))stop("this test requires package GNE.")


try(funSSR("x"))

try(funSSR(rep(1, 2) , 1, grobj=function(x) x))

try(funSSR(rep(1, 2), rep(1, 2), grobj=function(x) x))

try(funSSR(rep(1, 2), rep(1, 2), grobj=function(x) x, compl=phiFB))

funSSR(rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB)

try(funSSR(rep(1, 2), rep(1, 2), grobj=function(x, i, j, arg) arg$x, compl=phiFB))

funSSR(rep(1, 2), rep(1, 2), grobj=function(x, i, j, arg) arg$x, arggrobj=list(x=1), compl=phiFB)

try(funSSR(rep(1, 2), rep(1, 2), rep(1, 3), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x) x))

try(funSSR(rep(1, 2), rep(1, 2), rep(1, 3), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x) x, grconstr=function(x) 1))

try(funSSR(rep(1, 2), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x) x, grconstr=function(x) 1))

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x) x, grconstr=function(x) 1))

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x) 1))

funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j) 1)

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i, arg) arg$x, grconstr=function(x, i, j) 1))

funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i, arg) arg$x, argconstr=list(x=1), grconstr=function(x, i, j) 1)

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j, arg) arg$x, 
	arggrconstr=list(x=1)))

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j) 1, 
	joint=function(x) x))

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j) 1, 
	joint=function(x) x, grjoint=function(x) 1, dimmu=1:2))

try(funSSR(rep(1, 4), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j) 1, 
	joint=function(x) x, grjoint=function(x) 1, dimmu=1))

funSSR(rep(1, 5), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	constr=function(x, i) x, grconstr=function(x, i, j) 1, 
	joint=function(x) x, grjoint=function(x, j) 1, dimmu=1)

try(funSSR(rep(1, 5), rep(1, 2), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	joint=function(x) x, grjoint=function(x, j) 1, dimmu=1))

funSSR(rep(1, 3), rep(1, 2), grobj=function(x, i, j) x, compl=phiFB,
	joint=function(x) x, grjoint=function(x, j) 1, dimmu=1)
	
	

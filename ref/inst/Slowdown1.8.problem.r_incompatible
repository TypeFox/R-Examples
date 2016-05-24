

r-devel@stat.math.ethz.ch

dramatic slowdown under R 1.8 compared to previous versions


Dear R-developers,

I am about to submit a package "ref" which implements simple references for R. Since it slowed down Kurt Hornik's PC (sorry Kurt)
I discovered, that part of the functionality is by order of magnitudes slower under R 1.8 compared to previous versions (1.7, 1.6, 1.5).
My tests say the same slow down happens to Henrik Bengtsson's approach described at http://www.maths.lth.se/help/R/ImplementingReferences.

It would be great if someone with knowledge of the R internals could find out why this happens and whether this can be avoided in future versions.

Thanks for this great software
and have a nice weekend


Jens Oehlschlägel


# Slow down example
nmatrix <- 1000
nloop <- 10
m <- matrix(nrow=nmatrix, ncol=nmatrix)
rm <- as.ref(m)
# the following is by order of magnitudes (factor 20 to 80) slower under R 1.8 compared to R 1.7, 1.6, 1.5
system.time(
      for(i in 1:nloop)
        deref(rm)[1,1] <- i
)





# where deref and the other functions are defined as

"deref<-" <-
function(ref, value)
{
  assign(ref$name, value, envir=ref$loc)
  ref
}

deref <-
function (ref)
{
    get(ref$name, envir = ref$loc)
}

ref <-
function (name, loc = parent.frame())
{
    temp <- list(name = name, loc = loc)
    class(temp) <- "ref"
    temp
}

as.ref <-
function (obj)
{
    obj.name <- substitute(obj)
    obj.loc <- parent.frame()
    if (!is.name(obj.name))
        stop("obj must be a named object")
    obj.name <- deparse(obj.name)
    if (is.ref(obj)) {
        obj
    }
    else {
        ref(obj.name, obj.loc)
    }
}

is.ref <- function(x)
{
  inherits(x, "ref")
}



> version
         _
platform i386-pc-mingw32
arch     i386
os       mingw32
system   i386, mingw32
status
major    1
minor    8.1
year     2003
month    11
day      21
language R



Dear r-developers,

on January 23rd I reported a dramatic performance loss of R1.8ff compared to previous versions.
As I did not receive any answer for this important observation, I try again with some more recent results.

The slowdown happens when manipulating parts of objects in other environments.
This severely impedes usage of Henrik Bengtsons oo extensions and of my package ref.
It would be great if it would be possible to recover the great performance R had until version 1.7.
Identifying or solving this problem is far beyond my understanding of the interpreters internals, so please forgive me just sending this report without offering any patches.

You find my test code below, in case you like to have a look on it. I would appreciate your comments very much.

Best regards



Jens Oehlschlägel


# This example shows the slowdown (where as.ref and deref are defined as below)
nmatrix <- 1000
nloop <- 100
m <- matrix(nrow=nmatrix, ncol=nmatrix)
rm <- as.ref(m)
# the following is by order of magnitudes (factor 20 to 80) slower
# under R 1.8, R 1.9           (11.45 seconds)
# compared to R 1.7, 1.6, 1.5  ( 0.12 seconds)
system.time(
      for(i in 1:nloop)
        deref(rm)[1,1] <- i
)
# May be this helps finding the reason:
# You got the same slow behaviour under 1.7.0 if you turned
# deref() into method deref.ref() for a generic deref(ref, value)UseMethod("deref")
# Under 1.8 every object has a class. Does that mean every function has method dispatch under 1.8?
# If yes, can we bypass that by making deref() an internal function?


# What always works is using substitute directly,
# but this requires explicit programming of what the interpreter could do automatically.
# I don't think it is possible to write deref() such that it takes care about embedding subsetting.
# under R 1.8, R 1.9           (0.02 seconds)
# compared to R 1.7, 1.6, 1.5  (0.15 seconds)
system.time(
      for(i in 1:nloop)
        eval(substitute(m[1,1] <- i, list(i=i)), rm$loc)
)


# I also tried whether having deref use the the new 1.9.0 operators for environments helps,
# but it doesn't
# under R 1.9           (8.84 seconds)
system.time(
      for(i in 1:nloop)
        deref2(rm)[1,1] <- i
)

# Again, coding with this 1.9.0 feature explicitely gives the real speed
# under R 1.9           (0.03 seconds)
system.time(
      for(i in 1:nloop)
        rm$loc[["m"]][1,1] <- i
)






# Function definitions

"deref2<-" <-
function(ref, value)
{
  ref$loc[[ref$name]] <- value
  ref
}

deref2 <-
function (ref)
{
    ref$loc[[ref$name]]
}



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


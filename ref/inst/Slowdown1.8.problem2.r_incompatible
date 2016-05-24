

Dear Henrik,

congrats on your work on R objects, especially about "Implementing support for references in [R]".
A while ago I have submitted a related small package "ref" implementing references for R, which uses similar underlying techniques.

May be you have noted, that our appraoches give much slower performance under R1.8ff than before. I have asked r-devel once with no result.
The new version 1.9.0 has new operators for environments, which - by themselves - allow very fast read and write acces to objects in other environments.
However, using these operators instead of get() assign() in our context still suffers from slowdown.

Do you have any idea why that happens or who should be asked about that?
If not, I plan to submit the question to r-devel once again.

My test code is below, in case you like to have look on it.
It would be nice to hear about your experiences.

Best regards


Jens Oehlschlägel




# Slow down example
nmatrix <- 1000
nloop <- 100
m <- matrix(nrow=nmatrix, ncol=nmatrix)
rm <- as.ref(m)
# the following is by order of magnitudes (factor 20 to 80) slower under R 1.8 compared to R 1.7, 1.6, 1.5
system.time(
      for(i in 1:nloop)
        deref(rm)[1,1] <- i
)
# now using the new 1.9.0 operators for environments is also that slow
system.time(
      for(i in 1:nloop)
        deref2(rm)[1,1] <- i
)
# however, in principle the speed should be about that fast
system.time(
      for(i in 1:nloop)
        rm$loc[["m"]][1,1] <- i
)





# where deref and the other functions are defined as


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


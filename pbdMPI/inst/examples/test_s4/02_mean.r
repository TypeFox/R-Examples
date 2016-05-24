suppressMessages(library(pbdMPI, quietly = TRUE))

### Define method.
mean.simple <- function(x, ...){
  ret <- allreduce(c(sum(x@a), length(x@a)), op = "sum")
  ret[1] / ret[2]
}
.a <- setClass(Class = "simple",
               representation = representation(a = "numeric"))
.a <- setGeneric(name = "mean",
                 useAsDefault = mean)
.a <- setMethod(f = "mean",
                signature = signature(x = "simple"),
                definition = mean.simple)

### Run.
init()
x <- new("simple", a = comm.rank())
y <- mean(x)
comm.print(y, all.rank = TRUE)
finalize()

library(microbenchmark)
library(likeLTD)
library(ggplot2)

caseName = 'hammer'
datapath = file.path(system.file("extdata", package="likeLTD"), caseName)

args = list(
  databaseFile = NULL,
  cspFile    = file.path(datapath, 'hammer-CSP.csv'),
  refFile      = file.path(datapath, 'hammer-reference.csv'),
  nUnknowns    = 1,
  doDropin     = TRUE,
  ethnic       = "NDU1",
  adj          = 1.0,
  fst          = 0.02,
  relatedness  = c(0, 0)/4
)

# Create hypothesis for defence and prosecution.
prosecutionHyp = do.call(prosecution.hypothesis, args)
defenceHyp     = do.call(defence.hypothesis, args)

bench.any <- function(hypothesis, times=100L, ...) {

  objective <- create.likelihood.vectors(hypothesis, addAttr=TRUE, ...)
  funcs <- attr(objective, "functions") 
  arguments <- initial.arguments(hypothesis, ...)
  arguments$localAdjustment = 1
  microbenchmark( D3=do.call(funcs$D3, arguments), 
                  vWA=do.call(funcs$vWA, arguments), 
                  D16=do.call(funcs$D16, arguments), 
                  D2=do.call(funcs$D2, arguments), 
                  D8=do.call(funcs$D8, arguments), 
                  D21=do.call(funcs$D21, arguments), 
                  D18=do.call(funcs$D18, arguments), 
                  D19=do.call(funcs$D19, arguments), 
                  TH01=do.call(funcs$TH01, arguments), 
                  FGA=do.call(funcs$FGA, arguments), times=times )

}
bench.prosecution <- function(...) bench.any(prosecutionHyp, ...)
bench.defence <- function(...) bench.any(defenceHyp, ...)

result = list() 
result[["one"]] = bench.prosecution(nUnknowns=1, doDropin=TRUE)
result[["two"]] = bench.prosecution(nUnknowns=2, doDropin=TRUE)
result[["three"]] = bench.prosecution(nUnknowns=3, doDropin=TRUE, times=10)

nbthreads = .Call(.cpp.nbthreads, PACKAGE="likeLTD")


cat(sprintf("NUMBER OF THREADS: %d\n", nbthreads))

cat("\n\nPROSECUTION -- unknown=1, doDropin=TRUE\n")
print(result[["one"]])
ggplot(result[["one"]], aes(x=expr, y=time)) +
  geom_boxplot()   + 
  scale_y_log10()  +
  labs(title="Timing of likehood-per-locus functions with 1 unprofiled contributor",
       x="Locus", y="timings")

cat("\n\nPROSECUTION -- unknown=2, doDropin=TRUE\n")
print(result[["two"]])
ggplot(result[["two"]], aes(x=expr, y=time)) +
  geom_boxplot()   + 
  scale_y_log10()  +
  labs(title="Timing of likehood-per-locus functions with 2 unprofiled contributor",
       x="Locus", y="timings")

cat("\n\nPROSECUTION -- unknown=3, doDropin=TRUE\n")
print(result[["three"]])
ggplot(result[["three"]], aes(x=expr, y=time)) +
  geom_boxplot()   + 
  scale_y_log10()  +
  labs(title="Timing of likehood-per-locus functions with 3 unprofiled contributor",
       x="Locus", y="timings")

if(FALSE) {
  path = sprintf("likeLTD.timings.%d.C", nbthreads)
  save(result, file=path)
}

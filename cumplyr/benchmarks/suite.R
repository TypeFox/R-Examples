library('cumplyr')

commit <- system("git log --pretty=format:'%H' -n 1", intern = TRUE)
time <- Sys.time()

library(rbenchmark)

source(file.path('benchmarks', '01.R'))
source(file.path('benchmarks', '02.R'))

for (N in c(5, 10, 25, 50, 100))
{
  benchmarks <- benchmark(benchmark01(N),
                          benchmark02(N),
                          columns = c("test",
                                      "replications",
                                      "elapsed"),
                          replications = 10)
  benchmarks <- transform(benchmarks, Commit = commit)
  benchmarks <- transform(benchmarks, Time = time)
  benchmarks <- transform(benchmarks, ProblemSize = N)
  
  write.table(benchmarks,
              file = file.path('benchmarks', 'benchmarks.tsv'),
              append = TRUE,
              quote = FALSE,
              sep = '\t',
              col.names = FALSE,
              row.names = FALSE)
}

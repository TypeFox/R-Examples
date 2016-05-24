library('lubridate')

benchmarks <- read.csv(file.path('benchmarks', 'benchmarks.tsv'), sep = '\t')
benchmarks <- transform(benchmarks, Time = ymd_hms(Time))

# Generate separate plots for each benchmark.
ggplot(subset(benchmarks, Benchmark == 'benchmark02(N)'),
       aes(x = Time, y = Elapsed, group = ProblemSize, color = ProblemSize)) +
  geom_line() +
  scale_y_log10()
ggplot(subset(benchmarks, Benchmark == 'benchmark02(N)'),
       aes(x = ProblemSize, y = Elapsed, group = Commit, color = Commit)) +
  geom_line()

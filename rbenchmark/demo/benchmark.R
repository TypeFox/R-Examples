library(rbenchmark)

# Example 1
# Benchmarking the allocation of one 10^6-element numeric vector,
# by default replicated 100 times
benchmark(1:10^6)

# simple test functions used in subsequent examples
random.array = function(rows, cols, dist=rnorm) 
                  array(dist(rows*cols), c(rows, cols))
random.replicate = function(rows, cols, dist=rnorm)
                      replicate(cols, dist(rows))

# Example 2
# Benchmarking an expression multiple times with the same replication count,
# output with selected columns only
benchmark(replications=rep(100, 3),
          random.array(100, 100),
          random.array(100, 100),
          columns=c('test', 'elapsed', 'replications'))

# Example 3
# Benchmarking two named expressions with three different replication
# counts, output sorted by test name and replication count,
# with additional column added after the benchmark
within(benchmark(rep=random.replicate(100, 100),
                 arr=random.array(100, 100),
                 replications=10^(1:3),
                 columns=c('test', 'replications', 'elapsed'),
                 order=c('test', 'replications')),
       { average = elapsed/replications })

# Example 4
# Benchmarking a list of arbitrary predefined expressions
tests = list(rep=expression(random.replicate(100, 100)), 
             arr=expression(random.array(100, 100)))
do.call(benchmark,
        c(tests, list(replications=100,
                      columns=c('test', 'elapsed', 'replications'),
                      order='elapsed')))


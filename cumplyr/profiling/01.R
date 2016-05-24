library('cumplyr')

source(file.path('benchmarks', '01.R'))

Rprof(file.path('profiling', 'iddply.out'))
benchmark01(10000)
Rprof(NULL)

summaryRprof(file.path('profiling', 'iddply.out'))

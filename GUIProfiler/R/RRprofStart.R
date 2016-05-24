RRprofStart<-function(filename = "RRprof.out", interval = 0.02, numfiles = 100L, bufsize = 10000L){
  Rprof(line.profiling = TRUE,filename = filename,interval = interval, numfiles = numfiles, bufsize = bufsize)
  
}

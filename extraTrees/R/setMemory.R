setJavaMemory <- function(memoryInMB) {
  .jinit("ExtraTrees.jar", 
         parameters = sprintf("-Xmx%dm", memoryInMB), 
         force.init = TRUE)
}

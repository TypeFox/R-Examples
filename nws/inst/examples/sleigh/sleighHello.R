# This simple R sleigh example executes 10 tasks on the sleigh
# that return a friendly greeting

# initialize
library(nws)
s = sleigh()
ntasks = 10

# define the worker function
worker = function(x) {
  Sys.sleep(1)
  paste("hello", x, "from worker", SleighRank)
}

# execute the tasks and print the result list
print(eachElem(s, worker, 1:ntasks))

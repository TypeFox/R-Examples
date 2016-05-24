library(nws)

worker <- function() {
  id <- SleighRank
  nwsStore(SleighUserNws, sprintf("foo_%d", id), Sys.info()['nodename'])
  Sys.sleep(1)
  sprintf('%d %s', id, Sys.info()['nodename'])
}

# change launch if you add nodeList parameter
s <- sleigh()

# execution 1 (non-blocking eachWorker)
eo <- list(blocking=FALSE)
errmsg <- function(e) {
 cat("SUCCESS: generated a sleigh occupied error message\n\n")
}
sp <- eachWorker(s, worker, eo=eo)
tryCatch(eachWorker(s, worker, eo=eo), error=errmsg)
r <- waitSleigh(sp)
cat("Results from execution 1:\n")
print(r)

# execution 2 (blocking eachWorker)
r <- eachWorker(s, worker)
cat("Results from execution 2:\n")
print(r)

# execution 3
r <- eachWorker(s, worker)
cat("Results from execution 3:\n")
print(r)

# execution 4
r <- eachWorker(s, worker)
cat("Results from execution 4:\n")
print(r)

# execution 5
r <- eachWorker(s, worker)
cat("Results from execution 5:\n")
print(r)

cat("Loop over about 1 minute of work\n")
for (i in 1:60)
  r <- eachWorker(s, Sys.sleep, 1)

cat("Variable listing for sleigh workspace:\n")
write(nwsListVars(s@nws), stdout())

stopSleigh(s)

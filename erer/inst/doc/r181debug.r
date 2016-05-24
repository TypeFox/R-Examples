# A. Two bugs in a function: x has to be numeric; w is missing
miss <- function(x) {
  y <- x * 2
  lost <- function(y) {z <- y + w; return(z)}
  m <- lost(y)
  return(m + 100)
}
ua <- miss(x = "abc")
ub <- miss(x = 10)

# B. Special tools: need to change source codes --------------------------
# browser(): unconditional or conditional
miss2 <- function(x) {
  browser(text = "x is character!", expr = inherits(x, "character"))
  y <- x * 2
  if (inherits(x, "numeric")) browser()
  lost <- function(y) {z <- y + w; return(z)}
  m <- lost(y); return(m + 100)
}
uc <- miss2(x = "abc")  # Invoke the 1st conditional browser()
ud <- miss2(x = 10)     # Invoke the 2nd conditional browser()

# C. Special tools: no need to change source codes -----------------------
# C1. debug(), undebug(), debugonce(), isdebugged()
debug(fun = miss, text = "me"); isdebugged(miss)
ue <- miss(x = "abc")
undebug(miss); isdebugged(miss)

debugonce(miss) 
uf <- miss(x = 10)  # can debug(lost) in the middle of single-stepping
isdebugged(miss) 

# C2. trace(), untrace()
as.list(body(miss)); as.list(body(lm))
trace(what = miss, tracer = browser, at = c(4, 5))
ug <- miss(x = 10)
untrace(miss)

# C3. setBreakpoint(), unrace()
# The miss() function at the beginning is saved in a file.
source("C:/aErer/r181missFunction.r")  
setBreakpoint(srcfile = "r181missFunction.r", line = 3, clear = FALSE)
ug <- miss(x = 10)
untrace(miss)

# C4: Post-morten debugging: traceback()
uf <- miss(x = 10)
traceback()

# C5: Post-morten debugging: debugger()
options(error = dump.frames)
miss(x = 10)
debugger()

# C6: Post-morten debugging: recover
options(error = recover); miss(x = 10)
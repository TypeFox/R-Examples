## ---- include = FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = FALSE)

## ----check_gcc-----------------------------------------------------------
Sys.which("gcc")

## ----create_SFO_SFO------------------------------------------------------
library("mkin")
SFO_SFO <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))

## ----benchmark_SFO_SFO---------------------------------------------------
library("microbenchmark")
library("ggplot2")
mb.1 <- microbenchmark(
  "deSolve, not compiled" = mkinfit(SFO_SFO, FOCUS_2006_D, 
                                    solution_type = "deSolve", 
                                    use_compiled = FALSE, quiet = TRUE),
  "Eigenvalue based" = mkinfit(SFO_SFO, FOCUS_2006_D, 
                               solution_type = "eigen", quiet = TRUE),
  "deSolve, compiled" = mkinfit(SFO_SFO, FOCUS_2006_D, 
                                solution_type = "deSolve", quiet = TRUE),
  times = 3, control = list(warmup = 0))

smb.1 <- summary(mb.1)
print(mb.1)
autoplot(mb.1)

## ------------------------------------------------------------------------
rownames(smb.1) <- smb.1$expr
smb.1["median"]/smb.1["deSolve, compiled", "median"]

## ----benchmark_FOMC_SFO--------------------------------------------------
FOMC_SFO <- mkinmod(
  parent = mkinsub("FOMC", "m1"),
  m1 = mkinsub( "SFO"))

mb.2 <- microbenchmark(
  "deSolve, not compiled" = mkinfit(FOMC_SFO, FOCUS_2006_D, 
                                    use_compiled = FALSE, quiet = TRUE),
  "deSolve, compiled" = mkinfit(FOMC_SFO, FOCUS_2006_D, quiet = TRUE),
  times = 3, control = list(warmup = 0))
smb.2 <- summary(mb.2)
print(mb.2)
smb.2["median"]/smb.2["deSolve, compiled", "median"]
autoplot(mb.2)

## ----sessionInfo, echo = FALSE-------------------------------------------
cat(capture.output(sessionInfo())[1:3], sep = "\n")
if(!inherits(try(cpuinfo <- readLines("/proc/cpuinfo")), "try-error")) {
  cat(gsub("model name\t: ", "CPU model: ", cpuinfo[grep("model name", cpuinfo)[1]]))
}


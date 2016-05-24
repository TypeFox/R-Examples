###compute effective dimension when using tensor products (TP)

EDofTP <- function() {
  
#  cat(paste("\nComputing effective dimension . . .\n"))
  gc();
  t2 <- Sys.time()
#  if (taucs) {
#    BB.dgC <- as(BB, "dgCMatrix")  
#    ED <- (.C("hattraceopt", as.integer(nrow(S)), as.double(S@x),
#              as.integer(S@i), as.integer(S@p), as.double(BB.dgC@x),
#              as.integer(BB.dgC@i), as.integer(BB.dgC@p), as.double(0),
#              PACKAGE = "hattraceopt"))[[8]]
#  } else {
    Senv <- dsC2env(S)
    .C("myInversion", Senv$diag, Senv$env, Senv$xenv, Senv$N, Senv$Nenv,
       Senv$Nplus1, DUP = FALSE, PACKAGE = "svcm")  
    BBenv <- dsC2env(BB)
    ED <- as.double(numeric(1))
    .C("myTrace", Senv$diag, Senv$env, Senv$xenv, Senv$N, Senv$Nenv,
       Senv$Nplus1, BBenv$diag, BBenv$env, BBenv$xenv, BBenv$N, BBenv$Nenv,
       BBenv$Nplus1, ED = ED, DUP = FALSE, PACKAGE = "svcm")
#  }
  t3 <- Sys.time()
#  cat(paste(if (taucs) {"Taucs"} else {"Hutchinson/deHoog"}, "needed", t3-t2,
#            attributes(t3 - t2)$units, "and computed effective dimension =",
#            ED, ".\n"))
  assign("ED", ED, env = svcm.env)

}

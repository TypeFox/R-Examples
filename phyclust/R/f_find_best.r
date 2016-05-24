# This file contains functions to run all init method and init procedure.

#.init.procedure <- c("exhaustEM", "emEM", "RndEM", "RndpEM")
#.init.method <- c("randomMu", "NJ", "randomNJ", "PAM", "K-Medoids", "manualMu")
find.best <- function(X, K, EMC = .EMC, manual.id = NULL, byrow = TRUE,
    init.procedure = .init.procedure, init.method = .init.method,
    file.tmp = NULL, visible = FALSE, save.all = FALSE){
  if(K <= 0){
    stop("K > 0.")
  }
  if(byrow){
    X.bycol <- t(X)
  } else{
    X.bycol <- X
  }
  org.logL <- -Inf
  ret <- NULL
  new.ret <- NULL

  if(! exists(".Random.seed")){
    set.seed(12345)
  }

  if(save.all){
    ret.all <- NULL
    count.all <- 1
  }

  for(init.proc in init.procedure){
    for(init.meth in init.method){
      if(((init.proc != "exhaustEM")) &&
          (init.meth %in% c("NJ", "PAM", "manualMu")) ||
         ((K == 1) && (init.proc %in% c("NJ", "randomNJ"))) ||
         ((init.meth == "manualMu") && is.null(manual.id))
        ){
        next
      }

      EMC$init.procedure <- init.proc
      EMC$init.method <- init.meth

      seed.start <- .Random.seed
      if(!is.null(file.tmp)){
        save(list = c("EMC", "seed.start",
                      "new.ret", "ret", "org.logL"), file = file.tmp)
      }

      if(visible){
        cat("Run: ", init.proc, " and ", init.meth, "\n", sep = "")
      }

      new.ret <- try(phyclust(X.bycol, K, EMC, manual.id = manual.id,
                              byrow = FALSE), silent = TRUE)

      if(visible){
        print(new.ret)
        cat("\n")
      }

      seed.end <- .Random.seed
      if(!is.null(file.tmp)){
        save(list = c("EMC", "seed.start", "seed.end",
                      "new.ret", "ret", "org.logL"), file = file.tmp)
      }

      if(save.all){
        if(class(new.ret) == "try-error"){
          new.ret$init.procedure <- init.proc
          new.ret$init.method <- init.meth
        }

        new.ret$seed.start <- seed.start
        new.ret$seed.end <- seed.end

        ret.all[[count.all]] <- new.ret
        count.all <- count.all + 1
      }

      if(class(new.ret) != "try-error"){
        if(is.finite(new.ret$logL) && new.ret$logL > org.logL){
          org.logL <- new.ret$logL
          ret <- new.ret
        }
      }
    }
  }

  if(save.all){
    ret$save.all <- save.all
  }
  ret
} # End of find.best().


########################################################################
# PV-t-test MAIN FUNCTIONS
########################################################################
svyPVttest <- function(formula, 
                       design, 
                       placeholder = 1:10){
  # ERROR MESSAGE "NO PVs"
  if(length(grep("\\.\\.",formula)) == 0){
    stop("NO PVS -> use svyglm in package survey")}
  if(as.character(formula)[3] %in% c("0","1")){
    # prepare model t-test arguments
    as.char.forms <- lapply(placeholder,
                            function(x) 
                              (gsub("\\.\\.", x, formula)))
    char.calls <-  lapply(as.char.forms, 
                          function(x)
                            (paste(x[2],x[1], x[3])))
    formulas <- lapply(char.calls, as.formula)
    # compute test results
    mods <- lapply(lapply(formulas,as.formula), 
                   t.test.aux, 
                   design)
    # get test results & compute end results
    LIST <- lapply(mods, function(x)(x$DESC))
    TAB <- do.call(rbind,LIST)
    DESC.est <- TAB[,1]
    DESC.sds <- TAB[,2]
    DESC.var <- DESC.sds^2
    DESC.M <- mean(DESC.est)
    DESC.sampling.v <- mean(DESC.var)
    DESC.imputation.v <- var(DESC.est) 
    DESC.tot.se <- sqrt(DESC.sampling.v + 
                          ((1+1/length(placeholder)) * 
                             DESC.imputation.v))
    DESC.TAB.END.RES <- data.frame(mean=DESC.M, 
                                   se=DESC.tot.se)
    t.value <- DESC.M/DESC.tot.se  
    t.value.d <- abs(t.value)
    Pr.t <- pt(t.value.d, 
               degf(design), 
               lower.tail = FALSE)*2
    TEST.RES <- data.frame(t.value, 
                           degf(design), 
                           format.pval(Pr.t))
    colnames(TEST.RES) <- c("t.value", "degf", "Pr.t")
    rownames(TEST.RES) <- c("")
    N <- mods[[1]]$N
    FINAL.RES <- list(DESC = DESC.TAB.END.RES, 
                      TEST = TEST.RES, 
                      N=N)
    return(FINAL.RES)
  }
  else{
    # prepare model t-test arguments
    as.char.forms <- lapply(placeholder,
                            function(x) 
                              (gsub("\\.\\.", x, formula)))
    char.calls <-  lapply(as.char.forms, 
                          function(x)
                            (paste(x[2],x[1], x[3])))
    formulas <- lapply(char.calls, as.formula)
    # compute test results
    mods <- lapply(lapply(formulas,as.formula), 
                   t.test.aux, 
                   design)
    DESC.TAB <- lapply(mods, function(x)(x$DESC))
    # get test results & compute end results
    DESC.TAB.4.END.RES <- DESC.TAB[[1]]
    DESC.est <- do.call(rbind, 
                        lapply(DESC.TAB, 
                               function(x)
                                 (t(as.data.frame(x)[,2]))
                        )
    )
    DESC.sds <- do.call(rbind, 
                        lapply(DESC.TAB, 
                               function(x)
                                 (t(as.data.frame(x)[,3]))
                        )
    )
    DESC.var <- DESC.sds^2
    DESC.M <- apply(DESC.est, 2, mean)
    DESC.sampling.v <- apply(DESC.var, 2, mean)
    DESC.imputation.v <- apply(DESC.est, 2, var) 
    DESC.tot.se <- sqrt(DESC.sampling.v + 
                          ((1+1/length(placeholder)) * DESC.imputation.v))
    DESC.TAB.4.END.RES[,2] <- DESC.M
    DESC.TAB.4.END.RES[,3] <- DESC.tot.se
    DESC.TAB.END.RES <- DESC.TAB.4.END.RES
    
    TEST.TAB <- lapply(mods, function(x)(x$TEST))
    TEST.est <- do.call(rbind, 
                        lapply(TEST.TAB, 
                               function(x)
                                 (t(as.data.frame(x)[1,]))
                        )
    )
    TEST.sds <- do.call(rbind, 
                        lapply(TEST.TAB, 
                               function(x)
                                 (t(as.data.frame(x)[2,]))
                        )
    )
    TEST.var <- TEST.sds^2
    TEST.M <- mean(TEST.est)
    TEST.sampling.v <- mean(TEST.var)
    TEST.imputation.v <- var(TEST.est) 
    TEST.tot.se <- sqrt(TEST.sampling.v + 
                          ((1+1/length(placeholder)) * TEST.imputation.v))
    t.value <- TEST.M/TEST.tot.se  
    t.value.d <- abs(t.value)
    Pr.t <- pt(t.value.d, (degf(design)-1), lower.tail = FALSE)*2
    TEST.RES <- data.frame(t.value, (degf(design)-1), format.pval(Pr.t))
    colnames(TEST.RES) <- c("t.value", "degf", "Pr.t")
    rownames(TEST.RES) <- c("")
    N <- mods[[1]]$N
    FINAL.RES <- list(DESC = DESC.TAB.END.RES, 
                      TEST = TEST.RES, 
                      N=N)
    return(FINAL.RES)
  }
}



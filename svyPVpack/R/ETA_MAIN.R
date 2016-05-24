########################################################################
# ETA MAIN FUNCTIONS
########################################################################
svyPVeta <- function(formula, 
                     design, 
                     placeholder = 1:10){
  # prepare args
  # which var is factor? define y as factor
  char.args <- as.character(formula)[2:3]
  form.args <- all.vars(formula)
  # checking no pv args if factor is included
  # if no factor is included --> error message 
  # "NO factor included"
  # Index for no pv vars
  I.char <- grep("\\.\\.", char.args)
  I.form <- grep("\\.\\.", form.args)
  # Error message, when both vars are pvs
  if(length(I.char) == 2){stop("No factor included")}
  # args with no pvs included (have to include factor
  # otherwise --> ERROR message)
  fac.char.arg <- char.args[-I.char]
  fac.form.arg <- form.args[-I.form ]
  tar.form.arg <- form.args[I.form]
  # is a factor defined (in dataframe or formula)?
  YND <- is.factor(design$variables[[fac.form.arg]])
  YND <- as.numeric(YND)
  FACD <- length(grep("factor", fac.char.arg)) 
  FACD <- as.numeric(FACD)
  FACTORDUMMY <- YND + FACD
  # Error message when no factor is included
  if(FACTORDUMMY == 0){stop("No factor included")}
  x.name <- tar.form.arg 
  y.name <- fac.form.arg
  # create all PV names if PVs are included
  if(length(grep("\\.\\.", x.name))==1){
    x.name  <- lapply(placeholder,
                      function(x) gsub("\\.\\.", x, x.name))
  }
  x.name.4.df <- unlist(x.name)
  y.name.4.df <- y.name
  Y <- as.data.frame(design$variables[y.name.4.df])
  X <- as.data.frame(design$variables[x.name.4.df])
  # create NA-dummy to count cases/weights & for survey update
  NA.dummy.matrix <- cbind(X,Y)
  NA.dummy <- apply(NA.dummy.matrix, 
                    1, 
                    function(x){
                      any(is.na(x))}
  )
  n <- sum(design$pweights[NA.dummy==FALSE])
  N <- sum(NA.dummy==FALSE)
  updated.design <- update(design, NA.dummy = NA.dummy)
  # computation with no PV as input
  if(dim(as.matrix(X))[2] == 1){
    result <- cor.aux.func(X[NA.dummy == FALSE,], 
                           Y[NA.dummy == FALSE,], 
                           subset(updated.design, 
                                  NA.dummy == FALSE)
    )  
    finalres <- rbind(t(as.data.frame(result)), N, n) 
    rownames(finalres) <- c("COR", "SE", "number.of.cases",
                            "sum.of.weights")
    t(finalres)
  }
  else{
    Y.d <- matrix(data = unlist(rep(Y, dim(as.matrix(X))[2])), 
                  ncol = dim(as.matrix(X))[2])  
    PVs.d <- lapply(1:dim(as.matrix(Y.d))[2], 
                    function(n){
                      inp1 <- X[NA.dummy == FALSE,n]
                      inp2 <- Y.d[NA.dummy == FALSE,n]
                      eta.func(inp1, 
                               inp2, 
                               subset(updated.design, 
                                      NA.dummy == FALSE))}
    )
    # dataframe with coefficients
    unlPV.with.lab <- as.data.frame(PVs.d)
    # point estimators + SD
    ndummy1 <- seq(1,2*max(dim(as.matrix(X))[2], 
                           dim(as.matrix(Y))[2]),2)
    ndummy2 <- 1+seq(1,2*max(dim(as.matrix(X))[2],
                             dim(as.matrix(Y))[2]),2)
    PV.means <- unlist(unlPV.with.lab[ndummy1])
    PV.ses.d <- unlist(unlPV.with.lab[ndummy2])
    PV.ses <- PV.ses.d^2
    Means <- mean(PV.means)
    sampling.v <- mean(PV.ses)
    imputation.v <- var(PV.means) 
    tot.v <- sqrt(sampling.v + ((1+1/length(placeholder)) * 
                                  imputation.v))
    # compilation of final result object
    result <- rbind(Means,tot.v,N,n)
    rownames(result) <- c("ETA", "SE", "number.of.cases",
                          "sum.of.weights")
    # return final result object
    t(result)
  }
}
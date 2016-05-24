## ----include=FALSE-------------------------------------------------------
library(localsolver)

## ----include=FALSE-------------------------------------------------------
model <- "function model() {
  x[i in 1..4] <- bool();

  // weight constraint  
  knapsackWeight <- sum[i in 1..nbItems](itemWeights[i] * x[i]);
  constraint knapsackWeight <= knapsackBound;
  
  // maximize value
  knapsackValue <- sum[i in 1..nbItems](itemValues[i] * x[i]);
  maximize knapsackValue;
}"


## ------------------------------------------------------------------------
lsp <- ls.problem(model)
lsp

## ------------------------------------------------------------------------
set.params(lsp, lsTimeLimit=60)

## ------------------------------------------------------------------------
lsp <- set.params(lsp, lsTimeLimit=60)

## ------------------------------------------------------------------------
lsp <- set.params(lsp, lsTimeLimit=60, lsIterationLimit=250)
lsp <- set.params(lsp, lsTimeLimit=300, lsSeed=7)
lsp

## ------------------------------------------------------------------------
lsp <- reset.lsp.params(lsp)
lsp

## ------------------------------------------------------------------------
lsp <- add.output.expr(lsp, "x", 4)
lsp <- add.output.expr(lsp, "knapsackWeight")
lsp <- add.output.expr(lsp, "knapsackValue")
lsp

## ------------------------------------------------------------------------
lsp <- clear.output.exprs(lsp)
lsp

## ----include = FALSE-----------------------------------------------------
lsp <- add.output.expr(lsp, "x", 4)
lsp <- add.output.expr(lsp, "knapsackWeight")
lsp <- add.output.expr(lsp, "knapsackValue")
lsp <- set.params(lsp, lsTimeLimit=60, lsIterationLimit=250)

## ------------------------------------------------------------------------
newElement <- list()
newElement[[1]] <- c(1,2,3.14)
newElement[[4]] <- c(1.6,2.77, 5, 34, 1)
newElement[[6]] <- 3
newElement

## ------------------------------------------------------------------------
data <- list(nbItems=4L, itemWeights=c(1L,2L,3L,4L), 
             
             itemValues=c(5,6,7,8), knapsackBound = 9L)

## ------------------------------------------------------------------------
lsp
data
ls.solve(lsp, data)

## ----include=FALSE-------------------------------------------------------
model <- "function model() {
  x[i in 1..4] <- float(0,100);

  // time constraint
  productionTime <- sum[i in 1..4](time[i] * x[i]);
  constraint productionTime <= 200;
  
  // raw material constraint
  rawMaterials <- sum[i in 1..4](materialsR[i] * x[i]);
  constraint rawMaterials <= 300;

  // pre produced material constraint
  preMaterials <- sum[i in 1..4](materialsP[i] * x[i]);
  constraint  preMaterials <= 500;

  //maximize revenue
  revenue <- sum[i in 1..4](price[i] * x[i]);
  maximize revenue;

}"


## ------------------------------------------------------------------------
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=60, lsIterationLimit=250)
lsp <- add.output.expr(lsp, "x", 4)
lsp <- add.output.expr(lsp, "revenue")
lsp

data <- list(time=c(5L,2L,6L,3L), materialsR=c(10L,8L,9L,7L), 
             materialsP=c(10L,20L,25L,22L), price=c(5L,4L,6L,3L) )
ls.solve(lsp, data)

## ----include = FALSE-----------------------------------------------------
model <-
  "function model(){
    x[0..nbProcesses-1][0..nbMachines-1] <- bool();
    //some model formulation
  }
  
  function param(){
    for [p in 0..nbProcesses-1][m in 0..nbMachines-1]
      if (m == initialMachine[p]) setValue(x[p][m], true);
      else setValue(x[p][m], false);
    ltTimeLimit = 60;
  }"  


## ------------------------------------------------------------------------
lsp <- ls.problem(model)
lsp <- set.params(lsp, lsTimeLimit=300)
lsp

## ------------------------------------------------------------------------
lsp <- ls.problem("function model(){ ... }")
lsp
lsp <- set.temp.dir(lsp, path=file.path(tempdir(), '..'))
lsp


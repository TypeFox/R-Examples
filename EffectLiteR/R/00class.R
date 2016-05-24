######################## class definitions #####################

## structure of class effectlite
# - user input (class input)
# - parameter names (class parnames)
# - generated lavaansyntax (class lavsyntax)
# - obtained results (class results)

setClass("input", representation(
  vnames="list", ## variable names
  vlevels="list", ## variable levels (for x, k, kstar and cell)
  control="character",
  ng="integer", ## number of treatment groups
  nz="integer", ## number of z
  nk="integer", ## number of unfolded categories of K  
  data="data.frame", 
  measurement="character",
  add="character",
  fixed.cell ="logical",
  fixed.z ="logical",
  missing="character",
  observed.freq="numeric", ## observed group frequencies (fixed.cell only)
  sampmeanz="array", ## manifest sample means for continuous covriates
  se="character", ## lavaan standard errors
  bootstrap="numeric", ## number of bootstrap draws
  mimic="character", ## lavaan's mimic option
  interactions="character", ## type of interaction (all, 2-way, no)
  complexsurvey="list",
  homoscedasticity="logical",
  outprop="list" ## output from propensity score model
)
)

setClass("parnames", representation(
  alphas="array", 
  betas="array", 
  gammas="array", 
  gammalabels="array",
  cellmeanz="character",
  meanz="character",
  pk="character",
  px="character",
  Ezk="character",
  Pkgx="character", ## P(K=k|X=x)
  Pxgk="character", ## P(X=x|K=k)
  Ezgx="character", ## E(Z|X=x)
  Ezgk="character", ## E(Z|K=k)
  Ezkgx="character", ## E(Z*K|X=x)
  groupw="character",
  relfreq="character",
  Egx="character",
  Egxgx="character", ## E(gx|X=x)
  Egxgk="character", ## E(gx|K=k)
  Egxgxk="character", ## E(gx|X=x,K=k)
  adjmeans="character"
)
)

setClass("lavsyntax", representation(
  model="character", 
  hypotheses="list"
)
)


setClass("results", representation(
  lavresults="lavaan",
  hypotheses="data.frame",
  Egx="data.frame",
  Egxgx="data.frame",
  Egxgk="data.frame",
  Egxgxk="data.frame",
  gx="list",
  adjmeans="data.frame",
  condeffects="data.frame"
)
)


setClass("effectlite", representation(
  call="call",
  input="input",
  parnames="parnames",
  lavaansyntax="lavsyntax",
  results="results"
)
)

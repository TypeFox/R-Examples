setClass("SW", representation(Opt.W="matrix",First.C="vector",Second.C="vector"))
setClass("PAR",representation(fid="character",LB="numeric",UB="numeric",grid="numeric",ds="vector"))
setClass("OPT",representation(Par="PAR",Opt="matrix",Eff="numeric"))



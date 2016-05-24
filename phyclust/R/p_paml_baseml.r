### This file contains functions for "baseml" in PAML.

paml.baseml <- function(X, seqname = NULL, opts = NULL, newick.trees = NULL){
  ### Check arguments.
  if(is.null(X) || nrow(X) < 2){
    stop("Input sequences are not correct.")
  }
  if(is.null(opts)){
    opts <- paml.baseml.control()
  }
  if(is.null(newick.trees) && opts$runmode %in% 0:1){
    stop("Newick trees are requried for runmode 0 or 1.")
  }

  ### Create temporary directory and files.
  temp.dir <- tempfile("paml_baseml.")
  if(file.exists(temp.dir)){
    unlink(temp.dir, recursive = TRUE)
  } else{
    dir.create(temp.dir)
  }
  temp.file.tree <- paste(temp.dir, "/baseml.trees", sep = "")
  temp.file.nuc <- paste(temp.dir, "/baseml.nuc", sep = "")
  temp.file.control <- paste(temp.dir, "/baseml.ctl", sep = "")
  baseml.file.names <- c("baseml.trees", "baseml.nuc", "baseml.ctl")

  temp.file.stdout <- paste(temp.dir, "/stdout", sep = "")

  ### Dump trees to the temporary file if available.
  if(! is.null(newick.trees)){
    opts$treefile <- "baseml.trees"

    write(paste(length(newick.trees), " ", nrow(X), sep = ""),
          file = temp.file.tree)
    write(unlist(newick.trees, use.names = FALSE),
          file = temp.file.tree, sep = "\n", append = TRUE)
  }

  ### Dump sequences to the temporary file.
  write.paml(X, temp.file.nuc, seqname = seqname, code.type = .code.type[1])

  ### Dump control to the temporary file.
  opts$seqfile <- "baseml.nuc"
  opts$outfile <- "mlb"
  for(i.opts in names(opts)){
    cat(i.opts, " = ", opts[[i.opts]], "\n", file = temp.file.control,
        sep = "", append = TRUE)
  }

  ### Change path and run.
  current.dir <- getwd()
  setwd(temp.dir)
  argv <- c("baseml")
  tmp <- try(.Call("R_paml_baseml_main", argv, temp.file.stdout,
                   PACKAGE = "phyclust"),
             silent = TRUE)
  setwd(current.dir)

  ### Read in output.
  ret <- NULL
  if(class(tmp) == "try-error"){
    ret$error.baseml <- tmp
  }
  output.file.names <- list.files(temp.dir)
  for(i.file in output.file.names){
    if(i.file %in% baseml.file.names){
      next
    }

    ret[[i.file]] <- scan(file = paste(temp.dir, i.file, sep = "/"),
                     what = "character", sep = "\n", quiet = TRUE)
  }

  ### Find the best tree.
  tmp <- try(read.tree(text = ret$mlb[length(ret$mlb)]), silent = TRUE)
  if(class(tmp) == "try-error"){
    ret$error.tree <- tmp
  } else{
    id <- length(ret$mlb)
    ret$logL <- gsub(".*lnL: *", "", ret$mlb[id - 1])
    ret$best.tree <- ret$mlb[id]
  }
  class(ret) <- "baseml"

  ### Remove temporary directory and return.
  unlink(temp.dir, recursive = TRUE)
  return(ret)
} # End of paml.baseml().

print.baseml <- function(x, ...){
  if(!is.null(x$error.baseml)){
    print(x$error.baseml)
  }
  if(!is.null(x$error.tree)){
    print(x$error.tree)
  } else{
    cat("logL: ", x$logL, "\n", sep = "")
    cat(x$best.tree, "\n", sep = "")
  }
} # End of baseml().

paml.baseml.control <- function(...){
  default <- list(
    noisy = 0,         # 0,1,2,3: how much rubbish on the screen
    verbose = 0,       # 1: detailed output, 0: concise output
    runmode = 2,       # 0: user tree;  1: semi-automatic; 2: automatic
                       # 3: StepwiseAddition; (4,5):PerturbationNNI
    model = 0,         # 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                       # 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
    Mgene = 0,         # 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
    clock = 0,         # 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0,     # 0: estimate kappa; 1: fix kappa at value below
    kappa = 5,         # initial or fixed kappa
    fix_alpha = 0,     # 0: estimate alpha; 1: fix alpha at value below
    alpha = 0.5,       # initial or fixed alpha, 0:infinity (constant rate)
    Malpha = 0,        # 1: different alpha's for genes, 0: one alpha
    ncatG = 5,         # no. of categories in the dG, AdG, or nparK models of rates
    nparK = 0,         # rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK
    nhomo = 0,         # 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
    getSE = 0,         # 0: don't want them, 1: want S.E.s of estimates
    RateAncestor = 0,  # (0,1,2): rates (alpha>0) or ancestral states
    Small_Diff = 7e-6,
    cleandata = 1,     # remove sites with ambiguity data (1:yes, 0:no)?
    method = 0         # Optimization method 0: simultaneous; 1: one branch a time
  )
  default.names <- names(default)

  ret <- default
  control <- list(...)
  for(i.control in names(control)){
    if(i.control %in% c("seqfile", "treefile", "outfile")){
      next
    }
    ret[[i.control]] <- control[[i.control]]
  }

  return(ret)
} # End of paml.baseml.control().

paml.baseml.show.default <- function(){
  cat("Default options:
noisy = 0          # 0,1,2,3: how much rubbish on the screen
verbose = 0        # 1: detailed output, 0: concise output
runmode = 2        # 0: user tree;  1: semi-automatic; 2: automatic
                   # 3: StepwiseAddition; (4,5):PerturbationNNI
model = 0          # 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                   # 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
Mgene = 0          # 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
clock = 0          # 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
fix_kappa = 0      # 0: estimate kappa; 1: fix kappa at value below
kappa = 5          # initial or fixed kappa
fix_alpha = 0      # 0: estimate alpha; 1: fix alpha at value below
alpha = 0.5        # initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0         # 1: different alpha's for genes, 0: one alpha
ncatG = 5          # no. of categories in the dG, AdG, or nparK models of rates
nparK = 0          # rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK
nhomo = 0          # 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
getSE = 0          # 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0   # (0,1,2): rates (alpha>0) or ancestral states
Small_Diff = 7e-6 
cleandata = 1      # remove sites with ambiguity data (1:yes, 0:no)?
method = 0         # Optimization method 0: simultaneous; 1: one branch a time
")
} # End paml.baseml.show.default().


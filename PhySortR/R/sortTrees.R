#' Sorts Phylogenetic Trees using Taxa Identifiers
#' 
#' Reads phylogenetic trees from a directory and sorts them based on the presence 
#' of Exclusive and Non-Exclusive clades containing a set of given target leaves at   
#' a desired support value. Can interpret trees in both Newick and extended Newick format.
#' @param target.groups         a set of one or more terms that represent the target leaves whose membership will 
#'                              be tested in each clade during sorting. Multiple terms are to be separated by a comma 
#'                              (\code{"Taxon1,Taxon2"}). This process is case sensitive and uses strict  
#'                              string-matching, so the taxa identifiers must be unique i.e.
#'                              \code{"plantae"} and \code{"Viridiplantae"} might not be appropriate
#'                              as the first is a subset of the second.
#' @param min.support           the minimum support (i.e. between 0-1 or 0-100) of a clade (Default = 0). Support
#'                              values missing from phylogenetic trees are interpreted as zero.
#' @param min.prop.target       the minimum proportion (between 0.0-1.0) of target leaves to be present in a
#'                              clade out of the total target leaves in the tree (Default = 0.7).
#' @param in.dir                directory containing the phylogenetic trees to be sorted (Default = current working directory).
#' @param out.dir               directory to be created within \code{in.dir} for the trees identified during  
#'                              sorting. If \code{out.dir} is omitted the default of \code{Sorted_Trees/}
#'                              will be used. 
#' @param mode                  option to \code{"m"} (move), \code{"c"} (copy) or \code{"l"} (list) trees identified
#'                              during sorting. In \code{"l"} mode (default) a list of the sorted trees is returned, in the 
#'                              \code{"m"} and \code{"c"} modes a list is returned and the identified trees are 
#'                              moved/copied to the \code{out.dir}.
#' @param clades.sorted         option to control if the function will sort for Exclusive (\code{"E"}) and/or 
#'                              Non-Exclusive (\code{"NE"}) clades. Specify both options by comma separation \code{"E,NE"} (Default).
#'                              Exclusive clades are also sorted into a sub-group of All Exclusive trees. 
#' @param extension             the file extension of the tree files to be analyzed (Default = \code{".tre"}).
#' @param clade.exclusivity     the minimum proportion (0.0 <= x < 1.0) of target leaves to interrupting leaves allowed 
#'                              in each non-exclusive clade (Default = 0.9).
#' @return                      Will always return a list containing the names of the trees identified during sorting, 
#'                              irrespective of the \code{mode} argument.
#' @import ape
#' @import phytools
#' @export
#' @examples
#'  ### Load data ###
#'  extdata <- system.file("extdata", package="PhySortR")
#'  file.copy(dir(extdata, full.names = TRUE), ".")
#'  dir.create("Algae_Trees/")
#'  file.copy(dir(extdata, full.names = TRUE), "Algae_Trees/")
#'  
#'  ### Examples ###
#'  # (1) Sorting using 3 target terms, all other parameters default. 
#'  sortTrees(target.groups = "Rhodophyta,Viridiplantae")
#'  
#'  # The function will search in the users current working directory for files 
#'  # with the extension ".tre" and check them (using default min.support, 
#'  # min.prop.target and clade.exclusivity) for Exclusive, All Exclusive or 
#'  # Non-Exclusive clades. A list will be returned with the names of the trees 
#'  # identified during sorting. 
#'  
#'  
#'  
#'  # (2) Sorting with a target directory and an out directory specified.
#'  sortTrees(target.groups = "Rhodophyta,Viridiplantae",
#'    in.dir= "Algae_Trees/", 
#'    out.dir="Sorted_Trees_RVG/", 
#'    mode = "c")
#'    
#'  # The function will search in "Algae_Trees/" for files with the extension
#'  # ".tre" and check them (using default min.support, min.prop.target, 
#'  # clade.exclusivity) for Exclusive, All Exclusive or Non-Exclusive clades. 
#'  # The function will both (a) return a list of the trees identified during 
#'  # sorting and (b) copy the files into their respective subdirectories of
#'  # "Algae_Trees/Sorted_Trees_RVG/Exclusive/", 
#'  # "Algae_Trees/Sorted_Trees_RVG/Exclusive/All_Exclusive/" and 
#'  # "Algae_Trees/Sorted_Trees_RVG/Non_Exclusive/".
#'  
#'  
#'  
#'  # (3) Sorting with in/out directories, min.prop.target and min.support specified.
#'  sortTrees(target.groups = "Rhodophyta,Viridiplantae",
#'    min.prop.target = 0.8,
#'    min.support = 90,
#'    in.dir= "Algae_Trees/",
#'    out.dir="Sorted_Trees_RVG_95/",
#'    mode = "c",
#'    clades.sorted = "NE",
#'    clade.exclusivity = 0.95)
#'    
#'  # The function will search in "Algae_Trees/" for files with the 
#'  # extension ".tre" and check them for only Non-Exclusive clades. 
#'  # A clade will only be defined if it has support >= 90 and contains at least
#'  # 80% of the total target leaves in the tree. A Non-Exclusive clade must also
#'  # be composed of >= 95% target taxa (i.e. < 5% non-target taxa).
#'  # The function will (a) return a list of the trees identified during 
#'  # sorting and (b) copy the trees identified during sorting to the out 
#'  # directory "Algae_Trees/Sorted_Trees_RVG/Non_Exclusive/".
#'  
#'  ### Clean up ###
#'  unlink("Algae_Trees", recursive=TRUE)
#'  unlink("Sorted_Trees.log")
#'  unlink(dir(".", ".*.tre$"))

sortTrees <- function(target.groups, min.support = 0, min.prop.target = 0.7, in.dir = ".", out.dir = "Sorted_Trees", mode = "l", clades.sorted = "E,NE", extension = ".tre", clade.exclusivity = 0.90) {
  ############################
  #### VALIDATE ARGUMENTS ####
  ############################
  # Check status of arguments. Raise fatal error if problem found.
  
  ## Test target.groups argument.
  if (missing(target.groups)) {
    err <- simpleError(call = "target.groups", message =  "ERROR: 'target.groups' is missing!")
    stop(err)
  }
  if (! is.character(target.groups) || is.list(target.groups) || length(target.groups) > 1) {
    err <- simpleError(call = "target.groups", message = "ERROR: 'target.groups' is not a string!")
    stop(err)
  }
  
  
  ## Test min.support argument.
  if (! is.numeric(min.support) || is.list(min.support) || length(min.support) > 1) {
    err <- simpleError(call = "min.support", message = "ERROR: 'min.support' is not an integer!")
    stop(err)
  }
  
  
  ## Test min.prop.target argument.
  if (! is.numeric(min.prop.target) || is.list(min.prop.target) || length(min.prop.target) > 1 || ! min.prop.target <= 1 || ! min.prop.target >= 0) {
    err <- simpleError(call = "min.prop.target", message = "ERROR: 'min.prop.target' is not an integer between 0.0-1.0!")
    stop(err)
  }
  
  
  ## Test in.dir argument.
  if (! is.character(in.dir) || is.list(in.dir) || length(in.dir) > 1) {
    err <- simpleError(call = "in.dir", message = "ERROR: 'in.dir' is not a string!")
    stop(err)
  }
  
  
  ## Test out.dir argument.
  if (! is.character(out.dir) || is.list(out.dir) || length(out.dir) > 1) {
    err <- simpleError(call = "out.dir", message = "ERROR: 'out.dir' is not a string!")
    stop(err)
  }
  
  
  ## Test mode argument.
  if (! is.character(mode) || is.list(mode) || length(mode) > 1 || (! "l" %in% mode && ! "c" %in% mode && ! "m" %in% mode)) {
    err <- simpleError(call = "mode", message = "ERROR: 'mode' not valid; either 'l', 'c' or 'm'")
    stop(err)
  }
  
  
  ## Test clades.sorted argument.
  if (! is.character(clades.sorted) || is.list(clades.sorted) || length(clades.sorted) > 1 || !any(grepl("^E$|^NE$", strsplit(clades.sorted, ",")[[1]]))) {
    err <- simpleError(call = "clades.sorted", message = "ERROR: 'clades.sorted' is not 'E' or 'NE' or 'E,NE'")
    stop(err)
  }
  
  
  ## Test extension argument.
  if (! is.character(extension) || is.list(extension) || length(extension) > 1) {
    err <- simpleError(call = "extension", message = "ERROR: 'extension' is not a string!")
    stop(err)
  }
  
  
  ## Test clade.exclusivity argument.
  if (! is.numeric(clade.exclusivity) || is.list(clade.exclusivity) || length(clade.exclusivity) > 1 || ! clade.exclusivity < 1.0 || ! clade.exclusivity >= 0.0) {
    err <- simpleError(call = "clade.exclusivity", message = "ERROR: 'clade.exclusivity' is not a number between 0.0 <= x < 1.0!")
    stop(err)
  }
  
  
  ##########################
  #### HANDLE VARIABLES ####
  ##########################
  # Split targets by comma
  target.groups <- strsplit(target.groups, ",")[[1]]
  # Remove leading and/or trailing whitespace
  target.groups <- gsub("^\\s+|\\s+$", "", target.groups)
  
  # Create regex for file extension.
  extension.regex <- paste("\\", extension, "$", sep="")
  # Get absolute path for all trees in in.dir.
  file.fullNames <- dir(in.dir, pattern=extension.regex, full.names=TRUE)
  ## Fatal Error if no trees found.
  if (length(file.fullNames) == 0) {
    stop(simpleError(call = "EmptyDirectory", message = "No trees found!"))
  }
  
  ## Handle the different modes
  clades.sorted.split <- strsplit(clades.sorted, ",")[[1]]
  # Remove leading and/or trailing whitespace
  clades.sorted.split <- gsub("^\\s+|\\s+$", "", clades.sorted.split)
  
  # Sort Exclusive trees.
  if ("E" %in% clades.sorted.split) {
    sort.exclusive <- TRUE
  } else {
    sort.exclusive <- FALSE
  }
  # Sort Non-Exclusive trees. 
  if ("NE" %in% clades.sorted.split) {
    sort.nonexclusive <- TRUE
  } else {
    sort.nonexclusive <- FALSE
  }
  
  # To copy (TRUE) or move (FALSE) trees.
  copy.trees <- TRUE
  if (mode == "m") {
    copy.trees <- FALSE
  }
  # To copy/move trees or just return list.
  copy.trees.onoff <- TRUE
  if (mode == "l") {
    copy.trees.onoff <- FALSE
  }
  

  # If trees are to be copied or moved. 
  if (copy.trees.onoff) {
    # Create out.dir directory if specified.
    final.outpath <- in.dir
    final.outpath <- file.path(final.outpath, out.dir)
    # Create directory; warn if already exists.
    tryCatch({dir.create(final.outpath)}, 
             warning = function(war) {warning(simpleWarning(call = "make.out.dir", message = war[1]))},
             error = function(err) {stop(simpleError(call = "make.out.dir", message = err[1]))})
    # Check if directorie created. 
    if (!file.exists(final.outpath)) {
      stop(simpleError(call="make.out.dir", message=paste("Could not create", final.outpath)))
    }
    
    # Create output files, writes WARNING if path already exists.
    # Create Exclusive out locations if required
    if (sort.exclusive) {
      # Exclusive
      tryCatch({dir.create(file.path(final.outpath, "Exclusive"))}, 
               warning = function(war) {warning(simpleWarning(call = "make.out.dir", message = war[1]))},
               error = function(err) {stop(simpleError(call = "make.out.dir", message = err[1]))})
      # Check if directorie created. 
      if (!file.exists(file.path(final.outpath, "Exclusive"))) {
        stop(simpleError(call="make.out.dir", message=paste("Could not create", file.path(final.outpath, "Exclusive"))))
      }
      
      # All-Exclusive
      tryCatch({dir.create(file.path(final.outpath, "Exclusive/All_Exclusive"))}, 
               warning = function(war) {warning(simpleWarning(call = "make.out.dir", message = war[1]))},
               error = function(err) {stop(simpleError(call = "make.out.dir", message = err[1]))})
      # Check if directorie created. 
      if (!file.exists(file.path(final.outpath, "Exclusive/All_Exclusive"))) {
        stop(simpleError(call="make.out.dir", message=paste("Could not create", file.path(final.outpath, "Exclusive/All_Exclusive"))))
      }
    }
    # Create Non-Exclusive out location if required
    if (sort.nonexclusive) {
      tryCatch({dir.create(file.path(final.outpath, "Non_Exclusive"))}, 
               warning = function(war) {warning(simpleWarning(call = "make.out.dir", message = war[1]))},
               error = function(err) {stop(simpleError(call = "make.out.dir", message = err[1]))})
      # Check if directorie created. 
      if (!file.exists(file.path(final.outpath, "Non_Exclusive"))) {
        stop(simpleError(call="make.out.dir", message=paste("Could not create", file.path(final.outpath, "Non_Exclusive"))))
      }
    }
  }
  
  # Open log file.
  LOG <- file(file.path(in.dir, paste(gsub("^/+|/+$", "", out.dir), ".log", sep="")), "w")
  ##########################
  
  # Set variables and counters.
  number.trees <- length(file.fullNames)
  
  # Create names list based off of selected clades.sorted parameters. 
  if (sort.exclusive & sort.nonexclusive) {
    tree.exclusivity <- list("All.Exclusive" = c(), "Exclusive" = c(), "Non.Exclusive" = c())
  } else if (sort.exclusive & ! sort.nonexclusive) {
    tree.exclusivity <- list("All.Exclusive" = c(), "Exclusive" = c())
  } else if (! sort.exclusive & sort.nonexclusive) {
    tree.exclusivity <- list("Non.Exclusive" = c())
  }
  
  # For each tree found in file. 
  for (i in 1:number.trees) {
    # Get name of tree for later use.
    tmp <- strsplit(file.fullNames[i], "/")[[1]]
    tree.name <- tmp[[length(tmp)]]
    # Write header for tree to LOG file.
    write("===========", LOG, append=TRUE)
    write(tree.name, LOG, append=TRUE)
    
    # Read in tree file. 
    tree.txt <- NULL
    possibleError <- tryCatch({
      tree.raw <- file(file.fullNames[i], "r")
      tree.txt <- readLines(tree.raw)
      close(tree.raw)}, 
    warning = function(war) {return(war)},
    error = function(err) {return(err)})
    
    # Move to next tree if file reading failed for any reason.  
    if (inherits(possibleError, "error") | inherits(possibleError, "warning")) {
      print (possibleError)
      warning(simpleError(call="read.tree.file", message=paste(possibleError$call, possibleError$message, "Moving to next file!")))
      next
    }
    # Check if readLines found an empty file (move to next file).
    if (is.null(tree.txt) | length(tree.txt) == 0) {
      warning(simpleError(call = "empty.tree.file", message = paste("File", file.fullNames[i], "appears to be empty. Moving to next file!")))
      next
    }
    
    # Take only 1st line in file. Ignore all other lines. 
    tree.txt <- tree.txt[[1]]
    
    # Check if it is in Newick format or extended Newick format before reading with ape.
    # If extended format, convert with convert.eNewick(). 
    if (grepl("\\[.*\\]", tree.txt)) {
      possibleError <- tryCatch({
        tree.txt <- PhySortR::convert.eNewick(tree.txt)},
        warning = function(war) {return(war)},
        error = function(err) {return(err)})
      
      # Check if convert.eNewick had problems. 
      if (inherits(possibleError, "error") | inherits(possibleError, "warning")) {
        warning(simpleError(call = "convert.eNewick", message = paste(possibleError$call, possibleError$message, "Moving to next file!")))
        next
      }
    }

    # Import tree using ape function. Will return NULL if tree can not be resolved.
    tree <- NULL
    tryCatch({tree <- ape::read.tree(text=tree.txt)}, 
      warning = function(war) {warning(simpleError(call = "ape.read.tree", message = paste(war$call, war$message)))},
      error = function(err) {warning(simpleError(call = "ape.read.tree", message = paste(err$call, err$message, "Moving to next file!"))); next})
    
    # If still NULL or not of class 'phylo' print error message and move on to next tree. 
    # Also check that ape::read.tree has returned an object with the required headers. 
    if (class(tree) != 'phylo' | !all(c("edge", "Nnode", "tip.label") %in% names(tree))) {
      warning(simpleError(call = "ape.read.phylo", message = "ape::read.tree failed to return an object of class 'phylo'. Moving to next file!"))
      next
    }
    
    ##################################
    #### GENERATE DATA STRUCTURES ####
    ##################################
    
    # Get size variables.
    Nnodes <- tree[["Nnode"]]
    all.tips <- tree[["tip.label"]]
    Ntips <- length(all.tips)
    
    # Logic matrix indicating if a particular leaf matches a particular target group.
    #   Each column represents a target terms
    #   Each rows represent a tip from the whole tree
    # Each value of matrix indicates if target matches that tip.
    tip.target.matrix <- sapply(target.groups, function(x) grepl(x, all.tips))
    
    # Get number of target leaves in whole tree over all target.groups
    leaves.match.tree <- sum(tip.target.matrix)
    
    # Get descendants (both nodes and tips) for all nodes in tree. 
    # getDescendants() returns the tip and node numbers used by class "phylo" (APE Package) 
    node.desc <- sapply((Ntips+1):(Ntips+Nnodes), function(x) phytools::getDescendants(tree, x))
    # Take only the tip numbers i.e. numbers <= Ntips
    node.desc.leaves <- sapply(1:Nnodes, function(x) node.desc[[x]][node.desc[[x]] <= Ntips])
    # Logic vector indicating if a tip is in a node (invert for rooted nodes!).
    #   Each column represents a node
    #   Each rows represents a tip
    logic.of.tips <- sapply(1:Nnodes, function(x) c(1:Ntips) %in% node.desc.leaves[[x]])
    
    # List of nodes and there associated tips.
    list.of.tips <- lapply(1:Nnodes, function(x) all.tips[logic.of.tips[ ,x]])
    list.of.tips.rev <- lapply(1:Nnodes, function(x) all.tips[!logic.of.tips[ ,x]])
    
    # Number of targets in each node.
    number.targets.in.nodes <- sapply(1:Nnodes, function(x) sum(logic.of.tips[ ,x] & tip.target.matrix))
    # Number of target in rooted nodes is the difference between the total 
    # targets in the tree and the targets in the unrooted node.
    number.targets.in.nodes.rev <- c(leaves.match.tree - number.targets.in.nodes)
    
    # Returns a logic vector with each element indicating if all the target.groups have been found in that node.
    all.targets.in.nodes <- sapply(1:Nnodes, function(x) all(colSums(logic.of.tips[ ,x] & tip.target.matrix) >=1))
    # Invert logic.of.tips for rooted nodes.
    all.targets.in.nodes.rev <- sapply(1:Nnodes, function(x) all(colSums(!logic.of.tips[ ,x] & tip.target.matrix) >=1))
    
    # Number of leaves in each node.
    number.leaves.in.nodes <- sapply(1:Nnodes, function(x) length(list.of.tips[[x]]))
    # Rooted nodes is the difference between the unrooted nodes and the total tips in the tree.
    number.leaves.in.nodes.rev <- c(Ntips - number.leaves.in.nodes)
    
    # Check if ape returned 'node.label' list. 
    if (!"node.label" %in% names(tree)) {
      # Assume if 'node.label' is absent make nodes zero. 
      tree.nodes.support.numeric <- rep(0, Nnodes)
    } else {
      # Get vector with each entry as support for that node.
      tree.nodes.support.numeric <- c((as.numeric(tree$node.label[1:Nnodes])))
      # Convert 'NA' values to zero. 
      tree.nodes.support.numeric[is.na(tree.nodes.support.numeric)] <- 0
    } 
    
    # Gets logic vector of nodes with support above threshold.
    tree.nodes.support <- (!is.na(tree.nodes.support.numeric[1:Nnodes]) & 
                             tree.nodes.support.numeric[1:Nnodes] >= min.support)
    
    # Screen both the unrooted and rooted nodes.
    #   If node has support >= min.support &
    #     If all targets represented in node &
    #       (Total target leaves in node / total target leaves in whole tree) >= cut off
    rules.check <- c(tree.nodes.support & all.targets.in.nodes & 
                       ((number.targets.in.nodes / leaves.match.tree) >= min.prop.target))
    rules.check.rev <- c(tree.nodes.support & all.targets.in.nodes.rev & 
                           ((number.targets.in.nodes.rev / leaves.match.tree) >= min.prop.target))
    
    # Set final variables.
    whole.tree <- FALSE
    all.queries.in.one.clade.only <- FALSE
    exclusive.monophyly.found <- FALSE
    nonexclusive.monophyly.found <- FALSE
    
    # Check if whole tree (first node) is populated entirely by targets.
    if (number.targets.in.nodes[[1]] == length(list.of.tips[[1]]) & all.targets.in.nodes[[1]]) {
      whole.tree <- TRUE
      all.queries.in.one.clade.only <- TRUE
      exclusive.monophyly.found <- TRUE
      if (sort.exclusive) {
        write(paste("** Whole Tree; Targetgroup in tree:", 
                    number.targets.in.nodes[[1]]), LOG, append=TRUE)
        write(list.of.tips[[1]], LOG, append=TRUE)
      }
    }
    
    #########################
    #### CHECH EACH NODE ####
    #########################
    # Store extra output information for Non-Exclusive trees.
    nelines <- character()
    
    # Check each node for compliance with cut-offs.
    for (x in (1:Nnodes)[rules.check | rules.check.rev]) {
      if (rules.check[[x]]) {
        # If all leaves in node are queries
        if (number.targets.in.nodes[[x]] == number.leaves.in.nodes[[x]]) {
          if (sort.exclusive) {
            all.queries.in.one.clade.only <- TRUE
            exclusive.monophyly.found <- TRUE
            write(paste("** Support: ", tree.nodes.support.numeric[[x]], 
                        "; Targetgroup in clade: ", number.targets.in.nodes[[x]], 
                        "; Targetgroup in tree: ", leaves.match.tree, sep=""), LOG, append=TRUE)
            write(list.of.tips[[x]], LOG, append=TRUE)
          }
          # If NOT all leaves in node are queries.
        } else {
          all.queries.in.one.clade.only <- TRUE
          if ((number.targets.in.nodes[x] / number.leaves.in.nodes[[x]]) >= clade.exclusivity) {
            if (sort.nonexclusive) {
              nonexclusive.monophyly.found <- TRUE
              write(paste("** NE Support: ", tree.nodes.support.numeric[[x]], 
                          "; Targetgroup in clade: ", number.targets.in.nodes[[x]], 
                          "; Species in clade: ", number.leaves.in.nodes[[x]], 
                          "; Total in tree: ", leaves.match.tree, sep=""), LOG, append=TRUE)
              write(list.of.tips[[x]], LOG, append=TRUE)
              nelines <- c(nelines, paste("** NE Support: ", tree.nodes.support.numeric[[x]], 
                                          "; Targetgroup in clade: ", number.targets.in.nodes[[x]], 
                                          "; Species in clade: ", number.leaves.in.nodes[[x]], 
                                          "; Total in tree: ", leaves.match.tree, sep=""))
              nelines <- c(nelines, list.of.tips[[x]])
            }
          }
        }
      }
      if (rules.check.rev[[x]]) {
        # If all leaves in node are queries.
        if (number.targets.in.nodes.rev[[x]] == number.leaves.in.nodes.rev[[x]]) {
          if (sort.exclusive) {
            all.queries.in.one.clade.only <- TRUE
            exclusive.monophyly.found <- TRUE
            write(paste("** Support: ", tree.nodes.support.numeric[[x]], 
                        "; Targetgroup in clade: ", number.targets.in.nodes.rev[[x]], 
                        "; Targetgroup in tree: ", leaves.match.tree, sep=""), LOG, append=TRUE)
            write(list.of.tips.rev[[x]], LOG, append=TRUE)
          }
          
          # If NOT all leaves in node are queries.
        } else {
          all.queries.in.one.clade.only <- TRUE
          if ((number.targets.in.nodes.rev[[x]] / number.leaves.in.nodes.rev[[x]]) >= clade.exclusivity) {
            if (sort.nonexclusive) {
              nonexclusive.monophyly.found <- TRUE
              write(paste("** NE Support: ", tree.nodes.support.numeric[[x]], 
                          "; Targetgroup in clade: ", number.targets.in.nodes.rev[[x]], 
                          "; Species in clade: ", number.leaves.in.nodes.rev[[x]], 
                          "; Total in tree: ", leaves.match.tree, sep=""), LOG, append=TRUE)
              write(list.of.tips.rev[[x]], LOG, append=TRUE)
              nelines <- c(nelines, paste("** NE Support: ", tree.nodes.support.numeric[[x]], 
                                          "; Targetgroup in clade: ", number.targets.in.nodes.rev[[x]], 
                                          "; Species in clade: ", number.leaves.in.nodes.rev[[x]], 
                                          "; Total in tree: ", leaves.match.tree, sep=""))
              nelines <- c(nelines, list.of.tips.rev[[x]])
            }
          }
        }
      }
    }
    
    ################
    #### RETURN ####
    ################
    # If all query leaves have been found in ONLY one node.
    if (all.queries.in.one.clade.only) {
      # If NOT all leaves in node are queries.
      if (exclusive.monophyly.found & sort.exclusive) {
        # If monophyly is whole tree.
        if (whole.tree) {
          # If returning trees with exclusive clades.  
          if (copy.trees.onoff) {
            # Action to take.
            if (copy.trees) {
              file.copy(file.fullNames[i], file.path(final.outpath, "Exclusive/All_Exclusive"))
            } else {
              file.copy(file.fullNames[i], file.path(final.outpath, "Exclusive/All_Exclusive"))
              file.remove(file.fullNames[i])
            }
          }
          tree.exclusivity$All.Exclusive <- c(tree.exclusivity$All.Exclusive, tree.name)
          
          # If monophyly is not whole tree.
        } else {
          # If returning trees with exclusive clades.  
          if (copy.trees.onoff) {
            # Action to take.
            if (copy.trees) {
              file.copy(file.fullNames[i], file.path(final.outpath, "Exclusive"))
            } else {
              file.copy(file.fullNames[i], file.path(final.outpath, "Exclusive"))
              file.remove(file.fullNames[i])
            }
          }
          tree.exclusivity$Exclusive <- c(tree.exclusivity$Exclusive, tree.name)
        }
        
        # If MOST leaves in node are queries and is not exclusive monophyly.
      } else if (! exclusive.monophyly.found & nonexclusive.monophyly.found & sort.nonexclusive) {
        # If returning trees with non-exclusive clades.  
        if (copy.trees.onoff) {
          # Write nelines to file.
          NEOUT <- file(file.path(final.outpath, "Non_Exclusive", paste(tree.name, ".sortgroup.txt", sep='')), "w")
          writeLines(nelines, NEOUT, sep="\n")
          close(NEOUT)
          # Action to take.
          if (copy.trees) {
            file.copy(file.fullNames[i], file.path(final.outpath, "Non_Exclusive"))
          } else {
            file.copy(file.fullNames[i], file.path(final.outpath, "Non_Exclusive"))
            file.remove(file.fullNames[i])
          }
        }
        tree.exclusivity$Non.Exclusive <- c(tree.exclusivity$Non.Exclusive, tree.name)
        
      # Else if all queries not in one clade only.
      } else {
        write("Negative.", LOG, append=TRUE)
      }
    } else {
      write("Negative.", LOG, append=TRUE)
    }
  }
  # Close log file
  close(LOG)
  
  return (tree.exclusivity)
}

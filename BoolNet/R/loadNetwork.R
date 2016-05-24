# Load a network in a specified rule description language
# from file <file>.
# <bodySeparator> is the character that separates targets and factors
# <lowercaseGenes> specifies whether gene names are converted to lower case
loadNetwork <- function(file, bodySeparator=",", lowercaseGenes=FALSE, symbolic=FALSE) 
{
  func <- readLines(file,-1)
  
  # strip comments
  func <- gsub("#.*", "", trim(func))
  func <- func[nchar(func) > 0]
  
  # check header
  if (length(func) == 0)
    stop("Header expected!")
  
  header <- func[1]
  header <- tolower(trim(strsplit(header, bodySeparator)[[1]]))
  
  if (length(header) < 2 || header[1] != "targets" || !(header[2] %in% c("functions","factors")) ||
      (length(header) == 3 && header[3] != "probabilities"))
      stop(paste("Invalid header:", func[1]))

  func <- func[-1]      
  
  if (lowercaseGenes)
    func <- tolower(func)
  
  # Replace all invalid characters to be able to load most networks
  func <- gsub("[^\\[\\]a-zA-Z0-9_\\|\\&!\\(\\) \t\\-+=.,]+","_", func, perl=TRUE)

  tmp <-  unname(lapply(func,function(x)
  # split strings at separators that are NOT
  # inside a bracket block
  {
    bracketCount <- 0
    lastIdx <- 1

    chars <- strsplit(x,split="")[[1]]
    
    res <- c()
    
    if (length(chars) > 0)
    {
      for (i in seq_along(chars))
      {
        if (chars[i] == "(")
          bracketCount <- bracketCount + 1
        else
        if (chars[i] == ")")
          bracketCount <- bracketCount -1
        else
        if (chars[i] == bodySeparator && bracketCount == 0)
        {
          res <- c(res, trim(paste(chars[lastIdx:(i-1)],collapse="")))
          lastIdx <- i + 1
        }
      }
      res <- c(res, trim(paste(chars[lastIdx:length(chars)],collapse="")))
    }
    return(res)
  }))
  targets <- sapply(tmp,function(rule)trim(rule[1]))
  
  for (target in targets)
  {
    if (regexec("^[a-zA-Z_][a-zA-Z0-9_]*$", target)[[1]] == -1)
      stop(paste("Invalid gene name:",target))
  }
   
  
  factors <- sapply(tmp,function(rule)trim(rule[2]))
  
  temporal <- length(grep("timeis|timelt|timegt|\\[|\\]", factors, ignore.case=TRUE) > 0)
  if (temporal && !symbolic)
  {
    warning("The network contains temporal elements. This requires loading the model with symbolic=TRUE!")
    symbolic <- TRUE
  }
  
  
  probabilities <- sapply(tmp,function(rule)
                          {
                            if (length(rule) >= 3)
                              as.numeric(trim(rule[3]))
                            else
                              1.0
                          })
  
  factors.tmp <- lapply(factors,matchNames)

  # create list of all gene names in both sides of the functions
  genes <- unique(c(targets,unname(unlist(factors.tmp))))
  isProbabilistic <- (length(unique(targets)) < length(targets))
                          
  if (symbolic)
  {
    if (isProbabilistic)
      stop("Probabilistic networks cannot be loaded with symbolic=TRUE!")
                                
    interactions <- lapply(factors, function(rule)parseBooleanFunction(rule, genes))
    
    onlyInputs <- setdiff(genes,targets)
    if(length(onlyInputs) > 0)
    # some genes are only used as inputs, but are not involved in the network
    # -> create dummy input and function
    {
      for(gene in onlyInputs)
      {
        warning(paste("There is no transition function for gene \"",
                       gene,"\"! Assuming an input!",sep=""))
        
        interactions[[gene]] = parseBooleanFunction(gene, genes)
      }
    }
    
    delays <- apply(sapply(interactions,maxTimeDelay,genes=genes),1,max)
    names(interactions) <- genes
    
    fixed <- as.integer(sapply(interactions, function(int)
                        {
                          if (int$type == "constant")
                            int$value
                          else
                            -1L  
                        }))
    names(fixed) <- genes
    
    res <- list(genes = genes, interactions=interactions, fixed=fixed)
    res$internalStructs <- .Call("constructNetworkTrees_R",res)
    res$timeDelays <- delays
    
    class(res) <- "SymbolicBooleanNetwork"
    return(res)
  }
  else
  {
    # extract "real" gene names from the list, drop constants
    #suppressWarnings(genes <- genes[is.na(as.integer(genes))])

    fixed <- rep(-1,length(genes))
    names(fixed) <- genes

    interactions <- list()

    for(i in seq_along(targets))
    {
      target <- targets[i]
      interaction <- generateInteraction(factors[i], genes);

      if (isProbabilistic)
      {
       interaction$probability <- probabilities[i]
       interactions[[target]][[length(interactions[[target]]) + 1]] <- interaction
      }
      else
      {
        if (length(interaction$func) == 1)
        # this is a constant gene => fix it
        {
          fixed[target] <- interaction$func
        }
        interactions[[target]] <- interaction
      }
    }

    onlyInputs <- setdiff(genes,targets)
    if(length(onlyInputs) > 0)
    # some genes are only used as inputs, but are not involved in the network
    # -> create dummy input and function
    {
      for(gene in onlyInputs)
      {
        warning(paste("There is no transition function for gene \"",
                       gene,"\"! Assuming an input!",sep=""))
        if (isProbabilistic)
          interactions[[gene]][[1]] = list(input = length(interactions)+1,func=c(0,1),
          									  expression = gene, probability=1.0)
        else
          interactions[[gene]] = list(input = length(interactions)+1,func=c(0,1),
          							expression = gene)
      }
    }
    
    if (isProbabilistic)
    {
      wrongProb <- sapply(interactions,function(interaction)
                   {
                      abs(1.0-sum(sapply(interaction,function(func)func$probability))) > 0.0001
                   })
      if (any(wrongProb))
        stop(paste("The probabilities of gene(s) ",paste(genes[wrongProb],collapse=", ")," do not sum up to 1!",sep=""))
    }  

    res <- list(interactions = interactions,
                genes = genes,
                fixed = fixed)
    
    if (isProbabilistic)
      class(res) <- "ProbabilisticBooleanNetwork"
    else
      class(res) <- "BooleanNetwork"
    return(res)

  }
}

matchNames <- function(rule)
{
  regexpr <- "([_a-zA-Z][_a-zA-Z0-9]*)[,| |\\)|\\||\\&|\\[]"
  rule <- paste(gsub(" ", "", rule, fixed=TRUE)," ",sep="")
  res <- unique(unname(sapply(regmatches(rule,gregexpr(regexpr, rule))[[1]], 
                    function(m) 
                    {
                      sapply(regmatches(m,regexec(regexpr,m)),function(x)x[2])
                    })))
  
  # remove operators
  isOp <- sapply(res, function(x)
          {
            tolower(x) %in% c("all", "any", 
                              "sumis", "sumgt", "sumlt", 
                              "maj", "timegt", "timelt", "timeis")
          })
                      
  return(res[!isOp])
}



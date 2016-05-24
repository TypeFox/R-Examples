dctrl <-
function(d,correct,items)
{
  problemlist <- list()  
  a <- 1  
  probl <- FALSE
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ dataframe
  if(!is.data.frame(d))
  {
    problemlist[[a]] <- "data object is not a data.frame!"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ minimum
  minimum <- sapply(d,function(x)min(as.numeric(x),na.rm=T) == 0)
#   if(!any(minimum))
#   {
#     problemlist[[a]] <- "minimum != 1 in at least one column!"
#     a <- a + 1
#     probl <- TRUE
#   }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ number of cat
  
  un <- sapply(d,function(x)length(table(x))> 1)
  if(!any(un))
  {
    problemlist[[a]] <- "there are variables containing only 1s!"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ correct_1
  
  leco <- length(correct) == ncol(d)
  if(!leco)
  {
    problemlist[[a]] <- "length(correct) != ncol(d)"
    a <- a + 1
    probl <- TRUE
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~ correct_2
  
  leco2 <- mapply(function(A,B) B %in% A ,A=d,B=correct)
  if(!any(leco2))
  {
    problemlist[[a]] <- "the category of correct response was not observed (in one or more variables)!"
    probl <- TRUE
  }
  
  

  list(problemlist=problemlist,probl=probl)
}



################## more controls #######################


dctrl2 <-
  function(d_new,da,items)
  {
    
    
# checks the number of obs in each category!    
    testdim <- lapply(d_new,function(x)
      {
        apply(x,2,function(aa3) sort(unique(aa3)))  
      })
    
    if(!all(sapply(testdim[-1],function(td) identical(td,testdim[[1]])))){stop("At least in one item there are categories with no observations in one of the groups!")}
    

    ## all successive
    if(!all(sapply(da[,items],function(zzq) all(diff(sort(unique(zzq)),lag=1)==1) ))){stop("Category levels must be consecutive numbers!")} # thx to CH for nice error message! see pcIRT :-)
      


  }


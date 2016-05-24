summary.catIrt <-
function(object, group = TRUE, ids = "none", ... ){
  
# For the summary we want to calculate theta/SE/Items stuff #
  x <- object
  
# First, making sure that if "ids" is all, we should include all summaries.
  if( any( ids == "all" ) )
    ids <- seq_along(x$cat_theta)
    
# Next, checking the "ids" argument.
  if( any( ids != "none" & ids != "all" ) )
    if( missing(ids) | ( any( ids < 1 | ids > length(x$cat_theta) ) ) )
      stop( "ids must be 'all', 'none', or an integer vector containing 1 -- length(thetas)" )
  
# We might need to reclass the model after removing the first column.
  full_params <- x$full_params
  
## Doing basic extraction and estimation ##
  cat_theta  <- x$cat_theta
  cat_sem    <- x$cat_sem
  tot_theta  <- x$tot_theta
  tot_sem    <- x$tot_sem
  true_theta <- x$true_theta
    
  cat_length <- x$cat_length
    

###############
# Theta Stuff #
###############

## Doing basic calculation on theta ##
  tot_theta.mn   <- mean(tot_theta, na.rm = TRUE)
  tot_theta.sd   <- sd(tot_theta, na.rm = TRUE)
  tot_theta.min  <- unique(min(tot_theta, na.rm = TRUE))
  tot_theta.max  <- unique(max(tot_theta, na.rm = TRUE))
  cat_theta.mn   <- mean(cat_theta, na.rm = TRUE)
  cat_theta.sd   <- sd(cat_theta, na.rm = TRUE)
  cat_theta.min  <- unique(min(cat_theta, na.rm = TRUE))
  cat_theta.max  <- unique(max(cat_theta, na.rm = TRUE))
  
  if( length(true_theta) >= 1 & any( !is.na(true_theta) ) ){
    true_theta.mn  <- mean(true_theta, na.rm = TRUE)
    true_theta.sd  <- sd(true_theta, na.rm = TRUE)
    true_theta.min <- unique(min(true_theta, na.rm = TRUE))
    true_theta.max <- unique(max(true_theta, na.rm = TRUE))
  } else{
    true_theta.mn <- true_theta.sd <- true_theta.min <- true_theta.max <- NULL
  } # END ifelse STATEMENT
  
# Combining those calculations into a matrix.
  theta <- list()
  theta$sum <- matrix( c(tot_theta.mn, tot_theta.sd, tot_theta.min, tot_theta.max,
                         cat_theta.mn, cat_theta.sd, cat_theta.min, cat_theta.max),
                       nrow = 2, byrow = TRUE )
  theta$sum <- rbind( c(true_theta.mn, true_theta.sd, true_theta.min, true_theta.max),
                      theta$sum )
                    
  colnames(theta$sum) <- c("mean", "sd", "min", "max")
  rownames(theta$sum) <- { if( !is.null(true_theta.mn) ){
  	                         c("true", "tot", "cat")
  	                       } else{
  	                         c("tot", "cat")
  	                       } # END ifelse STATEMENT
  	                     }
  	                         
  
# Doing comparison calculations on cat_theta --> tot_theta.
  theta.n      <- length(cat_theta)
  theta.cord   <- cor(cat_theta, tot_theta, use = "complete.obs")^2
  theta.asd    <- mean(tot_theta - cat_theta, na.rm = TRUE)
  theta.aad    <- mean(abs(tot_theta - cat_theta),   na.rm = TRUE )
  theta.msd    <- mean(   (tot_theta - cat_theta)^2, na.rm = TRUE )
  theta.rmsd   <- sqrt(theta.msd)
  
# Combining those calculations into a vector.
  theta$d_stat <- c(n = theta.n, rsq = theta.cord, asd = theta.asd,
                    aad = theta.aad, msd = theta.msd, rmsd = theta.rmsd)
  
# Doing comparison calculations on cat_theta --> true_theta.
  if( length(true_theta) >= 1 & any( !is.na(true_theta) ) ){
  	
    theta.core   <- cor(cat_theta, true_theta, use = "complete.obs")^2
    theta.ase    <- mean(true_theta - cat_theta, na.rm = TRUE)
    theta.aae    <- mean(abs(true_theta - cat_theta),   na.rm = TRUE )
    theta.mse    <- mean(   (true_theta - cat_theta)^2, na.rm = TRUE )
    theta.rmse   <- sqrt(theta.mse)
    
    theta$e_stat <- c(n = theta.n, rsq = theta.core, ase = theta.ase,
                      aae = theta.aae, mse = theta.mse, rmse = theta.rmse)
                          
  } # END if STATEMENT

    
#############
# SEM Stuff #
#############

## Doing basic calculation on SEM ##
  tot_sem.mn  <- mean(tot_sem, na.rm = TRUE)
  tot_sem.sd  <- sd(tot_sem, na.rm = TRUE)
  tot_sem.min <- unique(min(tot_sem, na.rm = TRUE))
  tot_sem.max <- unique(max(tot_sem, na.rm = TRUE))
  cat_sem.mn  <- mean(cat_sem, na.rm = TRUE)
  cat_sem.sd  <- sd(cat_sem, na.rm = TRUE)
  cat_sem.min <- unique(min(cat_sem, na.rm = TRUE))
  cat_sem.max <- unique(max(cat_sem, na.rm = TRUE))
  
# Combining those calculations into a matrix.
  sem <- list()
  sem$sum <- matrix( c(tot_sem.mn, tot_sem.sd, tot_sem.min, tot_sem.max,
                       cat_sem.mn, cat_sem.sd, cat_sem.min, cat_sem.max),
                     nrow = 2, byrow = TRUE )
                    
  colnames(sem$sum) <- c("mean", "sd", "min", "max")
  rownames(sem$sum) <- c("tot", "cat")


# Doing comparison calculations on sem.
  sem.cord <- cor(cat_sem, tot_sem, use = "complete.obs")^2
  sem.asd  <- mean(tot_sem - cat_sem, na.rm = TRUE)
  sem.aad  <- mean(abs(tot_sem - cat_sem),   na.rm = TRUE)
  sem.msd  <- mean(   (tot_sem - cat_sem)^2, na.rm = TRUE )
  sem.rmsd <- sqrt(sem.msd)
  
# Combining those calculations into a vector.
  sem$d_stat <- c(n = theta.n, rsq = sem.cord, asd = sem.asd,
                  aad = sem.aad, msd = sem.msd, rmsd = sem.rmsd)
  
###############
# Items Stuff #
###############

# Doing calculations on items.
  cat_length.mn  <- mean(cat_length, na.rm = TRUE)
  cat_length.sd  <- sd(cat_length, na.rm = TRUE)
  cat_length.min <- unique(min(cat_length, na.rm = TRUE))
  cat_length.max <- unique(max(cat_length, na.rm = TRUE))
  
# Combining those items into a list.
  items <- list()
  items$sum <- c(mean = cat_length.mn, sd = cat_length.sd, min = cat_length.min, max = cat_length.max)
  
  items$bank <- c(n.items = dim(x$full_params)[1], n.params = dim(x$full_params)[2] - 2)
  
##########################
# Total Bank Information #
##########################
  info <- list()
  
# Figuring out the realistic minimum and maximum of our data:

# First: find the min and max of the theta.
  min_theta <- min( unique(min(x$tot_theta)), unique(min(x$full_params[ , -1])) )
  max_theta <- max( unique(max(x$tot_theta)), unique(max(x$full_params[ , -1])) )
  
# Second: build a sequence from the lowest to the highest.
  info$eval <- seq(min_theta, max_theta, length.out = 1000)
  
# Third: pull out parameters and indicate class of those parameters.
  full_params        <- x$full_params[ , -c(1, ncol(x$full_params))]
  class(full_params) <- c(x$mod$mod, "matrix")

# Fourth: calculate info on those parameters across all people.
  y            <- FI(params = full_params, theta  = info$eval, type = "expected")
  info$info    <- y$test
  info$sem     <- y$sem

# And returning values so people can summarize for themselves:
  ans <- list( theta = theta, sem = sem, items = items,
               freq  = list(categ = table(as.character(x$cat_categ)), term = table(x$cat_term)),
               bank_info  = info,
               group = group, ids = ids,
               cat_run = x )
               
  if(length(ans$categ) == 0)
    ans$categ <- NULL
  
  class(ans) <- "summary.catIrt"
  
  return(ans)
  
} # END summary.catIRT FUNCTION
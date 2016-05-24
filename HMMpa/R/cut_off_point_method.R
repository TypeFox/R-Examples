cut_off_point_method <-
function(x, cut_points, names_activity_ranges = NA, hidden_PA_levels = NA, bout_lengths = NULL, plotting= 0) 
{   

################################################################################################################################################################################################################################# Check on arguments ##########################################################################################################
################################################################################################################################################################################  	
 if (is.null(bout_lengths))
 {
 	stop("Set variable 'bout_lengths' to use this function. See help-manual for further information. For example: bout_lengths=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,12,13,20,21,40,41,60,61,80,81,120,121,240,241,480,481,1440,1,1440)")
 }
	
	
################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  	
  cut_points <- c(0, cut_points)
  number_of_activity_ranges <- length(cut_points)
  
  activity_ranges <- array(data = NA, dim = c(1,2,number_of_activity_ranges))
  colnames(activity_ranges) <- c("[lower boundary,","upper boundary)")
  
  ############### Function for bout-extraction (from the result of the classification) #######################################################################################

  bouts <- function(counts, sequence_activity_range, cut_points, bout_lengthss)
  {
    
    for_plotting_a <- sequence_activity_range
    
    for (i in 1:length(cut_points))
    {
      for_plotting_a[for_plotting_a == i] <- cut_points[i]
    }
    for_plotting_a[sequence_activity_range == (length(cut_points) + 1)] <- max(counts)
    
    
    for_plotting_b <- sequence_activity_range
    for_plotting_b <- for_plotting_b - 1
    for (i in 1:length(cut_points))
    {
      for_plotting_b[for_plotting_b == i] <- cut_points[i]
    }
    for_plotting_b[for_plotting_b == 0] <- 0
    
    
    
    dimension_columns_of_next_arry <- max(length(bout_lengthss), length(x), max(bout_lengthss))
    if (dimension_columns_of_next_arry == Inf)
    {
      dimension_columns_of_next_arry <- max(bout_lengthss[!bout_lengthss == Inf], length(x))
    }
    
    bout_periods <- array(0, dim = c((length(cut_points) + 1), dimension_columns_of_next_arry, 1))
    
    if (all(!is.na(names_activity_ranges))) 
    {
      dimnames(bout_periods) <- list(activity_range = names_activity_ranges, Countbouts = 1: dimension_columns_of_next_arry)
    } else {
    dimnames(bout_periods) <- list(activity_range = 1:(length(cut_points) + 1), Countbouts = 1: dimension_columns_of_next_arry)
    }
    temp_sequence_activity_range <- c(sequence_activity_range, -Inf)
    
    for (j in 1:(length(cut_points) + 1))
    {
      temp_vec <- 0
      for (i in 1:length(temp_sequence_activity_range)) 
      {
        if (temp_sequence_activity_range[i] == j) 
        {
          temp_vec <- temp_vec + 1
        } else {
          temp_which_bout <- temp_vec
          bout_periods[j,temp_which_bout,1] <- bout_periods[j,temp_which_bout,1] + 1
          
          temp_vec <- 0  
        }
        
        if (temp_vec == length(temp_sequence_activity_range))
        {
          temp_which_bout <- temp_vec
          bout_periods[j,temp_which_bout,1] = bout_periods[j,temp_which_bout,1] + 1
          
          temp_vec <- 0
        }
      }
    }
     
    quantity_of_bouts <- sum(bout_periods)
    
    return(list(quantity_of_bouts = quantity_of_bouts,
                bout_periods = bout_periods,
                for_plotting_a = for_plotting_a, 
                for_plotting_b = for_plotting_b)) 
  }

################################################################################################################################################################################################################################# Assigning an activity range to the time-series of counts ####################################################################
################################################################################################################################################################################
  	  	  
  for (i in 1:number_of_activity_ranges)
  {
    activity_ranges[1,1,i] <- cut_points[i]
    activity_ranges[1,2,i] <- cut_points[i+1]
  }
  activity_ranges[1,2,number_of_activity_ranges] <- Inf
  
  classification <- rep(NA, times = length(x))
  
  if (any(is.na(hidden_PA_levels)))
  {
    for (i in 1:length(x))
    {	
      for (j in 1:number_of_activity_ranges)
      {
        if (x[i] >= activity_ranges[1,1,j] & x[i] < activity_ranges[1,2,j])
        {	
          classification[i] <- c(j)
        }
      }
    }
    classification_per_activity_range <- pairlist(NA)
    for (i in 1:number_of_activity_ranges)
    {
      classification_per_activity_range[[i]] <- x
      for (j in 1:length(classification_per_activity_range[[i]]))
      {
        if (!classification[j] == i) 
        {
          classification_per_activity_range[[i]][j] <- NA
        }
      }
    }
  }
  ############### If (hidden) Pa-levels are given, then the classification is conducted via the total magnitude of the hidden PA-level of a count ##############################
  if (any(!is.na(hidden_PA_levels)))
  {   
    for (i in 1:length(x))
    {	
      for (j in 1:number_of_activity_ranges)
      {
        if (hidden_PA_levels[i] >= activity_ranges[1,1,j] & hidden_PA_levels[i] < activity_ranges[1,2,j])
        {	
          classification[i] <- c(j)
        }
      }
    }
    classification_per_activity_range <- pairlist(NA)
    for (i in 1:number_of_activity_ranges)
    {
      classification_per_activity_range[[i]] <- x
      for (j in 1:length(classification_per_activity_range[[i]]))
      {
        if (!classification[j] == i)
        {
          classification_per_activity_range[[i]][j] <- NA
        }
      }
    }
  }

################################################################################################################################################################################################################################# Statistics of the classification        #####################################################################################
################################################################################################################################################################################
   
  freq_acitvity_range <- table(c(seq(1, number_of_activity_ranges), classification))
  freq_acitvity_range <- freq_acitvity_range - 1
  rel_freq_acitvity_range <- freq_acitvity_range / sum(table(classification))
  
  if (all(!is.na(names_activity_ranges)))
  {
    rownames(freq_acitvity_range) <- names_activity_ranges
    rownames(rel_freq_acitvity_range) <- names_activity_ranges
    colnames(activity_ranges[1,,]) <- names_activity_ranges
    names(classification_per_activity_range ) <- names_activity_ranges
  }
  
  ############### Bout-extraction (from the result of the classification) ###################################################################################################

  bout_extraction <- bouts(counts = x, sequence_activity_range = classification, cut_points = cut_points[c(-1)], bout_lengthss = bout_lengths)
  
  quantity_of_bouts <- bout_extraction$quantity_of_bouts
  
  bout_periods <- bout_extraction$bout_periods
  
  
  bout_classes <- matrix(bout_lengths, nrow = 2, byrow = FALSE)
  names = NULL
  for (i in seq(1, dim(bout_classes)[2]))
  {
    if (bout_classes[1,i] == bout_classes[2,i])
    {
      names <- c(names, toString(bout_classes[1,i]))
    } else {
      names <- c(names, paste(toString(bout_classes[1,i]), '-', toString(bout_classes[2,i])) )
    }
  }
  #	if (bout_lengths[[length(bout_lengths)]] == Inf)
  #	{
  #		names[length(bout_lengths) / 2] <- paste('>',toString(bout_lengths[length(bout_lengths) - 2]))
  #	}
  colnames(bout_classes) <- names
  rownames(bout_classes) <- c("from", "to")
  
  
  abs_freq_bouts_el <- rep(0, times = length(names))
  for (i in 1:(length(names)))
  {
    abs_freq_bouts_el[i] <- sum(bout_periods[,seq(bout_classes[1,i], bout_classes[2,i]),1])
  }
  

  
  names(abs_freq_bouts_el) <- names
  
  abs_freq_bouts_el <- pairlist(all_activity_ranges = abs_freq_bouts_el)
  
  for (j in 1: number_of_activity_ranges)	    
  {	 
    abs <- rep(0, times = length(names))
    for (i in 1:(length(names)))
    {
      abs[i] <- sum(bout_periods[j,seq(bout_classes[1,i], bout_classes[2,i]),1])
    }
      
    names(abs) <- names
    
    if (all(!is.na(names_activity_ranges))) 
    {
      abs_freq_bouts_el[[names_activity_ranges[j]]] <- abs	
    } else {
      abs_freq_bouts_el[[j+1]] <- abs	
    }
  }
  
################################################################################################################################################################################################################################# Plotting        ###########################################################################################################
################################################################################################################################################################################  
  
  if (!is.na(plotting))
  {	
  	############### output plot 1-5 ########################################################################################################################################## 
    if (plotting == 0)
    {   
      par(mfrow = c(2,2))
      
      if (any(!is.na(hidden_PA_levels))) 
      {
      	plot(x, main = "classification (HMM)", xlab = "time", ylab= "counts", col = "white")
      } else {
         plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	
      
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i], '(', round(100 * rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i, '(', round(100 * rel_freq_acitvity_range[i], digits=2), "%",')')))
        }
      }	
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      } else {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }	
  
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
   	  		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }		      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      if (any(!is.na(hidden_PA_levels)))
      {
     		 barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[1]], names.arg = names, main = "bouts in relation to their length (HMM)")      
      } else {
     		 barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[1]], names.arg = names, main = "bouts in relation to their length")      
      }		  
      

   
      par(mfrow = c(1,1))
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }		 
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }
      
          
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      
      par(mfrow = c(2,1))
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i]  , '(', round(100 * rel_freq_acitvity_range[i], digits=2), "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '(', round(100 * rel_freq_acitvity_range[i], digits=2), "%", ')')))
        }
      }		
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    	barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      } else {
      			barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i]  , '( # ', round(freq_acitvity_range[i], digits=2), ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '( # ', round( freq_acitvity_range[i], digits=2), ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		barplot(freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "absolute frequency", names.arg = foo)	
      } else {
      		barplot(freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "absolute frequency", names.arg = foo)		
      }		
      
      
      par(mfrow = c(number_of_activity_ranges / 2 + 1, 2))
      
      for (i in 1:(number_of_activity_ranges + 1))
      {   
      	
      	if (any(!is.na(hidden_PA_levels)))
      	{
      		    barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[i]], names.arg = names, main = "bouts in relation to their length (HMM)")
      	} else {
      			barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[i]], names.arg = names, main = "bouts in relation to their length")
      	}
      

        if (all(!is.na(names_activity_ranges)))
        {
          legend("topright", bg = "white", c("all", names_activity_ranges)[i], fill = c("gray"))
        } else {
          legend("topright", bg = "white", c("all",seq(1, number_of_activity_ranges))[i], fill = c("gray"))	
        }
      }
      par(mfrow = c(1,1))  
    }
    
    
    ############### summary of all results ####################################################################################################################################  	
    if (plotting == 1)
    {   
      par(mfrow = c(2,2))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
         plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i]  , '(', round( 100* rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '(', round( 100* rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      } else {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }		
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		plot(x, main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }		
      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      
            
      if (any(!is.na(hidden_PA_levels)))
      {
 			barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])),
              abs_freq_bouts_el[[1]], names.arg = names, main = "bouts in relation to their length (HMM)")
      } else {
 			barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])),
              abs_freq_bouts_el[[1]], names.arg = names, main = "bouts in relation to their length")
      }		    
      
      par(mfrow = c(1,1))    
    }	
    
    
    ############### time series of activity counts, classified into activity ranges #########################################################################################
    if (plotting == 2)
    {   
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	      plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
      	plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	      
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
    }	
    
    
    ############### time series of bouts (and, if available, the sequence of the estimated hidden physical activity levels, extracted by decoding a trained HMM) #############
    if (plotting == 3)
    {   
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }
      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
    }
    
    
    ############### barplots of absolute and relative frequencies of time spent in different activity ranges #################################################################
    if (plotting == 4)
    {   
      par(mfrow = c(2,1))
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i]  , '(', round( 100* rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '(', round( 100* rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    	barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      } else {
      			barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i]  ,'( # ',round(freq_acitvity_range[i], digits = 2),')')  )  )
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  ,'( # ',round( freq_acitvity_range[i], digits = 2),')')  )  )
        }
      }		
      

      
      if (any(!is.na(hidden_PA_levels)))
      {
 			barplot(freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "absolute frequency", names.arg = foo)		
      } else {
 			barplot(freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "absolute frequency", names.arg = foo)		
      }
      
      par(mfrow = c(1,1))  	
    }	
    
    
    ############### barplots of absolute frequencies of different bout intervals (overall and by activity ranges ) ##########################################################    
    if (plotting == 5)
    {   
      par(mfrow = c(number_of_activity_ranges/2+1,2))
      
      for (i in 1:(number_of_activity_ranges+1))
      {
      	
      if (any(!is.na(hidden_PA_levels)))
      {
      	    barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[i]], names.arg = names, main = "bouts in relation to their length (HMM)")
      } else {
      		barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[i]], names.arg = names, main = "bouts in relation to their length")
      }
      
      
        if (all(!is.na(names_activity_ranges)))
        {
          legend("topright", bg = "white", c("all",names_activity_ranges)[i], fill =c("gray"))
        } else {
          legend("topright", bg = "white", c("all",seq(1, number_of_activity_ranges))[i], fill =c("gray"))	
        }
      }
      par(mfrow = c(1,1))  
    }
  }
  
############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################  

print(list(activity_ranges = activity_ranges,
           classification = classification, 
           rel_freq_acitvity_range = rel_freq_acitvity_range,
           quantity_of_bouts = quantity_of_bouts, 
           abs_freq_bouts_el = abs_freq_bouts_el$all_activity_ranges))
  
  
return(list(activity_ranges = activity_ranges,
            classification = classification, 
            classification_per_activity_range = classification_per_activity_range,
            freq_acitvity_range = freq_acitvity_range, 
            rel_freq_acitvity_range = rel_freq_acitvity_range,
            quantity_of_bouts = quantity_of_bouts, 
            bout_periods = bout_periods,
            abs_freq_bouts_el = abs_freq_bouts_el))
}

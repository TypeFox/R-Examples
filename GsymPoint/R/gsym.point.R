gsym.point <-
function(methods, data, marker, status, tag.healthy, categorical.cov = NULL, CFN = 1, CFP = 1, control = control.gsym.point(), confidence.level = 0.95, trace = FALSE, seed = FALSE, value.seed = 3)
{

  # Possible errors:
  if(missing(methods) || is.null(methods))
  {
      stop("'methods' argument required.", call.=FALSE)
  }

  if(any(!(methods %in% c("GPQ","EL"))))
  {
      stop ("You have entered an invalid method.", call. = FALSE)
  }

  if (missing(data)|| is.null(data))
  {
      stop("'data' argument required.", call. = FALSE)
  }

  if (missing(marker)|| is.null(marker))
  {
      stop("'marker' argument required.", call. = FALSE)
  }

  if (missing(status)|| is.null(status))
  {
      stop("'status' argument required.", call. = FALSE)
  }

  if (missing(tag.healthy)|| is.null(tag.healthy))
  {
      stop("'tag.healthy' argument required.", call. = FALSE)
  }

  if (control$c_sampling <= 0)
  {
      stop("the constant for resampling in EL method must be a positive value.", call. = FALSE)
  }

  if (control$c_F <= 0)
  {
      stop("the constant for estimating the distribution in EL method must be a positive value.", call. = FALSE)
  }

  if (control$c_ELq <= 0)
  {
      stop("the constant for estimating the empirical likelihood function in EL method must be a positive value.", call. = FALSE)
  }

  if (control$c_R <= 0)
  {
      stop("the constant for estimating the ROC curve in EL method must be a positive value.", call. = FALSE)
  }

  if (control$B <= 0)
  {
      stop("the number of samples in EL method must be a positive value.", call. = FALSE)
  }

  if (control$I <= 0)
  {
      stop("the number of samples in GPQ method must be a positive value", call. = FALSE)	    
  }

  if (confidence.level < 0 | confidence.level > 1 | length(confidence.level) != 1)
  {
      stop("'confidence.level' must be a single number between 0 and 1.", call. = FALSE)
  }	   

  if (is.logical(trace) == FALSE)
  {
      stop("'trace' must be a logical-type argument.", call. = FALSE)   
  }

  if (is.logical(seed) == FALSE)
  {
      stop("'seed' must be a logical-type argument.", call. = FALSE)   
  }

  if (length(value.seed) != 1)
  {
      stop("'value.seed' must be a single number.", call. = FALSE)
  }	

  if(!all(c(marker,status,categorical.cov) %in% names(data))) {
      stop("Not all needed variables are supplied in 'data'.", call. = FALSE)
  }

  # NA's deleted
  data <- na.omit(data[, c(marker,status,categorical.cov)])

  # A data frame with the results is created:
  res <- vector("list", length(methods))
  names(res) <- methods

  if(!is.null(categorical.cov)) {
    	if(!is.factor(data[, categorical.cov])) data[, categorical.cov] <- factor(data[, categorical.cov])
        	data[, categorical.cov] <- droplevels(data[, categorical.cov])
        	levels.cat <- levels(data[, categorical.cov])
        	for (i in 1: length(methods)) {
            	res[[i]] <- vector("list", length(levels.cat))
            	names(res[[i]]) <- levels.cat
        }
  }

  else {
    		levels.cat = 1
    		res[[1]] <- vector("list", 1)
    		names(res[[1]]) <- "Global"
 }

  # Calculate the needed parameters:
  n0 <- length(data[data[,status] == tag.healthy, marker])
  n1 <- length(data[data[,status] != tag.healthy, marker])

  # If you we want to fix the seed for generating the trials:
  if(seed == TRUE) 
  {
	set.seed (value.seed)
  }

  # Each method is called up:
 	for (i in 1: length(methods))
  {
      	for(j in 1:length(levels.cat))
        {
        		if(trace)
            		{
        			cat("*************************************************\n")
        			text <- paste("Method: ", methods[i], sep = "")

        			if(length(levels.cat) > 1) {
        		 		 text <- paste(text, ". Level: ", levels.cat[j], sep = "")
        			}

        			cat(text)
        			cat("\nAnalysing ...")
        			cat("\n*************************************************\n")
        		}

         	  data.m <- if(length(levels.cat) != 1) data[data[,categorical.cov] == levels.cat[j], ] else data           	
                  rho = CFP/CFN
        	  res[[i]][[j]] <- eval(parse(text = paste("function.", methods[i], sep = "")))(data = data.m, marker = marker, status = status, tag.healthy = tag.healthy, CFN = CFN, CFP = CFP, control = control, confidence.level = confidence.level)
	          	  
        }
   }

  res$methods <-  methods

  if(length(methods) > 1)  # Both methods are considered
  {
	if (is.null(categorical.cov))
        {
		if (methods[1] == "GPQ")
		{
                	if(!"normality.transformed" %in% names(res$GPQ$Global))
  			{
				cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
				cat(paste("considered normally distributed in both groups.\n"))
				cat(paste("Therefore, although the results of both methods will be shown,\n"))
				cat(paste("the GPQ method would be more suitable for this dataset.\n"))
				
			      cat("\n")

				cat("Shapiro-Wilk test p-values\n")
                        cat("\n")

				print(matrix(c(round(res$GPQ$Global$pvalue.healthy,4), round(res$GPQ$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))

			}

			else
			{ 
				if (res$GPQ$Global$normality.transformed == "yes") 
				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can not \n"))
					cat(paste("be considered normally distributed in both groups.\n"))
					cat(paste("However, after transforming the marker using the Box-Cox \n"))
					cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
					cat(paste("indicates that the transformed marker can be considered \n"))
					cat(paste("normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))
						

					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$GPQ$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$GPQ$Global$pvalue.healthy,4), round(res$GPQ$Global$pvalue.diseased,4), round(res$GPQ$Global$pvalue.healthy.transformed,4), round(res$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				}
				
				else if (res$GPQ$Global$normality.transformed == "no")
				{
					cat(paste("According to the Shapiro-Wilk normality test, the original marker \n"))
					cat(paste("can not be considered normally distributed in both groups.\n"))
					cat(paste("After transforming the marker using the Box-Cox transformation\n"))
					cat(paste("estimate the Shapiro-Wilk normality test indicates that the \n"))
					cat(paste("transformed marker can not be considered normally distributed \n"))
					cat(paste("in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
					cat(paste("You must use the EL method instead.\n"))
				
					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$GPQ$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$GPQ$Global$pvalue.healthy,4), round(res$GPQ$Global$pvalue.diseased,4), round(res$GPQ$Global$pvalue.healthy.transformed,4), round(res$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				}

			}
		}	
        
		if (methods[1] == "EL")
		{
                	if(!"normality.transformed" %in% names(res$EL$Global))
  			{
				cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
				cat(paste("considered normally distributed in both groups.\n"))
				cat(paste("Therefore, although the results of both methods will be shown,\n"))
				cat(paste("the GPQ method would be more suitable for this dataset.\n"))

			      cat("\n")

				cat("Shapiro-Wilk test p-values\n")
                        cat("\n")

				print(matrix(c(round(res$EL$Global$pvalue.healthy,4), round(res$EL$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
			}

			else
			{ 
				if (res$EL$Global$normality.transformed == "yes") 
				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can not \n"))
					cat(paste("be considered normally distributed in both groups.\n"))
					cat(paste("However, after transforming the marker using the Box-Cox \n"))
					cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
					cat(paste("indicates that the transformed marker can be considered \n"))
					cat(paste("normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))
						
					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$EL$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$EL$Global$pvalue.healthy,4), round(res$EL$Global$pvalue.diseased,4), round(res$EL$Global$pvalue.healthy.transformed,4), round(res$EL$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				}


				else if (res$EL$Global$normality.transformed == "no")
				{
					cat(paste("According to the Shapiro-Wilk normality test, the original marker \n"))
					cat(paste("can not be considered normally distributed in both groups.\n")) 
					cat(paste("After transforming the marker using the Box-Cox transformation \n"))
					cat(paste("estimate the Shapiro-Wilk normality test indicates that the \n"))
					cat(paste("transformed marker can not be considered normally distributed \n"))
					cat(paste("in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
					cat(paste("You must use the EL method instead.\n"))
				
					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$EL$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$EL$Global$pvalue.healthy,4), round(res$EL$Global$pvalue.diseased,4), round(res$EL$Global$pvalue.healthy.transformed,4), round(res$EL$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))	
				}

			}
		}	
	}
	else
        {
		for(i in 1:length(levels.cat))
		{
			if(!"normality.transformed" %in% names(res$EL[[i]]))
			{
				cat(paste(levels.cat[i],": \n"))
				cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
				cat(paste("considered normally distributed in both groups.\n"))
				cat(paste("Therefore, although the results of both methods will be shown,\n"))
				cat(paste("the GPQ method would be more suitable for this dataset.\n"))

			      cat("\n")

				cat("Shapiro-Wilk test p-values\n")
                        cat("\n")

				print(matrix(c(round(res$EL[[i]][["pvalue.healthy"]],4), round(res$EL[[i]][["pvalue.diseased"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
				cat("\n")
			}
			
			else
			{
				if (res$EL[[i]][["normality.transformed"]] == "yes")
				{
					cat(paste(levels.cat[i],":\n"))
					cat(paste("According to the Shapiro-Wilk normality test, the marker can not \n"))
					cat(paste("be considered normally distributed in both groups.\n"))
					cat(paste("However, after transforming the marker using the Box-Cox \n"))
					cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
					cat(paste("indicates that the transformed marker can be considered \n"))
					cat(paste("normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))
					
					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$EL[[i]][["lambda"]],4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$EL[[i]][["pvalue.healthy"]],4), round(res$EL[[i]][["pvalue.diseased"]],4), round(res$EL[[i]][["pvalue.healthy.transformed"]],4), round(res$EL[[i]][["pvalue.diseased.transformed"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					cat("\n")
				}

				else if (res$EL[[i]][["normality.transformed"]] == "no")	            
				{
					cat(paste(levels.cat[i],":\n"))
					cat(paste("According to the Shapiro-Wilk normality test, the original marker \n"))
					cat(paste("can not be considered normally distributed in both groups.\n")) 
					cat(paste("After transforming the marker using the Box-Cox transformation \n"))
					cat(paste("estimate the Shapiro-Wilk normality test indicates that the \n"))
					cat(paste("transformed marker can not be considered normally distributed \n"))
					cat(paste("in both groups.\n")) 
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
					cat(paste("You must use the EL method instead.\n"))
		
					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$EL[[i]][["lambda"]],4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$EL[[i]][["pvalue.healthy"]],4), round(res$EL[[i]][["pvalue.diseased"]],4), round(res$EL[[i]][["pvalue.healthy.transformed"]],4), round(res$EL[[i]][["pvalue.diseased.transformed"]],4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					cat("\n")
				}
		 	}
          	}
    	}

    }

    else
    {
	if (methods == "EL")
	{
		if (is.null(categorical.cov))
            {
			if(!"normality.transformed" %in% names(res$EL$Global))
			{
				cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
				cat(paste("considered normally distributed in both groups.\n"))
				cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))
				cat("\n")


				cat("Shapiro-Wilk test p-values\n")
                        cat("\n")

				print(matrix(c(round(res$EL$Global$pvalue.healthy,4), round(res$EL$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
			}
			else 
			{
				if (res$EL$Global$normality.transformed == "yes")		
				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can not \n"))
					cat(paste("be considered normally distributed in both groups.\n")) 
					cat(paste("However, after transforming the marker using the Box-Cox \n"))
					cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
					cat(paste("indicates that the transformed marker can be considered \n"))
					cat(paste("normally distributed in both groups.\n"))
					cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

					cat("\n")

					cat(paste("Box-Cox lambda estimate =", round(res$EL$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$EL$Global$pvalue.healthy,4), round(res$EL$Global$pvalue.diseased,4), round(res$EL$Global$pvalue.healthy.transformed,4), round(res$EL$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				}
			}
		}
		else
		{
			for(i in 1:length(levels.cat))
			{
				if(!"normality.transformed" %in% names(res[[1]][[i]]))
				{
					cat(paste(levels.cat[i],":\n"))
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
					cat(paste("considered normally distributed in both groups. \n"))
					cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(res[[1]][[i]][["pvalue.healthy"]],4), round(res[[1]][[i]][["pvalue.diseased"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
					cat("\n")
				}

				else 
				{
					if (res[[1]][[i]][["normality.transformed"]] == "yes") 
					{
						cat(paste(levels.cat[i],":\n"))
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not \n"))
						cat(paste("be considered normally distributed in both groups.\n")) 
						cat(paste("However, after transforming the marker using the Box-Cox \n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
						cat(paste("indicates that the transformed marker can be considered \n"))
						cat(paste("normally distributed in both groups.\n")) 
						cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(res[[1]][[i]][["lambda"]],4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(res[[1]][[i]][["pvalue.healthy"]],4), round(res[[1]][[i]][["pvalue.diseased"]],4), round(res[[1]][[i]][["pvalue.healthy.transformed"]],4), round(res[[1]][[i]][["pvalue.diseased.transformed"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
						cat("\n")
					}
				}
			}


		}

	}

	if (methods == "GPQ")
	{
		if (is.null(categorical.cov))
            	{
			if ("normality.transformed" %in% names(res$GPQ$Global))
			{
				if(res$GPQ$Global$normality.transformed == "no")
				{
					cat(paste("According to the Shapiro-Wilk normality test, the original marker \n"))
					cat(paste("can not be considered normally distributed in both groups.\n")) 
					cat(paste("After transforming the marker using the Box-Cox transformation \n"))
					cat(paste("estimate the Shapiro-Wilk normality test indicates that the \n"))
					cat(paste("transformed marker can not be considered normally distributed \n"))
					cat(paste("in both groups. \n"))
					cat(paste("Therefore, the results obtained with the GPQ method may not be \n"))
					cat(paste("reliable. You must use the EL method instead.\n"))
					cat("\n")


					cat(paste("Box-Cox lambda estimate =", round(res$GPQ$Global$lambda,4),"\n"))
					cat("\n")

					cat("Shapiro-Wilk test p-values\n")

					print(matrix(c(round(res$GPQ$Global$pvalue.healthy,4), round(res$GPQ$Global$pvalue.diseased,4), round(res$GPQ$Global$pvalue.healthy.transformed,4), round(res$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				}

			}
		}

		else
		{
			for(i in 1:length(levels.cat))
			{
                              if ("normality.transformed" %in% names(res[[1]][[i]]))
      				{

						if (res[[1]][[i]][["normality.transformed"]] == "no")
						{
							cat(paste(levels.cat[i],":\n"))
							cat(paste("According to the Shapiro-Wilk normality test, the original marker \n"))
							cat(paste("can not be considered normally distributed in both groups.\n"))
							cat(paste("After transforming the marker using the Box-Cox transformation \n"))
							cat(paste("estimate the Shapiro-Wilk normality test indicates that the \n"))
							cat(paste("transformed marker can not be considered normally distributed \n"))
							cat(paste("in both groups.\n"))
							cat(paste("Therefore, the results obtained with the GPQ method may not be \n"))
							cat(paste("reliable. You must use the EL method instead.\n"))
							cat("\n")

							cat(paste("Box-Cox lambda estimate =", round(res[[1]][[i]][["lambda"]],4),"\n"))
							cat("\n")

							cat("Shapiro-Wilk test p-values\n")

							print(matrix(c(round(res[[1]][[i]][["pvalue.healthy"]],4), round(res[[1]][[i]][["pvalue.diseased"]],4), round(res[[1]][[i]][["pvalue.healthy.transformed"]],4), round(res[[1]][[i]][["pvalue.diseased.transformed"]],4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
							cat("\n")
						}

					}
			}

		}

	}

  }


  if(length(levels.cat) != 1) res$levels.cat  <-  levels.cat
  res$call <- match.call()
  res$data <- data

  class(res) <- "gsym.point"

  invisible(res)
  res
}

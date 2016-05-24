print.summary.gsym.point <-
function(x, ...) {
           
	cat("\n*************************************************\n")
	cat("OPTIMAL CUTOFF: GENERALIZED SYMMETRY POINT")
	cat("\n*************************************************\n")
      cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
      cat("\n")

	if (is.null(x$levels.cat))  # no categorical covariate
      {
		if (length(x$methods)>1)  # Both methods are considered
		{
			if (x$methods[1] == "GPQ")
			{
                		if(!"normality.transformed" %in% names(x$GPQ$Global))
  				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
					cat(paste("considered normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))

			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4), round(x$GPQ$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))

				}

				else
				{ 
					if (x$GPQ$Global$normality.transformed == "yes") 
					{					
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
						cat(paste("be considered normally distributed in both groups.\n"))
						cat(paste("However, after transforming the marker using the Box-Cox\n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
						cat(paste("indicates that the transformed marker can be considered\n"))
						cat(paste("normally distributed in both groups.\n"))
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the GPQ method would be more suitable for this dataset.\n"))
						
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$GPQ$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4), round(x$GPQ$Global$pvalue.diseased,4), round(x$GPQ$Global$pvalue.healthy.transformed,4), round(x$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))

					}
					else if (x$GPQ$Global$normality.transformed == "no")
					{
						cat(paste("According to the Shapiro-Wilk normality test, the original marker\n"))
						cat(paste("can not be considered normally distributed in both groups.\n")) 
						cat(paste("After transforming the marker using the Box-Cox transformation\n"))
						cat(paste("estimate the Shapiro-Wilk normality test indicates that the\n"))
						cat(paste("transformed marker can not be considered normally distributed\n"))
						cat(paste("in both groups.\n")) 
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
						cat(paste("You must use the EL method instead.\n"))
				

						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$GPQ$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4), round(x$GPQ$Global$pvalue.diseased,4), round(x$GPQ$Global$pvalue.healthy.transformed,4), round(x$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					
					}

				}
			}

			if (x$methods[1] == "EL")
			{
                		if(!"normality.transformed" %in% names(x$EL$Global))
  				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be\n"))
					cat(paste("considered normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))

			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(x$EL$Global$pvalue.healthy,4), round(x$EL$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
				}
				else
				{ 
					if (x$EL$Global$normality.transformed == "yes") 
					{
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
						cat(paste("be considered normally distributed in both groups.\n")) 
						cat(paste("However, after transforming the marker using the Box-Cox\n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
						cat(paste("indicates that the transformed marker can be considered\n"))
						cat(paste("normally distributed in both groups.\n"))
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the GPQ method would be more suitable for this dataset.\n"))
						
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$EL$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$EL$Global$pvalue.healthy,4), round(x$EL$Global$pvalue.diseased,4), round(x$EL$Global$pvalue.healthy.transformed,4), round(x$EL$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					}
					
					else if (x$EL$Global$normality.transformed == "no")
					{
						cat(paste("According to the Shapiro-Wilk normality test, the original marker\n"))
						cat(paste("can not be considered normally distributed in both groups.\n")) 
						cat(paste("After transforming the marker using the Box-Cox transformation\n"))
						cat(paste("estimate the Shapiro-Wilk normality test indicates that the\n"))
						cat(paste("transformed marker can not be considered normally distributed\n"))
						cat(paste("in both groups.\n")) 
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
						cat(paste("You must use the EL method instead.\n"))
				
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$EL$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$EL$Global$pvalue.healthy,4), round(x$EL$Global$pvalue.diseased,4), round(x$EL$Global$pvalue.healthy.transformed,4), round(x$EL$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))	
					}

				}
			}	
		}
		else
		{
			if (x$methods == "GPQ")
			{
				if(!"normality.transformed" %in% names(x$GPQ$Global))
      			{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be\n"))
					cat(paste("considered normally distributed in both groups.\n"))
			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4),round(x$GPQ$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))

				}
		
				else  
				{
					if (x$GPQ$Global$normality.transformed == "yes")
					{
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
						cat(paste("be considered normally distributed in both groups. \n"))
						cat(paste("However, after transforming the marker using the Box-Cox \n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test \n"))
						cat(paste("indicates that the transformed marker can be considered \n"))
						cat(paste("normally distributed in both groups.\n"))
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$GPQ$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4), round(x$GPQ$Global$pvalue.diseased,4), round(x$GPQ$Global$pvalue.healthy.transformed,4), round(x$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))

					}

					else
					{
						cat(paste("According to the Shapiro-Wilk normality test, the original marker\n"))
						cat(paste("can not be considered normally distributed in both groups.\n")) 
						cat(paste("After transforming the marker using the Box-Cox transformation\n"))
						cat(paste("estimate the Shapiro-Wilk normality test indicates that the\n"))
						cat(paste("transformed marker can not be considered normally distributed\n"))
						cat(paste("in both groups.\n"))
						cat(paste("Therefore, the results obtained with the GPQ method may not be\n"))
						cat(paste("reliable. You must use the EL method instead.\n"))
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$GPQ$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$GPQ$Global$pvalue.healthy,4), round(x$GPQ$Global$pvalue.diseased,4), round(x$GPQ$Global$pvalue.healthy.transformed,4), round(x$GPQ$Global$pvalue.diseased.transformed,4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					}

				}
			}

      		if(x$methods == "EL") 
			{ 
				if(!"normality.transformed" %in% names(x$EL$Global))
      			{      
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be\n"))
					cat(paste("considered normally distributed in both groups.\n"))
					cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(x$EL$Global$pvalue.healthy,4), round(x$EL$Global$pvalue.diseased,4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))

      			}

				else
				{
					if (x$EL$Global$normality.transformed == "yes")
					{
      
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
						cat(paste("be considered normally distributed in both groups.\n"))
						cat(paste("However, after transforming the marker using the Box-Cox\n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
						cat(paste("indicates that the transformed marker can be considered\n"))
						cat(paste("normally distributed in both groups.\n"))
						cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$EL$Global$lambda,4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$EL$Global$pvalue.healthy,4), round(x$EL$Global$pvalue.diseased,4), round(x$EL$Global$pvalue.healthy.transformed,4), round(x$EL$Global$pvalue.diseased.transformed,4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))

					}
				}
      		}
		}
	}

	for(i in 1: length(x$p.table)) {
		if(!is.null(x$levels.cat)) {    # Categorical covariate is considered
			cat("\n*************************************************\n")
			cat(names(x$p.table)[i])
			cat("\n*************************************************\n")

			if (length(x$methods)>1)  # Both methods are considered
			{
				if(!"normality.transformed" %in% names(x$EL[[i]]))
				{
					cat(paste("According to the Shapiro-Wilk normality test, the marker can be\n"))
					cat(paste("considered normally distributed in both groups.\n"))
					cat(paste("Therefore, although the results of both methods will be shown,\n"))
					cat(paste("the GPQ method would be more suitable for this dataset.\n"))

			            cat("\n")

					cat("Shapiro-Wilk test p-values\n")
                              cat("\n")

					print(matrix(c(round(x$EL[[i]][["pvalue.healthy"]],4), round(x$EL[[i]][["pvalue.diseased"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))

				}
			
				else
				{
					if (x$EL[[i]][["normality.transformed"]] == "yes")
					{
						cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
						cat(paste("be considered normally distributed in both groups.\n"))
						cat(paste("However, after transforming the marker using the Box-Cox\n"))
						cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
						cat(paste("indicates that the transformed marker can be considered\n"))
						cat(paste("normally distributed in both groups.\n"))
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the GPQ method would be more suitable for this dataset.\n"))
						
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$EL[[i]][["lambda"]],4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$EL[[i]][["pvalue.healthy"]],4), round(x$EL[[i]][["pvalue.diseased"]],4), round(x$EL[[i]][["pvalue.healthy.transformed"]],4), round(x$EL[[i]][["pvalue.diseased.transformed"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))

					}
				
					else if (x$EL[[i]][["normality.transformed"]] == "no")	            
					{
						cat(paste("According to the Shapiro-Wilk normality test, the original marker\n"))
						cat(paste("can not be considered normally distributed in both groups.\n")) 
						cat(paste("After transforming the marker using the Box-Cox transformation\n"))
						cat(paste("estimate the Shapiro-Wilk normality test indicates that the\n"))
						cat(paste("transformed marker can not be considered normally distributed\n"))
						cat(paste("in both groups.\n"))
						cat(paste("Therefore, although the results of both methods will be shown,\n"))
						cat(paste("the results obtained with the GPQ method may not be reliable.\n"))
						cat(paste("You must use the EL method instead.\n"))
				
						cat("\n")

						cat(paste("Box-Cox lambda estimate =", round(x$EL[[i]][["lambda"]],4),"\n"))
						cat("\n")

						cat("Shapiro-Wilk test p-values\n")

						print(matrix(c(round(x$EL[[i]][["pvalue.healthy"]],4), round(x$EL[[i]][["pvalue.diseased"]],4), round(x$EL[[i]][["pvalue.healthy.transformed"]],4), round(x$EL[[i]][["pvalue.diseased.transformed"]],4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))

					}

		 		}
          		}
    	
			else
			{                  
                  	if(x$methods == "EL") 
				{ 
					if(!"normality.transformed" %in% names(x$EL[[i]]))
      				{      
						cat(paste("According to the Shapiro-Wilk normality test, the marker can be\n"))
						cat(paste("considered normally distributed in both groups.\n"))
						cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

			            	cat("\n")

						cat("Shapiro-Wilk test p-values\n")
                              	cat("\n")

						print(matrix(c(round(x$EL[[i]][["pvalue.healthy"]],4), round(x$EL[[i]][["pvalue.diseased"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
			
      				}

					else
					{
						if (x$EL[[i]][["normality.transformed"]] == "yes")
						{      
							cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
							cat(paste("be considered normally distributed in both groups.\n"))
							cat(paste("However, after transforming the marker using the Box-Cox\n"))
							cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
							cat(paste("indicates that the transformed marker can be considered\n"))
							cat(paste("normally distributed in both groups.\n"))
							cat(paste("Therefore the GPQ method would be more suitable for this dataset.\n"))

							cat("\n")

							cat(paste("Box-Cox lambda estimate =", round(x$EL[[i]][["lambda"]],4),"\n"))
							cat("\n")

							cat("Shapiro-Wilk test p-values\n")

							print(matrix(c(round(x$EL[[i]][["pvalue.healthy"]],4), round(x$EL[[i]][["pvalue.diseased"]],4), round(x$EL[[i]][["pvalue.healthy.transformed"]],4), round(x$EL[[i]][["pvalue.diseased.transformed"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
				
						}
					}

				}

				if (x$methods == "GPQ")
				{
					if(!"normality.transformed" %in% names(x$GPQ[[i]]))
      				{
						cat(paste("According to the Shapiro-Wilk normality test, the marker can be \n"))
						cat(paste("considered normally distributed in both groups.\n"))
			            	cat("\n")

						cat("Shapiro-Wilk test p-values\n")
                              	cat("\n")

						print(matrix(c(round(x$GPQ[[i]][["pvalue.healthy"]],4), round(x$GPQ[[i]][["pvalue.diseased"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker"),c("Group 0","Group 1"))))
					}

					else  
					{
						if (x$GPQ[[i]][["normality.transformed"]] == "yes")
						{
							cat(paste("According to the Shapiro-Wilk normality test, the marker can not\n"))
							cat(paste("be considered normally distributed in both groups.\n"))
							cat(paste("However, after transforming the marker using the Box-Cox\n"))
							cat(paste("transformation estimate, the Shapiro-Wilk normality test\n"))
							cat(paste("indicates that the transformed marker can be considered\n"))
							cat(paste("normally distributed in both groups.\n"))
							cat("\n")

							cat(paste("Box-Cox lambda estimate =", round(x$GPQ[[i]][["lambda"]],4),"\n"))
							cat("\n")

							cat("Shapiro-Wilk test p-values\n")

							print(matrix(c(round(x$GPQ[[i]][["pvalue.healthy"]],4), round(x$GPQ[[i]][["pvalue.diseased"]],4), round(x$GPQ[[i]][["pvalue.healthy.transformed"]],4), round(x$GPQ[[i]][["pvalue.diseased.transformed"]],4)), ncol=2, byrow=T, dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
						}

						else
						{
							cat(paste("According to the Shapiro-Wilk normality test, the original marker\n"))
							cat(paste("can not be considered normally distributed in both groups.\n"))
							cat(paste("After transforming the marker using the Box-Cox transformation\n"))
							cat(paste("estimate the Shapiro-Wilk normality test indicates that the\n"))
							cat(paste("transformed marker can not be considered normally distributed\n"))
							cat(paste("in both groups. \n"))
							cat(paste("Therefore, the results obtained with the GPQ method may not be\n"))
							cat(paste("reliable. You must use the EL method instead.\n"))
							cat("\n")

							cat(paste("Box-Cox lambda estimate =", round(x$GPQ[[i]][["lambda"]],4),"\n"))
							cat("\n")

							cat("Shapiro-Wilk test p-values\n")

							print(matrix(c(round(x$GPQ[[i]][["pvalue.healthy"]],4), round(x$GPQ[[i]][["pvalue.diseased"]],4), round(x$GPQ[[i]][["pvalue.healthy.transformed"]],4), round(x$GPQ[[i]][["pvalue.diseased.transformed"]],4)), ncol=2,dimnames=list(c("Original marker", "Box-Cox transformed marker"),c("Group 0","Group 1"))))
					
						}

					}
				}

			}
		}
          
		cat("\nArea under the ROC curve (AUC): ", x$p.table[[i]][["AUC"]], "\n\n")
            	if (length(x$p.table[[i]]) == 1) j = 1
		if (length(x$p.table[[i]]) > 1) {
			for (j in 1:(length(x$p.table[[i]]) - 1)) {
				cat(paste("METHOD: " , names(x$p.table[[i]])[j], sep = ""))
				cat("\n\n")
				if(length(x$p.table[[i]][[j]]) != 0) {
					for (k in 1:length(x$p.table[[i]][[j]])) {
						print(x$p.table[[i]][[j]][[k]], quote = FALSE, justify = "right", na.print = "-")
						cat("\n")
					}
				}
			}
		}
	}
	cat("\n\n")
   	invisible(x)
}

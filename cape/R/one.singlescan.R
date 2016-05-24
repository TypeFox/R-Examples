one.singlescan <-
function(phenotype.vector, genotype.mat, covar.vector = NULL, n.cores = 2) {

		g = NULL #for appeasing R CMD check

		#This function gets regression statistics with a
		#covariate table
		get.stats <- function(phenotype, genotype, covar){

				
				model <- lm(phenotype~cbind(genotype,covar))
				
		
				#take the last line of coefficients.
				model.coef <- summary(model)$coefficients
				slope <- model.coef[dim(model.coef)[1],1]
				se <- model.coef[dim(model.coef)[1],2]
				t.stat <- abs(model.coef[dim(model.coef)[1],3])
				p.val <- model.coef[dim(model.coef)[1],4]
							
			#put together all the statistics we want to keep
			#we keep the absolute value of the t statistic,
			#the p value, and the covariate flag

			table.row <- c(slope, se, t.stat, p.val)
			return(table.row)
			}
	
	
	#==========================================

		
		#take out the response variable
				
		#apply the modeling function to each marker column
		cl <- makeCluster(n.cores)
		registerDoParallel(cl)
		results.table <- foreach(g = genotype.mat, .combine = "rbind") %dopar% {
			get.stats(phenotype = phenotype.vector, genotype = g, covar = covar.vector)
			}
		stopCluster(cl)
		
		colnames(results.table) <- c("slope", "se", "t.stat", "p.val")
		rownames(results.table) <- colnames(genotype.mat)
		return(results.table)
	
	}

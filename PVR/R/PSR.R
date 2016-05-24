PSR <- function(x, trait = NULL, null.model = FALSE, Brownian.model = FALSE, times = 1000){
####---------Computing psr curve
		
		pvr <- x@Eigen
		pD <- x@phyDist
		psr <- data.frame(r.squared = 0, eigenvalues = 0)

		#---------loading species traits
		if(is.character(trait)){
			
			trait <- read.table(trait)
		}
		
		trait <- as.matrix(trait)
		nsp <- nrow(trait)
		Ntraits <- ncol(trait)
		areas <- data.frame(PSR.area = 0, null.p = NA, Brownian.p = NA, iterations = "none")
		Naxis <- ncol(pvr$vectors)
		binary <- FALSE
			
		if(max(trait[,1]) == 1 & min(trait[,1]) == 0 & length(unique(trait[ ,1])) == 2){
			
				warning("Binary traits: using Brownian expectancy as null model", "\n")
				null.model = FALSE
				binary <- TRUE
		}

		
		for(k in 1:Ntraits){
			
			pvr$vectors[ ,1:Naxis] <- t(t(pvr$vectors[ ,1:Naxis])/sqrt(pvr$values[1:Naxis]))
			relVal <- pvr$values/sum(pvr$values)
			
			for(i in 1:Naxis){
				
				psr[i, 1] <- summary(lm(trait[,k] ~ matrix(as.double(pvr$vectors[,1:i]), ncol = i)))$r.squared
				psr[i, 2] <- sum(relVal[1:i])	
		
			}
			
			#---------Computing psrarea
			coords = data.frame(x = c(0, psr[ ,2]), y = c(0, psr[ ,1]))
			while(round(coords[nrow(coords), 1],7) == round(1,7) & round(coords[nrow(coords), 2],7) == round(1,7)){
				
				coords <- coords[-nrow(coords),]
			}
			
			coords[(nrow(coords) + 1), ] <- 1
			psrarea <- .Area(coords)
			areas[1, 1] <- psrarea
			
			#---------Creating PSRarea null distribution
			if(null.model){
				
				nullDistribution <- matrix(0, nrow = Naxis, ncol = times)
				psrNull <- data.frame(r.squared = 0, eigenvalues = 0)
				nullPsrarea <- data.frame(a = 0)
				for(t in 1:times){
					
					#---------shuffling trait vector
					traitRand <- sample(trait[,k])
					
					#---------Computing psr random curve
					for(i in 1:Naxis){
						
						psrNull[i, 1] <- summary(lm(traitRand ~ matrix(as.double(pvr$vectors[,1:i]), ncol = i)))$r.squared
						psrNull[i, 2] <- sum(relVal[1:i])
					
					}
					
					nullDistribution[ ,t] <- psrNull[ ,1]
					#---------Computing psrarea
					coordsRand = data.frame(x = c(0, psrNull[ ,2]), y = c(0, psrNull[ ,1]))
					while(round(coordsRand[nrow(coordsRand), 1],7) == round(1,7) & round(coordsRand[nrow(coordsRand), 2],7) == round(1,7)){
						
						coordsRand <- coordsRand[-nrow(coordsRand),]
					}
					coordsRand[(nrow(coordsRand) + 1), ] <- 1
					
					nullPsrarea[t,1] <- .Area(as.matrix(coordsRand))
				}
				
				#---------Computing p
				distribCum <- ecdf(nullPsrarea[,1])
				nullp2 <- 1 - distribCum(psrarea)
				areas[1, 2] <- nullp2
				
				areas$iterations <- as.numeric(times)
			}
#			if(null.model){ used in version 0.1
#				
#				areas[k, 3] <- p2
#				nullareas[k, 1] <- paste("traits.set.", k, sep = "")
#				nullareas[k, 2] <- mean(nullPsrarea[ ,1])
#				nullareas[k, 3] <- var(nullPsrarea[ ,1])
#			}
			 
		}
			#---------Creating PSRarea neutral distribution
			if(Brownian.model){
				
				BrownianDistribution <- matrix(0, nrow = Naxis, ncol = times)
				psrBrownian <- data.frame(r.squared = 0, eigenvalues = 0)
				BrownianPsrarea <- data.frame(a = 0)
				phy <- x@phylo
				for(t in 1:times){
					
					#---------simulating trait vector

						traitRand <- rTraitCont(phy, mpdel = "BM")
						
					if(binary){
						
						traitRand <- rTraitCont(phy, mpdel = "BM")
						ones <- sum(trait[ ,k])
						zeros <- nsp - ones
						ind <- 1:nsp
						ind <- ind[order(traitRand)]
						traitRand <- c(rep(1, ones), rep(0, zeros))
						traitRand <- traitRand[order(ind)]
					}
						
					#---------Computing psr random curve
					for(i in 1:Naxis){
									
						psrBrownian[i, 1] <- summary(lm(traitRand ~ matrix(as.double(pvr$vectors[,1:i]), ncol = i)))$r.squared
						psrBrownian[i, 2] <- sum(relVal[1:i])
					}
					
					BrownianDistribution[ ,t] <- psrBrownian[ ,1]
					
					#---------Computing psrarea
					coordsBrownian = data.frame(x = c(0, psrBrownian[ ,2]), y = c(0, psrBrownian[ ,1]))
					while(round(coordsBrownian[nrow(coordsBrownian), 1],7) == round(1,7) & round(coordsBrownian[nrow(coordsBrownian), 2],7) == round(1,7)){
						
						coordsBrownian <- coordsBrownian[-nrow(coordsBrownian),]
					}
					
					coordsBrownian[(nrow(coordsBrownian) + 1), ] <- 1
					BrownianPsrarea[t,1] <- .Area(as.matrix(coordsBrownian))
				}
				
				#---------Computing p
				distribCum <- ecdf(BrownianPsrarea[,1])
				Brownianp2 <- 1 - distribCum(psrarea)
				areas[1, 3] <- Brownianp2
				
				areas$iterations <- as.numeric(times)
			}
		results <- new("PSR", x)
		results@PSRarea <- areas
		results@PSR <- data.frame(Cumul.eigen.values = psr[,2], r.squared = psr[,1])

		Expect.area.values <- list(
				
				meanNull = as.numeric(NA),
				varianceNull = as.numeric(NA),
				meanBrownian = as.numeric(NA),
				varianceBrownian = as.numeric(NA)
				)
		
		if(null.model){
			
			Expect.area.values$meanNull <- mean(nullPsrarea[ ,1]) 
			Expect.area.values$varianceNull <- var(nullPsrarea[ ,1])
			results@nullPSR <- nullDistribution
			
		}
		
		if(Brownian.model){
			
			Expect.area.values$meanBrownian <- mean(BrownianPsrarea[ ,1]) 
			Expect.area.values$varianceBrownian <- var(BrownianPsrarea[ ,1])
			results@BrownianPSR <- BrownianDistribution
		}
		
		results@Expect.area.values <- Expect.area.values
		if(binary){
			
			attr(results, "trait.type") <- "binary"
		} else
			attr(results, "trait.type") <- "continuous"

		return(results)
}
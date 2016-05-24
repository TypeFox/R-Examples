"mfit" <-

function(y, ..., dec = 3.)



{



	####  dec controls number of decimal places



	temp <- list(...)



	locall <- sys.call()



	nvar <- length(temp)



	omean <- mean(y)



	cat(" ", fill = TRUE)



	cat("For more decimals add ,dec=k after variables", fill = TRUE)



	cat(" ", fill = TRUE)



	cat("Overall Mean of Y variable", as.character(locall[2.]), " = ",

	

		round(omean, dec), fill = TRUE)



	cat(" ", fill = TRUE)



	cat("******************************", fill = TRUE)



	for(k in 1.:nvar) {



		cat("Fitted main Effect of Y variable", as.character(locall[2.]),



			"by X\nvariable", as.character(locall[k + 2.]), fill = TRUE)



		cat(" ", fill = TRUE)



		ord <- order(temp[[k]])



		temp.ord <- temp[[k]][ord]



		y.ord <- y[ord]



		tstat <- round(stats(y.ord, by = temp.ord)[1.:2.,  ], dec)



		tstat <- data.frame(N = tstat[1,  ], Main.Effect = 



			round(tstat[2,  ] - omean, dec))



		print(tstat)



		cat("******************************", fill = TRUE)



	}



	if(nvar > 2.)



		for(k in 1.:(nvar - 1.)) {



			for(j in (k + 1.):nvar) {



				cat("Table of 2-way Fitted Interaction Effects for",



					as.character(locall[2.]), fill = TRUE)



				cat("by X variables", as.character(locall[k + 2.]),



					"and", as.character(locall[j + 2.]), fill = TRUE)



				cat(" ", fill = TRUE)



				print(mfit.2way(y, temp[[k]], temp[[



					j]]), dec = dec)



				cat("******************************",



					fill = TRUE)



			}



		}



	else cat("Table of 2-way Fitted Interaction Effects for", 



			as.character(locall[2.]), fill = TRUE)



	cat("by X variables", as.character(locall[3.]), "and",

	

		as.character(locall[4.]), fill = TRUE)



	if(nvar == 2.)



		cat(" ", fill = TRUE)



	mfit.2way(y, temp[[1.]], temp[[2.]], dec = dec)



}


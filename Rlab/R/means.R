"means" <-

function(y, ..., dec = 3.)



{



	####  dec controls number of decimal places



	temp <- list(...)



	locall <- sys.call()



	nvar <- length(temp)



	cat(" ", fill = TRUE)



	for(k in 1.:nvar) {



		cat("Mean of Y variable", as.character(locall[2.]), "by X variable",



			as.character(locall[k + 2.]), fill = TRUE)



		cat(" ", fill = TRUE)



		ord <- order(temp[[k]])



		temp.ord <- temp[[k]][ord]



		y.ord <- y[ord]



		print(round(stats(y.ord, by = temp.ord)[



			1.:3.,  ], dec))



		cat("******************************", fill = TRUE)



	}



	cat("Note: in the following, if sample sizes are not equal\n")



	####  two-way means



	cat("row and col means will not be the average of cell means.\n\n"



		)



	if(nvar > 2.)



		for(k in 1.:(nvar - 1.)) {



			for(j in (k + 1.):nvar) {



				cat("Mean of Y Variable", as.character(locall[2.]),



					"by X variables", as.character(locall[k + 2.]),

					

					"and", as.character(locall[j + 2.]),



					fill = TRUE)



				cat(" ", fill = TRUE)



				print(means.2way(y, temp[[k]], temp[[



					j]]), dec = dec)



				cat("******************************",



					fill = TRUE)



			}



		}



	else cat("Mean of Y variable", as.character(locall[2.]), "by X variables",



			as.character(locall[3.]), "and", as.character(locall[4.]), fill = TRUE)



	if(nvar == 2.)



		means.2way(y, temp[[1.]], temp[[2.]], dec = dec)



}


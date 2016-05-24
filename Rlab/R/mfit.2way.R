"mfit.2way" <-

function(y, a, b, dec = 3.)



{



	omean <- mean(y)



	ord <- order(a, b)



	a <- a[ord]



	b <- b[ord]



	y <- y[ord]



	rlab <- unique(a)



	rlab <- as.character(rlab)



	clab <- unique(b)



	clab <- as.character(clab)



	temp <- stats(y, by = paste(a, b))



	temp2 <- temp[2.,  ]



	temp3 <- matrix(temp2, ncol = length(clab), nrow = length(



		rlab), byrow = TRUE)



	temp3 <- cbind(temp3, stats(y, by = a)[2.,  ])



	temp3 <- rbind(temp3, c(stats(y, by = b)[2.,



		], omean))



	for(j in 1:length(rlab)) {



		for(k in 1:length(clab)) {



			temp3[j, k] <- temp3[j, k] - temp3[j, (length(



				clab) + 1)] - temp3[(length(rlab) +



				1), k] + omean



		}



	}



	temp3 <- temp3[1:length(rlab), 1:length(clab)]



	dimnames(temp3) <- list(rlab, clab)



	round(temp3, dec)



}


"means.2way" <-

function(y, a, b, dec = 3.)



{



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



		], mean(y)))



	dimnames(temp3) <- list(c(rlab, "Col Mean"), c(clab, "Row Mean"



		))



	round(temp3, dec)



}


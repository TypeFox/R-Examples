"mplot" <-

function(y, ..., both = FALSE)



{



	temp <- list(...)



	locall <- sys.call()



	locall <- as.character(locall)



	nvar <- length(temp)



	if( (nvar == 2.) && !(both) )



		set.panel(1., 1.)



	else if( (nvar == 2.) && (both) )



		set.panel(2., 1.)



	else if( !(both) ) set.panel(choose(nvar,2), 1.)



	else if( (both) ) set.panel(choose(nvar,2), 2.)





	for(k in 1.:(nvar - 1.)) {



		for(j in (k + 1.):nvar) {



			interaction.plot(temp[[k]], temp[[j]],



				y, xlab = locall[k + 2.], 



				trace.label = locall[j + 2.],



				ylab = locall[2.])



			title("Mean of Y for values of X1 (x axis) and \n X2 (diff. line types)")



			if(both) {



				interaction.plot(temp[[j]], temp[[k]],



					y, xlab = locall[j + 2.], 



					trace.label = locall[k + 2.],



					ylab = locall[2.])



				title("Mean of Y for values of X2 (x axis) and \n X1 (diff. line types)")



			}



		}

			

	}



}


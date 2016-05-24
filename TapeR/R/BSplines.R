BSplines <-
function(knots=c(seq(0,1,0.1)), ord = 4, der = 0, x = c(seq(0,1,0.01)), ...){
#   ***********************************************************************************************

		## require(splines) 	# splineDesign {splines}

		TK = TransKnots(knots=knots,ord=ord)
		BSplines = splineDesign(knots=TK, x=x, ord = ord, derivs=c(rep(der,length(x))), outer.ok = T)

    }

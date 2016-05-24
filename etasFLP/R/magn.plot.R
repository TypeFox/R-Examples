magn.plot <-
function(catalog,...){mfreq=MLA.freq(catalog$magn1)
                        plot(mfreq$x,log(mfreq$back.cum.rel),
			    xlab="magnitude", ylab="Log-number events exceeding a magnitude value",
			    main="Transformed plot of magnitude frequencies",...)
                        lines(mfreq$x,log(mfreq$back.cum.rel) ,...)
	
			    }

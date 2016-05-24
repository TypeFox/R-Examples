"ex" <-

function(fun)



{



	temp <- substitute(fun)



	message2(get(paste(temp, ".ex", sep = "")))



	invisible()



}


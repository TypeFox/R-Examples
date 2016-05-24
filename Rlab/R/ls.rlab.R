"ls.rlab" <-

function(what="data")



{



        if (what == "all") what = "a"

        if (what == "data") what = "d"

        if (what == "function") what = "f"

        if (what == "functions") what = "f"

        if (what == "ex") what = "e"

        if (! (what == "a" || what == "d" || what == "f" || what == "e")) {



		stop("Accepted arguments are: data(default), functions, ex or all")



	}



	dpos <- NA

	fpos <- NA



	for(k in 1.:length(search())) {



		if(exists("Rlab.version", k, inherits=FALSE))



			dpos <- k



		if(exists("ls.rlab", k, inherits=FALSE))



			fpos <- k



	}



	if(is.na(dpos) && ((what == "d") ||(what == "e"))) {



		stop("R lab datasets have not been loaded")



	}



	if(is.na(dpos) && (what == "a")) {



		cat("Warning: R lab datasets have not been loaded\n")

		

		objects(pos=fpos)



	}



	else {

	

		if (what == "a") {

		

			if (dpos==fpos) objects(pos=dpos)

			

			else c(objects(pos=dpos),objects(pos=fpos))

			

			}

		

		else {



			if (what == "f") objs <- objects(pos=fpos)



			else objs <- objects(pos=dpos)



			ds <- c()



			for(i in 1.:length(objs)) {



				if ((! is.function(get(objs[i])) && ! is.character(get(objs[i])) && (what == "d")) ||

				    (is.function(get(objs[i])) && (what == "f")) ||

				    (length(grep(".ex",objs[i],fixed=TRUE))==1 && is.character(get(objs[i])) && (what == "e")))



					ds <- c(ds, objs[i])



			}



			ds



		}



	}



}


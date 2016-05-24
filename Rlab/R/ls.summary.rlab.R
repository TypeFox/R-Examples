"ls.summary.rlab" <-

function(what="data")



{



        if (what == "all") what = "a"

        if (what == "data") what = "d"

        if (what == "function") what = "f"

        if (what == "functions") what = "f"

        if (what == "ex") what = "e"



	objs <- c()

	

	objs <- ls.rlab(what)



	if (length(objs)==0) return()

	

	ds <- c()



	cls <- c()



	mde <- c()



	sze <- c()



	for(i in 1.:length(objs)) {



		ds <- c(ds, objs[i])



		cls <- c(cls, class(get(objs[i])))



		mde <- c(mde, mode(get(objs[i])))



		if (is.data.frame(get(objs[i])))



			sze <- c(sze, paste(dim(get(objs[i]))[1],'x',dim(get(objs[i]))[1]))



		else if (is.character(get(objs[i])))



			sze <- c(sze, paste(sum(nchar(get(objs[i]))),' ch (',length(get(objs[i])),' ln)',sep=''))



		else if (is.function(get(objs[i])))



			sze <- c(sze, paste(length(body(get(objs[i]))),'commands'))



		else

			

			sze <- c(sze, length(get(objs[i])))



	}



	result <- data.frame(ds,cls,mde,sze)



	names(result) <- c('Dataset','Class','Mode','Dimension')



	result



}


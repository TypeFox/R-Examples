"doConstrSuper" <-
function (X, Xlist, multimodel, thetalist, group)  
{   
    dset <- group[[1]][2]
    m <- multimodel@modellist[[dset]]
    colnames(X) <- m@compnames 
    if(multimodel@getXsuper) 
	X <- getXsuper(Xlist, multimodel@modellist, thetalist, group)
    t <- thetalist[[dset]]
    if(m@clpdep) {
	X <- doClpConstr(X,  clp_ind = group[[1]][1], 
                  	 clpCon = m@clpCon, clpequ = t@clpequ, 
		  	 num_clpequ = length(m@clpequspec), 
			 usecompnames0 = m@usecompnames0, 
			 usecompnamesequ = m@usecompnamesequ)
    }
    X
}


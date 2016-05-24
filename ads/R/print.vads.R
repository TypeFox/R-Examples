print.vads<-function(x,...) {
	UseMethod("print.vads")
}

print.vads.dval<-function(x,...) {
	cat("First-order local density values:\n")
	str(x)
	#cat("class: ",class(x),"\n")
    #cat("call: ")
	#print(x$call)
	#cat("sampling window :",x$window$type,"\n")
	#cat("\n")
	#sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	#sumry[1, ] <- c("$r", length(x$r), mode(x$r), "distance (r)")
	#class(sumry) <- "table"
    #print(sumry)
	#cat("\n")
	#sumry <- array("", c(3, 4), list(1:3, c("matrix", "nrow", "ncol", "content")))
    #sumry[1, ] <- c("$grid", nrow(x$grid), ncol(x$grid), "(x,y) coordinates of the sampling points A")
	#sumry[2, ] <- c("$count", nrow(x$count), ncol(x$count), "counting function NA(r)")
	#sumry[3, ] <- c("$density", nrow(x$dens), ncol(x$dens), "local density function nA(r)")
	#class(sumry) <- "table"
    #print(sumry)
}

print.vads.kval<-function(x,...) {
	cat("Univariate second-order local neighbourhood values:\n")
	str(x)
	#cat("class: ",class(x),"\n")
    #cat("call: ")
	#print(x$call)
	#cat("\n")
	#sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	#sumry[1, ] <- c("$r", length(x$r), mode(x$r), "distance (r)")
	#class(sumry) <- "table"
    #print(sumry)
	#cat("\n")
	#sumry <- array("", c(5, 4), list(1:5, c("matrix", "nrow", "ncol", "content")))
    #sumry[1, ] <- c("$coord", nrow(x$coord), ncol(x$coord), "(x,y) coordinates of points i")
	#sumry[2, ] <- c("$gi", nrow(x$gi), ncol(x$gi), "individual pair density values gi(r)")
	#sumry[3, ] <- c("$ni", nrow(x$ni), ncol(x$ni), "individual local neighbour density values ni(r)")
	#sumry[4, ] <- c("$ki", nrow(x$ki), ncol(x$ki), "individual Ripley's values Ki(r)")
	#sumry[5, ] <- c("$li", nrow(x$li), ncol(x$li), "modified individual Ripley's values Li(r)")
	#class(sumry) <- "table"
    #print(sumry)
}

print.vads.k12val<-function(x,...) {
	#verifyclass(x,"k12ival")
	cat("Bivariate second-order local neighbourhood values:\n")
	str(x)
	#cat("class: ",class(x),"\n")
    #cat("call: ")
	#print(x$call)
	#cat("mark1: ",x$marks[1],"\n")
	#cat("mark2: ",x$marks[2],"\n")
	#cat("\n")
	#sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	#sumry[1, ] <- c("$r", length(x$r), mode(x$r), "distance (r)")
	#class(sumry) <- "table"
    #print(sumry)
	#cat("\n")
	#sumry <- array("", c(5, 4), list(1:5, c("matrix", "nrow", "ncol", "content")))
    #sumry[1, ] <- c("$coord1", nrow(x$coord1), ncol(x$coord1), "(x,y) coordinates of points of mark 1")
	#sumry[2, ] <- c("$g12i", nrow(x$g12i), ncol(x$g12i), "individual pair density values g12i(r)")
	#sumry[3, ] <- c("$n12i", nrow(x$n12i), ncol(x$n12i), "individual local neighbour density values n12i(r)")
	#sumry[4, ] <- c("$k12i", nrow(x$k12i), ncol(x$k12i), "individual intertype values K12i(r)")
	#sumry[5, ] <- c("$l12i", nrow(x$l12i), ncol(x$l12i), "modified intrtype Ripley's values L12i(r)")
	#class(sumry) <- "table"
    #print(sumry)
}


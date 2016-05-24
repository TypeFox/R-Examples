summary.HTSClusterWrapper <-
function (object, ...) 
{
	x <- object
    	if (class(x) != "HTSClusterWrapper") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSClusterWrapper"), sep = ""), sep = "")
              }
	cat("*************************************************\n")
	cat("Selected number of clusters via ICL = ", ncol(x$ICL.results$lambda), "\n", sep = "")
	cat("Selected number of clusters via BIC = ", ncol(x$BIC.results$lambda), "\n", sep = "")
	cat("Selected number of clusters via Djump = ", ncol(x$Djump.results$lambda), "\n", sep = "")
	cat("Selected number of clusters via DDSE = ", ncol(x$DDSE.results$lambda), "\n", sep = "")
	cat("*************************************************\n")
      }


# a wrapper to get only the statistic and p value from the
# friedman.test function.

lw.FriedmanTest<-function(x) {
 ft<-friedman.test(x)
 return(c(ft$statistic,ft$p.value))
}

# calls subsidiary functions to fill the matrix of ranks x
# with values consistent or inconsistent with the existing
# complete rows of x.
# Return value:
# a list including:
# the name of the test used to determine the maximally consistent
# and inconsistent matrices created during the imputation process.
# two lists, each of which includes:
# the indices of the matrix that was selected
# the completed matrix (consistent or inconsistent)
# the statistic and p value of the test used in the screening

lwscreen<-function(x,scrtest="lw.FriedmanTest") {
 sum.is.na<-function(x) return(sum(is.na(x)))
 maxx<-x
 maxxna<-sum(is.na(x))
 while(maxxna > 0) {
  maxx<-listBuilder(maxx,fillArows,TRUE)
  # look for NAs
  maxxna<-listCrawler(maxx,sum.is.na)$value
 }
 maxxstat<-listCrawler(maxx,scrtest)
 minx<-x
 minxna<-sum(is.na(x))
 while(minxna > 0) {
  minx<-listBuilder(minx,fillArows,FALSE)
  # look for NAs
  minxna<-listCrawler(minx,sum.is.na)$value
 }
 minxstat<-listCrawler(minx,scrtest)
 lwstat<-list(test=unlist(strsplit(scrtest,"\\."))[2],
  maxx=maxxstat,minx=minxstat)
 class(lwstat)<-"lwstat"
 return(lwstat)
}

print.lwstat<-function(x,...) {
 cat("Lim-Wolfe test for incomplete rankings using",x$test,"\n")
 cat("Maximal consistency rank matrix\n")
 print(x$maxx$element)
 cat("Maximal consistency statistic = ",x$maxx$value[1],
  ", p = ",x$max$value[2],"\n",sep="")
 cat("\nMinimal consistency rank matrix\n")
 print(x$minx$element)
 cat("Minimal consistency statistic = ",x$minx$value[1],
  ", p = ",x$minx$value[2],"\n",sep="")
}

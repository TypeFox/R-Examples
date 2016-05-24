.joda.times <- matrix(c(
    # Days
    '^[0-9]{4}[0-9]{2}([-]|[ ]|[/])?[0-9]{2}$', '^yyyy([-]|[ ]|[/])?MM([-]|[ ]|[/])?dd$', '1988-03-02',
    '^[0-9]{2}([-]|[ ]|[/])[0-9]{2}([-]|[ ]|[/])?[0-9]{4}$', '^dd([-]|[ ]|[/])?MM([-]|[ ]|[/])?yyyy$', '02-03-1988',
    '^[0-9]{4}([-]|[ ]|[/])[A-Za-z]{3}([-]|[ ]|[/])?[0-9]{2}$', '^yyyy([-]|[ ]|[/])?MMM([-]|[ ]|[/])?dd$', '1988-Mar-02',
    '^[0-9]{2}([-]|[ ]|[/])[A-Za-z]{3}([-]|[ ]|[/])?[0-9]{4}$', '^dd([-]|[ ]|[/])?MMM([-]|[ ]|[/])?yyyy$', '02-Mar-1988',    
    # Months
    '^[0-9]{4}([-]|[ ]|[/])?[0-9]{2}$', '^yyyy([-]|[ ]|[/])?MM$', '1988-03',
    '^[0-9]{4}([-]|[ ]|[/])?[A-Za-z]{3}$', '^yyyy([-]|[ ]|[/])?MMM$', '1988-Mar',
    '^[A-Za-z]{3}([-]|[ ]|[/])?[0-9]{4}$', '^MMM([-]|[ ]|[/])?yyyy$', 'Mar-1988',
    '^[0-9]{2}([-]|[ ]|[/])?[0-9]{4}$', '^MM([-]|[ ]|[/])?yyyy$', '03-1988',
    # Years
    '^[0-9]{4}$', '^yyyy$', '1998'
  ), ncol=3, byrow=T
)
colnames(.joda.times) <- c('regex','format','example')
.joda.times <- as.data.frame(.joda.times, stringsAsFactors=F)

checkTimeFormat <- function(fmt) {
################################################################################
# Checks if the specified timeFormat is supported by DSPL
################################################################################
  result <- lapply(
    .joda.times[,2], 
    function(X,Y) {grep(pattern=X,x=Y)}, 
    Y=fmt
    )
  result <- sum(unlist(result)) == 1
  return(result)
}
#checkTimeFormat("yyyy")
#checkTimeFormat("yyyyMM")
#checkTimeFormat("yyyy-MM")
#checkTimeFormat("yyyy/MM")
#checkTimeFormat("yyyy, MM")
#.joda.times <- as.data.frame(joda.times, stringsAsFactors = F)
# Example vector
# x <- c(1998, '2010-01', '2010-01-01')
#test <- c('02-03-2001', '1992', 'Ene-2008')

#lapply(.joda.times[,1], function(X,Y) {grep(pattern=X,x=Y,value=T)}, Y=test)


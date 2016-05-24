
#################################################################
summary.mipo.nmi <- function(object, digits=4 ,  ...) {
	# adapted from summary.mipo from mice package
    x <- object
    table <- array( x$qbar, dim = c(length(x$qbar), 10) )
    dimnames(table) <- list(labels(x$qbar), 
			c("est", "se", "t", "df", "Pr(>|t|)", "lo 95", "hi 95", 
					"fmi" , "fmi_Betw" , "fmi_Within"))
    table[, 2] <- sqrt( diag(x$Tm) )
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- x$df
    table[, 5] <- if (all(x$df > 0)) 
        2 * (1 - stats::pt(abs(table[, 3]), x$df)) else NA
    table[, 6] <- table[, 1] - stats::qt(0.975, x$df) * table[, 2]
    table[, 7] <- table[, 1] + stats::qt(0.975, x$df) * table[, 2]
#    if (is.null(x$nmis) | is.null(names(x$qbar)))
#        table[, 8] <- NA else table[, 8] <- x$nmis[names(x$qbar)]
    table[, "fmi"] <- x$lambda
    table[, "fmi_Betw"] <- x$lambda_Between
	table[, "fmi_Within"] <- x$lambda_Within
	table <- as.data.frame(table)
	if ( is.na(table$se)[1] ){
			table$df <- NA
						}
	if ( ! is.null( object$u_NULL ) ){
		if ( object$u_NULL){
			table <- table[ , "est" , drop=FALSE ]		
						}
					}
						
						
	table0 <- table
	for (vv in seq(1 , ncol(table) ) ){
		table[,vv] <- round( table[,vv] , digits=digits )
							}
    print(table)
	invisible(table0)
}
#################################################################
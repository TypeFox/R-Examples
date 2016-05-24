AT.SPC.tapply <- function( spc, 
                           INDEX, 
                           FUN, 
                           mixed.field.arguments = list(E.MeV.u = "E.mid.MeV.u", fluence.cm2 = "N.per.primary", particle.no = "particle.no"),
						   additional.arguments = NULL, 
                           names.results = NULL)
{    
    ##############################
    # Get index columns and levels
    index.columns    <- which(is.element(names(spc), INDEX))
    if(length(INDEX) != length(index.columns)){
        cat("At least one index variable not found in spc data.\n")
        return(NULL)
    }

    index.variable   <- NULL
    for (i in 1:length(index.columns)){
        # DEBUG: i <- 1
        index.variable    <- paste(index.variable, spc[,index.columns[i]])
    }
    levels           <- unique(index.variable)

    ###########################
    # Get argument list for FUN
	# Match with mixed field args
    args.FUN         <- names(formals(FUN))
    args.list        <- "("
	for (i in 1:length(args.FUN)) {
		mixed.field.arguments.idx <- match(args.FUN[i], names(mixed.field.arguments))
		if (!is.na(mixed.field.arguments.idx)) {
			args.list <- paste(	args.list, 
								args.FUN[i], " = spc$", 
								mixed.field.arguments[[mixed.field.arguments.idx]], "[ii],", 
								sep = "")
		}
	}
	
    if(!is.null(additional.arguments)){
        for(j in 1:length(additional.arguments)){
             if(additional.arguments[[j]][3] == TRUE){
                 args.list    <- paste( args.list, 
                                        additional.arguments[[j]][1], 
                                        " = spc$",
                                       additional.arguments[[j]][2],
                                       "[ii],",
                                       sep = "")
             }else{
                 args.list    <- paste( args.list, 
                                        additional.arguments[[j]][1], 
                                        " = ",
                                       additional.arguments[[j]][2],
                                       ",",
                                       sep = "")
             }
        }
    }
    args.list        <- paste(substring(args.list, 1, nchar(args.list) - 1), ")")

    df.return        <- NULL
    for(cur.level in levels){
        # DEBUG: cur.level <- levels[2]
        ii            <- index.variable == cur.level

        res           <- eval( parse( text = paste( "FUN",
                                                    args.list,
                                                    sep = "")))
        df.cur.level  <-  cbind.data.frame( unique(data.frame( spc[ii,index.columns])), 
                                            res)
        if(is.null(df.return)){
            df.return    <- df.cur.level
        }else{
            df.return    <- rbind.data.frame( df.return,
                                              df.cur.level)
        }
    }
    row.names(df.return) <- 1:nrow(df.return)
    names(df.return)[1:length(index.columns)]  <- INDEX
    if(!is.null(names.results)){
        names(df.return)     <- c(names(df.return)[1:length(index.columns)], names.results)
    }
    return(df.return)
}
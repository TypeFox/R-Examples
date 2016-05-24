#.sortWidget----------------------------2012-12-20
.sortWidget <- function( d, hisname )
{
	#given a tcl_variable that represents a matrix, fetch it and convert into an R matrix
	.sortWidgetGetVal <- function( tcl_var )
	{
		tmp <- tclvalue( tcl("array", "get", tcl_var) )
		tmp <- .tclArrayToVector( tmp )
		
		row_max <- 0
		col_max <- 0
		for( i in seq(from=1, to=length(tmp), by=2 ) ) {
			key <- as.integer( strsplit(tmp[i], ",")[[1]] )
			row_max <- max( row_max, key[1] )
			col_max <- max( col_max, key[2] )
		}
		ret <- matrix( nrow=row_max, ncol=col_max )
		for( i in seq(from=1, to=length(tmp), by=2 ) ) {
			key <- as.integer( strsplit(tmp[i], ",")[[1]] )
			ret[ key[1], key[2] ] <- tmp[i+1]
		}
		return( ret )
	}

	#called when user pushes save button
	.done_sorting <- function( r, hisname, tt )
	{
		sorted_data <- .sortWidgetGetVal( r )
		sorted_index <- as.integer( sorted_data[,1] )
		tget(PBS.history)
		tmp <- PBS.history[[hisname]]
		j <- 2
		#eval(parse(text=".PBSmodEnv$PBS.history[[hisname]] <<- .PBSmodEnv$PBS.history[[hisname]][1]"))
		PBS.history[[hisname]] <- PBS.history[[hisname]][1]
			for (i in sorted_index) {
			if (!is.na(i)) {
				#eval(parse(text=".PBSmodEnv$PBS.history[[hisname]][[j]] <<- tmp[[i + 1]]"))
				PBS.history[[hisname]][[j]] <- tmp[[i + 1]]
				j <- j + 1
			}
		}
		tput(PBS.history)
		tkdestroy( tt )
		jumpHistory( hisname, 1 )
	}

	if( is.data.frame( d ) == FALSE ) stop( "d must be a data.frame" )
	r <- tclArray()
	for( i in 1:nrow(d) ) {
		#display initial row position as first column
		tcl( "set", paste(r,"(",i,",",1,")", sep=""), i )
		for( j in 1:ncol(d) )
			tcl( "set", paste(r,"(",i,",",j+1,")", sep=""), d[i,j] )
	}

	tt <- tktoplevel()
	#frame <- tkframe( tt )
	#tkgrid( tklabel( frame, text="Click and drag row items to adjust the history order" ), row=1, column=1, sticky="W" )
	#tkgrid( tkbutton( frame,text="Save Order",command=function() print("done!") ), row=1, column=2, sticky="E" )
	#tkpack( frame, fill="both", expand=0 )
	frame <- tkframe( tt )
	tkpack( tklabel( frame, text="Click and drag row items to adjust the history order" ), side="left" )
	tkpack( tkbutton( frame, text="Save Order",command=function() .done_sorting( r, hisname, tt ) ), side="right" )
	tkpack( frame )

	p <- tcl( "PBSmodelling::create", .Tk.subwin( tt ), r, paste( nrow(d), ncol(d)+1 ), c("history index", colnames(d) ) )
	tkpack( p, expand=1, fill="both" )
	tkfocus(tt)
	return( r )
}
#--------------------------------------.sortWidget


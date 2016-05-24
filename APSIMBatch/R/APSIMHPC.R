# * Author:    Bangyou Zheng (Bangyou.Zheng@csiro.au)
# * Created:   14/04/2011
# *
# * $Revision: 2369 $
# * $Id: APSIMHPC.R 2369 2011-07-19 10:42:28Z zhe00a $
# * $Author: zhe00a $
# * $Date: 2011-07-19 20:42:28 +1000 (Tue, 19 Jul 2011) $

#' Run APSIM with HPC
#' 
#' @param ... A list contain several simulations, 
#' or a list contain sim or RData files, See example for more 
#' detailed arguments.
#' @param apsim The relative and absolute path to apsim.exe
#' @param extra A data frame of extra identification for all simulation.
#' Row number must be the same as simulation number 
#' @return A list contain all simulation results.
#' @examples
#' # NO RUN
#' # Run all sim files
#' \dontrun{files <- list.files( "./simtest", pattern = "(.*)(\\.sim$)", full.names = TRUE )}
#' \dontrun{runAPSIM( files = files )}
#' # Run Rdata file
#' \dontrun{runAPSIM( files = "test.RData" )}
#' # Set apsim path
#' \dontrun{apsim <- "\"C:/Program Files/Apsim71/Model/apsim.exe\""}
#' \dontrun{runAPSIM( files = "test.RData", apsim = apsim )}

runAPSIM <- function( ..., extra = NULL, apsim = "../Apsim71/Model/apsim.exe" )
{
	all_sim <- as.list( NULL )
	other_args <- list( ... )
	if ( !is.null( other_args$sims ) )
	{
		all_sim <- other_args$sims  
	}
	if ( !is.null( other_args$files ) )
	{
		files <- other_args$files 
		for ( i in seq( along = files ) )
		{
			file_name <- basename( files[i] )
			file_ext <- tolower( gsub( "(.*)(\\..*$)", "\\2", file_name ) )
			file_name <- tolower( gsub( "(.*)(\\..*$)", "\\1", file_name ) )
			if ( file_ext == ".sim" )
			{
				all_sim[[file_name]] <- readLines( files[i], warn = FALSE )
			} else if (  file_ext == ".rdata" )
			{
				temp_env <- new.env()
				sim_vars <- load( files[i], envir = temp_env )
				for ( j in seq( along = sim_vars ) )
				{
					all_sim <- c( all_sim, get( sim_vars[j], envir = temp_env ) )
				}
			}
		}
		
	}

	sim_names <- names( all_sim )
	if ( !is.null( extra ) && nrow( extra ) != length( sim_names ) )
	{
		stop( "Row number of extra must be the same as simulation number." )
	}
	extra <- as.data.frame( extra )
	all_sim_res <- as.list( NULL )
	for ( j in seq( all_sim ) )
	{
		
		new_apsim <- all_sim[[j]]
		sim_file_name <- paste( sim_names[j], ".sim", sep = "" )
		write.table( new_apsim, file = sim_file_name, 
				quote = FALSE,
				row.names = FALSE,
				col.names = FALSE )
		# RUN APSIM to get results
		system( paste( apsim, " \"", sim_file_name, "\"", sep = "" ), 
				show.output.on.console = FALSE )
		
		# Count number of report
		c_sim_name <- findElement( "title", new_apsim )
		c_sim_name <- getElementVlaue( new_apsim[c_sim_name] )
		reports <- new_apsim[grep( "executable=.*Report.dll", new_apsim )]
		reports <- gsub( "(.*name=\")(.*)(\".*executable.*$)", "\\2", reports )
		reports_file <- paste( c_sim_name, " ", reports, ".out", sep = "" )
		sim_sim <- as.list( NULL )
		for ( i in seq( along = reports ) )
		{
			temp <- readLines( reports_file[i], n = 100 )
			temp <- temp[grep( "=", temp )]
			report_res <- as.data.frame( extra[j,] )
			names( report_res ) <- names( extra )
			for ( m in seq( along = temp ) )
			{
				com <- splitEqual( temp[m] )
				report_res[[com$name]] <- com$value
			}
			report_sim <- read.table( reports_file[i],
					sep = "", header = FALSE, as.is = TRUE, 
					skip = length( temp ) + 2,
					col.names = scan( reports_file[i], "", 
							sep = "", skip = length( temp ), 
							nline = 1 ) )
			report_res <- cbind( as.data.frame( report_res ), report_sim )
			sim_sim[[reports[i]]] <- report_res
		}
		all_sim_res[[sim_names[j]]] <- sim_sim
		# remove output files
		file.remove( c( sim_file_name, reports_file ) )
	}
	return( all_sim_res )
}


#' Generate simulations according several factors
#' 
#' @param template File path to template sim files 
#' @param factors A data frame which contained all combinations of factors.
#' The column names of data frame are parameter names in the sim file.
#' @return A list which contain all simulations. Row names of factors 
#' are used for simulation names.
generateSim <- function( template, factors )
{
	sim_names <- row.names( factors )
	factor_names <- names( factors )
	# read template files
	template_sim <- readLines( template, warn = FALSE )
	all_sim <- NULL
	for ( i in seq( length= nrow( factors ) ) )
	{
		new_sim <- template_sim
		
		# replace simulation names
		sim_title <- grep("(^.*simulation name.*=.*\")(.*)(\".*executable.*$)", 
				new_sim )
		new_sim[sim_title] <- sub( "(^.*simulation name.*=.*\")(.*)(\".*executable.*$)", 
				paste( "\\1", sim_names[i], "\\3", sep = "" ), 
				new_sim[sim_title] )
		sim_title <- grep("<title>.*</title>", new_sim )
		new_sim[sim_title] <- gsub( "<title>.*</title>", 
				paste( "<title>", sim_names[i], "</title>", sep = "" ), 
				new_sim[sim_title] )
		# replace each factor
		c_factor <- removeAttribure( factors[i,] )
		replace_num <- 0
		for ( j in seq( along = factor_names ) )
		{
			#<date type="text" description="Enter sowing date (dd-mmm) : ">16-apr</date>
			factor_row <- findElement( factor_names[j], new_sim ) 
			
			if ( length( factor_row ) == 0 )
			{
				warning( paste( "Can not find factor for ", c_factor[j], " in ",
								factor_names[j], ". Skip it.", sep = "" ) )
				next
			}
			
			# testing whether cultivar parameters
			base_cul_row <- grep( "<base_cultivar>", new_sim[1:factor_row[1]] )
			if ( length( base_cul_row ) == 0 )
			{
				if ( length( factor_row ) > 1 )
				{
					warning( paste( "There are several rows found for ", 
									c_factor[j], " in ",
									factor_names[j], 
									". Only first row is used.", sep = "" ) )
				}
				
				new_sim[factor_row[1]] <- replaceElementVlaue( c_factor[j], 
						new_sim[factor_row[1]] )
						
			} else
			{ # For cultivar parameters
				# find cultivar 
				cultivar <- findElement( "cultivar", new_sim ) 
				cultivar <- sub( "(^.*<cultivar.*>)(.*)(</cultivar>.*$)", 
						"\\2", new_sim[cultivar] )
				
				# Find starting and end row of cultivar
				cultivar_s_row <- grep( paste( "^.*<", cultivar,
								" .*cultivar.*=.*\"yes\".*>.*$", 
								sep = ""), new_sim )
				cultivar_e_row <- grep( paste( "^.*</", cultivar, ">.*$", 
								sep = ""), new_sim )
				if ( length( cultivar_s_row ) == 0 | 
						length( cultivar_e_row ) == 0 )
				{
					warning( paste( "Can not find parameters for Cultivar ", 
									cultivar, ". Skip it.", sep = "" ) )
					next
				}
				cul_xml <- new_sim[cultivar_s_row:cultivar_e_row]
				cul_factor_row <- findElement( factor_names[j], cul_xml )
				
				if ( length( cultivar_s_row ) > 0 )
				{ # Find attribute for this cultivar
					new_sim[cultivar_s_row + cul_factor_row - 1 ] <- 
							replaceElementVlaue( c_factor[j], cul_xml[cul_factor_row] )
				} else
				{ # Not find attribute for this cultivar, add a new row
					new_row <- paste( "          <", factor_names[j], ">", 
							c_factor[j], "</", factor_names[j], 
							">", sep = "" )
					new_sim <- c( new_sim[1:(cultivar_e_row-1)], new_row, 
							new_sim[(cultivar_e_row):length(new_sim)] )
				}
			}
			replace_num <- replace_num + 1
		}
		if ( replace_num == 0 )
		{
			warning( "There are nothing to generate new simulation. Skip it.")
		} else
		{
			all_sim[[sim_names[i]]] <- new_sim
		}
	}
	return( all_sim )
}


#' Find a row for element name
#' @param name name
#' @param xml xml
findElement <- function( name, xml )
{
	pos <- grep( paste( "(^.*<", name,
					")(>| .*>).*(</", name, ">.*$)",
					sep = "" ), xml )
	return( pos )
}


#' Replace a row for element name
#' @param value value
#' @param xml xml
replaceElementVlaue <- function( value, xml )
{
	value <- as.character( unlist( value ) )
	value <- gsub( "\\\\", "\\\\\\\\", value )
	newvalue <- sub( paste( "(^.*<.*>).*(</.*>.*$)", sep = "" ),
			paste( "\\1", value, "\\2", sep = "" ), xml )
	return( newvalue )
}

#' Remove all attribure of a vector
#' @param x x
removeAttribure <- function( x )
{
	x <- as.list( x )
	res <- NULL
	for ( i in seq( along = x ) )
	{
		res <- c( res, as.character( x[[i]] ) )
	}
	return( res )
}

#' Replace a row for element name
#' @param xml xml
getElementVlaue <- function( xml )
{
	newvalue <- sub( paste( "(^.*<.*>)(.*)(</.*>.*$)", sep = "" ),
			"\\2", xml )
	return( newvalue )
}

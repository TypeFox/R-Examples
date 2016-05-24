readXMLLines <- function(lines, line.init=1){

	rlist <- list()
	prev_line_type <- ''

	# READ EACH LINE
	i <- line.init
	while(i <= length(lines)){
	
		# SET LINE
		line <- lines[i]

		#cat(i, ':', line, '\n')
		
		# SKIP EMPTY LINES
		if(line == ""){
			i <- i + 1
			next
		}
		
		# CHECK FOR OBJECT OPEN/CLOSE TAG
		reg_expr_open <- regexpr(pattern='<[[:alnum:]|_.]+[ |>]', text=line)
		reg_expr_close <- regexpr(pattern='</[ ]*[[:alnum:]|_.]+[ ]*>', text=line)

		# SET LINE TYPE
		line_type <- 'fill'
		if(reg_expr_open > 0) line_type <- 'object_open'
		if(reg_expr_close > 0) line_type <- 'object_close'

		# GET OBJECT DATA
		if(line_type == 'object_open'){
			
			if(prev_line_type == 'object_open') {

				read_xml_2_list <- readXMLLines(lines, line.init=i)

				if(is.null(rlist)) rlist <- list()
				if(is.matrix(rlist) && nrow(rlist) == 0) rlist <- list()

				rlist[[read_xml_2_list$obj.name]] <- read_xml_2_list$rlist
				#print(rlist)
				
				i <- read_xml_2_list$line.init
				next

			}

			# GET OBJECT NAME
			object_name <- gsub("<|[\t]", "", substr(line, reg_expr_open, reg_expr_open+attr(reg_expr_open,"match.length")-2))

			# CHECK FOR OBJECT TYPE
			reg_expr <- regexpr(pattern='type=[\'|\"]?[[:alnum:]|.]+[\'|\"]?', text=line)

			if(reg_expr > 0){
				object_type <- gsub('[\'|\"]', "", substr(line, reg_expr+5, reg_expr+attr(reg_expr,"match.length")-1))
			}else{
				object_type <- 'vector'
			}
			
			# CHECK FOR AS NUMERIC
			if(grepl(pattern='as[.]?numeric=[\'|\"]?true[\'|\"]?', x=line, ignore.case=TRUE)){is_numeric <- TRUE}else{is_numeric <- FALSE}

			# CREATE EMTPY VECTOR
			if(object_type == 'vector') rlist <- c()

			if(object_type == 'matrix'){

				# GET NUMBER OF TABS IN FIRST ROW
				ntab <- length(strsplit(x=lines[i+1], split='\t')[[1]]) - 1

				# FIND IF INCLUDES ROWNAMES
				if(grepl(pattern='rownames=[\'|\"]?true[\'|\"]?', x=line, ignore.case=TRUE)){add_rownames <- TRUE}else{add_rownames <- FALSE}

				# FIND IF INCLUDES COLNAMES
				if(grepl(pattern='colnames=[\'|\"]?true[\'|\"]?', x=line, ignore.case=TRUE)){add_colnames <- TRUE}else{add_colnames <- FALSE}

				if(add_rownames){

					# CREATE EMPTY VECTOR FOR ROWNAMES
					row_names <- c()

					ncol <- ntab
				}else{
					ncol <- ntab + 1
				}

				if(add_colnames){
					ncol <- ntab + 1
					colnames_saved <- FALSE
				}else{
					colnames_saved <- TRUE
				}

				# CREATE EMTPY MATRIX
				rlist <- matrix(NA, nrow=0, ncol=ncol)
			}

			if(object_type == 'list'){

				# FIND IF INCLUDES NAMES
				if(grepl(pattern='names=[\'|\"]?true[\'|\"]?', x=line, ignore.case=TRUE)){add_names <- TRUE}else{add_names <- FALSE}

				# CREATE EMTPY LIST
				rlist <- list()
			}


		}else if(line_type == 'fill'){

			if(object_type == 'vector'){
			
				data_values <- line

				# SET AS NUMERIC
				if(is_numeric) data_values <- suppressWarnings(as.numeric(data_values))

				rlist <- c(rlist, data_values)

			}else if(object_type == 'matrix'){

				line_split <- strsplit(x=line, split='\t')[[1]]

				if(add_colnames && !colnames_saved){
					col_names <- line_split
					colnames_saved <- TRUE
				}else{
					if(add_rownames){
					
						# ADD ROWNAMES
						row_names <- c(row_names, line_split[1])
					
						# GET DATA VALUES
						data_values <- line_split[2:length(line_split)]
					}else{

						# GET DATA VALUES
						data_values <- line_split[1:length(line_split)]
					}
				
					# SET AS NUMERIC
					if(is_numeric) data_values <- suppressWarnings(as.numeric(data_values))

					# CHECK THAT NUMBER OF DATA VALUES MATCHES NUMBER OF COLUMNS IN MATRIX
					if(length(data_values) != ncol) stop(paste0("Inconsistent number of columns in row (line ", i, ")."))
				
					rlist <- rbind(rlist, data_values)
				}

			}else if(object_type == 'list'){
				
				line_split <- strsplit(x=line, split='\t')[[1]]
				
				if(add_names){
					
					# GET NAME
					name <- line_split[1]
					
					# GET DATA VALUES
					data_values <- line_split[2:length(line_split)]
				}else{

					# GET DATA VALUES
					data_values <- line_split[1:length(line_split)]
				}

				# SET AS NUMERIC
				if(is_numeric) data_values <- suppressWarnings(as.numeric(data_values))

				if(add_names){
					rlist[[name]] <- data_values
				}else{
					rlist[[length(rlist[[object_name]])+1]] <- data_values
				}
			}


		}else if(line_type == 'object_close'){

			if(is.matrix(rlist)){
			
				# ADD ROWNAMES
				if(add_rownames){
					rownames(rlist) <- row_names
				}else{
					rownames(rlist) <- NULL
				}

				if(add_colnames){
					colnames(rlist) <- col_names
				}else{
					colnames(rlist) <- NULL
				}
				
				if(nrow(rlist) == 0) rlist <- NULL
			}

			return(list(rlist=rlist, obj.name=object_name, line.init=i+1))
		}
		
		i <- i + 1
		prev_line_type <- line_type
	}
}
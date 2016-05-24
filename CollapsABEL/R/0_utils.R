
#' Create directory if it does not already exist
#' 
#' @param dir character. Path of directory to be created.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
dir.create2 = function(dir) {
	if(!file.exists(dir)) {
		dir.create(dir, recursive = TRUE)
	} else {
		TRUE
	}
}

#' Create file if it does not already exist
#' @param filename character. Path of file to be created.
#' @author Kaiyin Zhong, Fan Liu
#' @export
file.create2 = function(filename) {
	dir.create2(dirname(filename))
	if(!file.exists(filename)) {
		file.create(filename)
	} else {
		TRUE
	}
}

#' Check whether a file is a SQLite3 database.
#' 
#' @param filename character. Path to file to be checked.
#' @author Kaiyin Zhong, Fan Liu
#' @export
isSQLite3 = function(filename) {
	if(!file.exists(filename)) {
		return(FALSE)
	} 
	if(file.info(filename)$size < 100){
		return(FALSE)
	}
	con = file(filename, "rb")
	tryCatch({
				header = rawToChar(readBin(con, "raw", 15))
				return(header == "SQLite format 3")
			}, finally = {
				close(con)
			})
}


#' Distance with lag
#' 
#' Calculate the distance between each element in a numeric vector and 
#' the element that is \code{lag} positions after it. For the last 
#' \code{lag} elemnts, this distance does not exist, so NA is used as a placeholder.
#' The returned vector is of the same length as the input vector.
#' 
#' @param vec numeric.
#' @param lag integer.
#' @param reverse logical. Default to FALSE, i.e. calculate \code{vec[i+lag] - vec[i]}. When set to TRUE, calculate \code{vec[i] - vec[i+lag]}
#' @return numeric.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
lagDistance = function(vec, lag = 1, reverse = FALSE) {
	new_vec = c(vec[(1+lag):length(vec)], rep(NA, lag))
	stopifnot(length(new_vec) == length(vec))
	res = new_vec - vec
	if(reverse) 
		(-1 * res)
	else
		res
}




sqliteFileGcdh = function(tag, dbname) {
	fp = file.path(tag2Dir(tag, "gcdh"), 
			sprintf("%s.sqlite", dbname))
	file.create2(fp)
	fp
}

#' Call system command with format string
#' 
#' @param ... passed to \code{sprintf}
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
systemFormat = function(...) {
	system(sprintf(...))
}

#' Stop with format string
#' 
#' @param ... passed to \code{sprintf}
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
stopFormat = function(...) {
	stop(sprintf(...))
}

# like expand.grid, but reverse the columns order
expand.grid.rev = function(...) {
	ret = expand.grid(...)
	ret[, rev(colnames(ret))]
}

#' Correlation coefficient of column-pairs of two data frames
#' 
#' @param dat1 first data.frame
#' @param dat2 second data.frame
#' @return A vector of correlation coefficients.
#' 
#' @author Kaiyin Zhong
#' @export
colCors = function(dat1, dat2) {
	stopifnot(all(dim(dat1) == dim(dat2)))
	sapply(1:ncol(dat1), function(i) {
				cor(dat1[, i], dat2[, i])
			})
}

#' Default value for expression. 
#' 
#' When an expression evals to NULL, take the default value instead. Copied from dplyr source.
#' @param x expression to be evaled. 
#' @param y default value.
#' @name getOrElse-operator
#' @author Hadley Wickham
#' @export
"%||%" <- function(x, y) if(is.null(x)) y else x


#' Change extension names
#' 
#' @param filename character. File path 
#' @param ext_name character. New extension name
#' @importFrom tools file_path_sans_ext
#' 
#' @author Kaiyin Zhong
#' @export
chExt = function(filename, ext_name) {
	base_name = tools::file_path_sans_ext(filename)
	paste0(base_name, ".", ext_name)
}

# ctime of db must be later than mtime of plink files
dbUpToDate = function(dbname) {
	fam_file = chExt(dbname, "fam")
	bim_file = chExt(dbname, "bim")
	bed_file = chExt(dbname, "bed")
	frq_file = chExt(dbname, "frq")
	db_ctime = as.integer(file.info(dbname)$ctime)
	fam_mtime = as.integer(file.info(fam_file)$mtime)
	bim_mtime = as.integer(file.info(bim_file)$mtime)
	bed_mtime = as.integer(file.info(bed_file)$mtime)
	frq_mtime = as.integer(file.info(frq_file)$mtime)
	frqUpToDate(frq_file) && (db_ctime > fam_mtime) && (db_ctime > bim_mtime) && (db_ctime > bed_mtime) && (db_ctime > frq_mtime)
}

frqUpToDate = function(frq_file) {
	bed_file = chExt(frq_file, "bed")
	bed_ctime = as.integer(file.info(bed_file)$ctime)
	frq_ctime = as.integer(file.info(frq_file)$ctime)
	frq_ctime > bed_ctime
}


rmPlinkFiles = function(plink_stem) {
	unlink(
			paste0(plink_stem, 
					c(".bed", 
							".bim", 
							".fam", 
							".frq", 
							".nosex", 
							".log", 
							".assoc",
							".qassoc",
							".assoc.linear", 
							".assoc.logistic", 
							".sqlite" 
					))
	)
}

#' Clear up CollapsABEL workspace
#' 
#' The workspace folder is defined in \code{collenv$.collapsabel_dir}.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
collClear = function() {
	unlink(Sys.glob(file.path(collenv$.collapsabel_dir, "*", "*")), recursive = TRUE)
}


#' Print quoted expression then its value
#' 
#' @param expr expression to be evaluated.
#' @export 
eprint = function(expr) {
	message(substitute(expr))
	print(expr)
}

#' Generate a m by n data.frame from normal distribution
#' 
#' @param m integer. Number of rows.
#' @param n integer. Number of columns.
#' 
#' @author Kaiyin Zhong
#' @export
randNormDat = function(m, n) {
	as.data.frame(matrix(rnorm(m * n), m, n))
}


#' Read phenotype file
#' @param file character, path to phenotype file.
#' @return data.frame
#' 
#' @author kaiyin
#' @export
read.phe.table = function(file) {
	phe = read.table(file, header = TRUE)
	stopifnot(all(c("FID", "IID") %in% colnames(phe)))
	for(i in 1:ncol(phe)) {
		if(
				all(
						isBinary(phe[, i]) && sort(unique(na.omit(phe[, i]))) == c(0, 1)
				)
				) {
			phe[, i] = phe[, i] + 1
			# https://www.cog-genomics.org/plink2/input#pheno
			warning(sprintf("%s, column %d: plink will take 0 as missing in binary trait.", file, i))
		}
	}
	phe$FID = as.character(phe$FID)
	phe$IID = as.character(phe$IID)
	phe
}

#' Write a phenotype data.frame to file
#' 
#' @param phe data.frame
#' @param file character, path to phenotype file.
#' 
#' @author Kaiyin Zhong
#' @export
write.phe.table = function(phe, file) {
	stopifnot(all(c("FID", "IID") %in% colnames(phe)))
	write.table(phe, file = file, quote = FALSE, row.names = FALSE)
}

#' Collpase genotypes
#' 
#' @param g1 numeric, genotype vector 1.
#' @param g2 numeric, genotype vector 2.
#' @param collapse_matrix matrix of integers range from 0 to 3.
#' @return numeric, collapsed genotype of g1 and g2.
#' 
#' @author Kaiyin Zhong
#' @export
collapse = function(g1, g2, collapse_matrix = NULL) {
	if(is.null(collapse_matrix)) {
		collapse_matrix = matrix(c(
						0L, 0L, 0L, 0L, 
						0L, 1L, 1L, 1L, 
						0L, 1L, 0L, 3L, 
						0L, 1L, 3L, 3L
						), 4, 4)
	}
	convertToMachineCode = function(g) {
		g[g == 0] = 3 
		g[g == 2] = 0 
		g[g == 1] = 2 
		g[is.na(g)] = 1
		g
	}
	convertToGeno = function(code) {
		code[code == 1] = NA
		code[code == 2] = 1 
		code[code == 0] = 2 
		code[code == 3] = 0
		code
	}
	idx = cbind(convertToMachineCode(g1), convertToMachineCode(g2))
	code = collapse_matrix[idx+1]
	convertToGeno(code)
}

#' Collapse two genotype matrices, column by column
#' 
#' Each column is assummed to be the genotype for a SNP. The two genotype matrices should have the same size. 
#' 
#' @param m1 first genotype matrix
#' @param m2 second genotype matrix
#' @param collapse_matrix collapsed genotype matrix
#' @return collapsed genotyp matrix
#' 
#' @author kaiyin
#' @export
collapseMat = function(m1, m2, collapse_matrix = matrix(c(
						0L, 0L, 0L, 0L, 
						0L, 1L, 1L, 1L, 
						0L, 1L, 0L, 3L, 
						0L, 1L, 3L, 3L
						), 4, 4))  {
			if(!all(dim(m1) == dim(m2))) stop("m1 and m2 should have same size!")
			m = m1 
			for(i in 1:ncol(m1)) {
				m[, i] = collapse(m1[, i], m2[, i], collapse_matrix = collapse_matrix)
			}
			m
		}

#coll_mat = matrix(c(
#				0L, 0L, 0L, 0L, 
#				0L, 1L, 2L, 1L, 
#				0L, 2L, 0L, 2L, 
#				0L, 1L, 2L, 3L
#		), 4, 4)
#g1 = sample(0:2, 100, repl = TRUE)
#g2 = sample(0:2, 100, repl = TRUE)
#g1[sample(1:100, 15)] = NA
#g2[sample(1:100, 15)] = NA
#cbind(g1, g2, collapse(g1, g2, coll_mat))

#' Estimate percentage of variation explained
#' 
#' @param df Dataframe
#' @param yn Name of the independent variable, must be one of the columns of \code{df}
#' 
#' @author Fan Liu
#' @export
getr2 <- function(df,yn){
	if(!(yn %in% names(df))) {
		stop("Dependent var name not valid!")
	}
	
	# construct a names map
	old_names = colnames(df)
	new_names = paste0("x", 1:ncol(df))
	names_map = data.frame(old_names, new_names)
	df = setNames(df, new_names)
	
	old_yn = yn
	yn = new_names[old_names == yn]
	
	# xn: independent var names
	xn = names(df)[names(df) != yn]
	
	fo <- as.formula(paste(yn, "~",paste(xn,collapse="+")))
	# summary table, intercept excluded
	tb <- summary(glm(fo,data=df))$coef[-1,]
	# sort rownames by reverse order of t values
	xn <- rownames(tb)[order(-abs(tb[,3]))]
	# setup result
	r2 <- rep(NA,length(xn))
	names(r2) <- xn
	for (i in 1:length(xn)){
		fo <- as.formula(paste(yn, "~",paste(xn[1:i],collapse="+")))
		m <- lm(fo,data=df)
		r2[i] <- summary(m)$r.squared
	}
	r2 <- (r2-c(0,r2[-length(r2)]))*100
	
	df1 <- as.data.frame(tb)
	df1$Name<-rownames(tb)
	df2 <- data.frame(Name=names(r2),r2=r2)
	df3 <- merge(df1, df2, by="Name")
	df3 <- df3[order(-df3$r2),]
	
	rownames(names_map) = new_names
	df3$Name = names_map[df3$Name, ]$old_names
	df3
}


#' Extract one row or column of a data frame as a vector
#' @param dat data.frame
#' @param i row or column number
#' @param row Logical. If TRUE, then i is the row number, otherwise i is the column number
#' @return A vector.
#' 
#' @author kaiyin
#' @export
datToVec = function(dat, i, row=TRUE) {
	if(row) {
		t(dat[i, ])[, 1]
	} else {
		dat[, i]
	}
}
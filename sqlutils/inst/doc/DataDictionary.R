### R code from vignette source 'DataDictionary.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(sqlutils)
library(RSQLite)
library(xtable)

sqlfile <- paste(system.file(package='sqlutils'), '/db/students.db', sep='')
m <- dbDriver("SQLite")
conn <- dbConnect(m, dbname=sqlfile)
sqlPaths(paste(system.file(package='sqlutils'), '/sql', sep=''), replace=TRUE)
queries <- getQueries()
nlevels <- 10

sanitizeLatex <- function(str) {
	gsub('([#$%&~_\\^\\\\{}])', '\\\\\\1', str, perl=TRUE)
}


###################################################
### code chunk number 2: docs
###################################################
for(q in seq_along(queries)) {
	doc <- sqldoc(queries[q])

	cat('\\clearpage\n')
	
	cat(paste('\\section{', queries[q], '}\n', sep=''))
	cat(paste('\\label{', queries[q], '}\n\n', sep=''))
	
	if(!is.null(doc$introduction)) {
		cat(doc$introduction)
		cat("\n\n")
	}
	
	tryCatch( {
	
	sql <- suppressWarnings(getSQL(queries[q]))
	cat('\\subsection{Parameters}\n')
	if(length(getParameters(queries[q])) > 0) {
		cat('\\begin{description}\n')
		for(j in 1:nrow(doc$params)) {
			desc <- doc$params[j,]$desc
			desc <- ifelse(is.null(desc) | desc=='', ' not specified', desc)
			default <- doc$params[j,]$default
			default <- ifelse(is.null(default) | default=='', ' not specified', default)
			cat(paste('\\item[', doc$params[j,]$param, '] ', 
					  desc, 
					  '\n\n\\begin{lstlisting}\n', 
					  default, 
					  '\n\\end{lstlisting}\n', sep=''))
		}
		cat('\\end{description}\n\n')
	} else {
		cat('None\n\n')
	}
	
	start <- proc.time()
	results <- suppressWarnings(execQuery(queries[q], connection=conn))
	time <- proc.time() - start
	df <- data.frame(Variable=character(), Type=character(), Missing=numeric(), Levels=character(), stringsAsFactors=FALSE)

	for(i in 1:ncol(results)) {
		v = names(results)[i]
		if(class(results[,i])[1] == 'factor') {
			t = paste('Factor with ', length(levels(results[,i])), ' levels', sep='')
			if(length(levels(results[,i])) > nlevels) {
				l = paste(levels(results[,i])[1:nlevels], collapse='; ')
			} else {
				l = paste(levels(results[,i]), collapse='; ')		
			}
		} else {
			t = paste(class(results[,i]), collapse=", ")
			l = ''
		}
		m = length(which(is.na(results[,i]))) / nrow(results) * 100
		df = rbind(df, data.frame(Variable=v, Type=t, Missing=m, Levels=l))
	}
	
	cat('\\subsection{Results}\n')
	cat(paste('Returned ', nrow(results), ' rows and ', ncol(results), ' columns. Took ', format(time[1], digits=1), ' seconds to execute query.\n\n', sep=''))
	#x = xtable(df, caption=NULL, label=paste(queries[q], '-results', sep=''),
	#		   align=c('l','l','r','r','p{3.0in}'), digits=2)
	#print(x, include.rownames=FALSE)
	cat('\\begin{description}\n')
	for(r in 1:nrow(df)) {
		cat(paste('\n\n\\item[', sanitizeLatex(df[r,]$Variable), '] ', 
				  df[r,]$Type, 
				  ' (', format(df[r,]$Missing, digits=2), '\\% missing)\n',
				  sanitizeLatex(df[r,]$Levels),
			sep=''))
	}
	cat('\n\n\\end{description}\n\n')

	cat('\\subsection{SQL}\n\\begin{lstlisting}\n')
	cat(sql)
	cat('\n\\end{lstlisting}\n')
	
	}, error=function(e) { print(e) } )
}




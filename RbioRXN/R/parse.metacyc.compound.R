parse.metacyc.compound <-
function(datPath) {
  
  dat = readLines(datPath)
  
  ## remove HTML code
  
  htmlFrom = c("&alpha;", "&beta;", "&delta;", "&Delta;", "&gamma;", "&omega;", "&mdash;", "&Psi;", "&psi;", "&chi;", "&xi;", "&zeta;", "&harr;", "&pi;", "&tau;", "&lambda;", "&amp;", "&kappa;", "&iota;", "<sup>", "</sup>", "<sub>", "</sub>", "</sub>", "<SUP>", "</SUP>", "<SUB>", "</SUB>", "<i>", "<I>", "</I>", "</i>","<em>", "</em>", "<small>", "</small>", "&larr;", "&rarr;", "&epsilon;", '&prime;', '<b>', '</b>')
  
  htmlTo = c("alpha", "beta", "delta", "Delta", "gamma", "omega", "-", "Psi", "psi", "chi", "xi", "zeta", "<->", "pi", "tau", "lambda", "&", "kappa", "iota", "", "", "", "", "", "", "",  "", "",  "", "", "","","","","","","->","<-","epsilon", "'", '', '')
  
  for(i in 1:length(htmlFrom)) {
    dat = gsub(htmlFrom[i], htmlTo[i], dat)
  }
  
  cat('processing dat\n')
  # Process continuing lines
  con_close = grepl('^//$', dat)
  #con_not_close = grepl('^//$', dat) == FALSE
  ind_close = grep('^//$', dat)
  ind_not_close = (1:length(dat))[-ind_close]
  ind_comment = grep('^/', dat)
  ind_continue = intersect(ind_not_close, ind_comment)
  con_continue = logical(length(dat))
  con_continue[ind_continue] = TRUE
  
  dat2 = character(length(dat) - length(ind_continue))
  j=1
  for(i in 1:length(dat)) {
    if(con_continue[i]) {
      comment = sub('^/', '', dat[i])
      dat2[j-1] = paste(dat2[j-1], comment)
    } else {
      dat2[j] = dat[i]
      j = j+1
    }
  }
  
  # Remove weird characters
  dat2 = gsub('\U3e32393c', "'", dat2)
  dat2 = gsub('\U3e31393c', "'", dat2)
  dat2 = gsub('\U3e36393c', '-', dat2)
  dat2 = gsub('\U3e63663c', 'u', dat2)
  dat2 = gsub('\U3e37653c', 'c', dat2)
  dat2 = gsub('\U3e33393c', '"', dat2)
  dat2 = gsub('\U3e34393c', '"', dat2)
  dat2 = gsub('\U3e36663c', 'o', dat2)
  dat2 = gsub('\U3e39653c', 'e', dat2)
  dat2 = gsub('\U3e30623c', '&deg;', dat2)
  dat2 = gsub('\t', '', dat2)
  
  # Split with ' - '
  regexp = '([[:graph:]]+)( - )([[:print:]]+)'
  db = sub(regexp, '\\1', dat2[grep(regexp, dat2)])
  id = sub(regexp, '\\3', dat2[grep(regexp, dat2)])
  fields = unique(db)
  
  # Indexing start and end
  con_uniqueId = grepl('UNIQUE-ID', db)
  
  # integrate same field
  tmp_id = numeric(length(db))
  j = 0
  for(i in 1:length((tmp_id))) {
    if(con_uniqueId[i]) {
      j = j+1
    }
    tmp_id[i] = j
  }
  
  cat('Processing hash table\n')
	table = data.table(data.frame(cbind(tmp_id, db, id), stringsAsFactors=FALSE))
	table$tmp_id = as.numeric(table$tmp_id)
  table2 = table[,paste(id, collapse='///'),by="tmp_id,db"]
	table2$db = gsub('-', '.', table2$db)
	setkey(table2, tmp_id)
  
  # Create empty data frame
  metacyc = data.frame(matrix(ncol=length(fields), nrow=length(which(con_uniqueId))))
  colnames(metacyc) = gsub('-', '.', fields)
  metacyc = data.table(metacyc, tmp_id = 1:nrow(metacyc))
  metacyc = metacyc[,lapply(.SD, as.character)]
  metacyc$tmp_id = as.numeric(metacyc$tmp_id)
  setkey(metacyc,tmp_id)
    
  cat("parsing", length(dat2), 'lines\n')

  for(i in 1:nrow(metacyc)) {
    
    if(i == floor(nrow(metacyc)/10)) {
      cat('10% finished\n')
    } else if(i == floor(nrow(metacyc)/5)) {
      cat('20% finished\n')
    } else if(i == floor(nrow(metacyc)/2)) {
      cat('50% finished\n')
    }

    db = table2[list(i)][['db']]
    id = table2[list(i)][['V1']]

    metacyc[i, (db):= as.list(id),with=F]
  }

  # Post-process
  for(i in names(metacyc)) {
    metacyc[[i]] = sub('^///', '', metacyc[[i]])
  }
  metacyc[is.na(metacyc)] = ''
  
  # Post-process DBLINK
  metacyc2 = data.frame(metacyc, stringsAsFactors=FALSE)
  metacyc_dblink = build.subtable(metacyc2, 'UNIQUE.ID', 'DBLINKS', '///')
  regexp = '(\\()(.+)( ")(.+)(".*\\))'
  db = sub(regexp, '\\2', metacyc_dblink$DBLINKS)
  id = sub(regexp, '\\4', metacyc_dblink$DBLINKS)
	metacyc_dblink = data.table(metacyc_dblink, db, id)
	metacyc_dblink = metacyc_dblink[,paste(id, collapse='///'),by="UNIQUE.ID,db"]
	setnames(metacyc_dblink, 'V1', 'id')
 
  UNIQUE.ID = '' # To bypass 'R CMD check' no binding global variable
	setkey(metacyc_dblink, UNIQUE.ID)
  setkey(metacyc, UNIQUE.ID)
 
  cat('processing DBLINKS\n')

	for(i in metacyc$UNIQUE.ID) {
		db_list = unlist(metacyc_dblink[i, 'db', with=F])
		id_list = unlist(metacyc_dblink[i, 'id', with=F])
		metacyc[i, db_list := as.list(id_list), with=F]
	}

  metacyc[is.na(metacyc)] = ''
  if(length(which(is.na(names(metacyc)))) > 0) {
    metacyc[[which(is.na(names(metacyc)))]] = NULL
  }
	

	# Edit chemical formula
	metacyc$CHEMICAL.FORMULA = gsub('\\)///\\(', '', metacyc$CHEMICAL.FORMULA)
	metacyc$CHEMICAL.FORMULA = gsub('[\\( \\)]', '', metacyc$CHEMICAL.FORMULA)

	## Two letters chemical --> both capital (e.g., BR instead of Br)
	i=1
	while(i <= dim(metacyc)[1]) {
  	second_char = str_locate(metacyc$CHEMICAL.FORMULA[i], '[A-Z]{2}')[2]
  	if(is.na(second_char) == F) {
			substring(metacyc$CHEMICAL.FORMULA[i], second_char, second_char) = tolower(substring(metacyc$CHEMICAL.FORMULA[i], second_char, second_char))
  	}
  	if(grepl('[A-Z]{2}', metacyc$CHEMICAL.FORMULA[i]) == F) {
    	i = i+ 1
  	}
	}
  
  metacyc3 = data.frame(metacyc, stringsAsFactors=FALSE)
  metacyc3$tmp_id = NULL
  metacyc3$DBLINKS = NULL
  cat('Parsed compounds:', nrow(metacyc3), '\n')
  return(metacyc3)
}

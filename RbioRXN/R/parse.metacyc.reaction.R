parse.metacyc.reaction <-
function(datPath) {

  dat = readLines(datPath)
	
  ## remove HTML code

	htmlFrom = c("&alpha;", "&beta;", "&delta;", "&Delta;", "&gamma;", "&omega;", "&mdash;", "&Psi;", "&psi;", "&chi;", "&xi;", "&zeta;", "&harr;", "&pi;", "&tau;", "&lambda;", "&amp;", "&kappa;", "&iota;", "<sup>", "</sup>", "<sub>", "</sub>", "</sub>", "<SUP>", "</SUP>", "<SUB>", "</SUB>", "<i>", "<I>", "</I>", "</i>","<em>", "</em>", "<small>", "</small>", "&larr;", "&rarr;", "&epsilon;", '&prime;')

	htmlTo = c("alpha", "beta", "delta", "Delta", "gamma", "omega", "-", "Psi", "psi", "chi", "xi", "zeta", "<->", "pi", "tau", "lambda", "&", "kappa", "iota", "", "", "", "", "", "", "",  "", "",  "", "", "","","","","","->","epsilon", "'")

	for(i in 1:length(htmlFrom)) {
		dat = gsub(htmlFrom[i], htmlTo[i], dat)
	}
  
  cat('processing dat\n')
  # Process continuing lines
  con_close = grepl('^//$', dat)
  con_not_close = grepl('^//$', dat) == FALSE
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
    
  # Split with ' - '
  dat_split = strsplit(dat2, ' - ')
  
  # Indexing start and end
  con_close = grepl('^//$', dat2)
  con_not_close = grepl('^//$', dat2) == FALSE
  start = max(grep('^#$', dat2)) + 1
  con_skip = grepl('^\\^', dat2) == FALSE
	ind_coefficient = grep('\\^COEFFICIENT', dat2)
	ind_compartment = grep('\\^COMPARTMENT', dat2)
	ind_right = grep('^RIGHT', dat2)
	ind_left = grep('^LEFT', dat2)
	ind_participant = c(ind_right, ind_left)
	con_participant = logical(length(dat2))
	con_participant[ind_participant] = TRUE
  
  
  # Process coefficient, compartment
  cat('processing coefficient, compartment\n')
  for(i in ind_coefficient) {
    coefficient = dat_split[[i]][2]
    dat_split[[i-1]][2] = paste(coefficient, dat_split[[i-1]][2])
  }
  
  for(i in ind_compartment) {
    compartment = dat_split[[i]][2]
    if(con_participant[i-1]) {
      dat_split[[i-1]][2] = paste(dat_split[[i-1]][2], '(', compartment, ')', sep='')
    } else {
      dat_split[[i-1]][2] = paste(dat_split[[i-2]][2], '(', compartment, ')', sep='')
    }
  }
     
  metacyc = data.frame()
  metacycRow = list()
  
  cat("parsing", length(dat_split), 'lines\n')
  for(i in start:length(dat2)) {
        
    if(i == floor(length(dat_split)/10)) {
      cat('10% finished\n')
    } else if(i == floor(length(dat_split)/5)) {
      cat('20% finished\n')
    } else if(i == floor(length(dat_split)/2)) {
      cat('50% finished\n')
    }
    
    if(con_not_close[i] && con_skip[i]) {
      field = dat_split[[i]][1]
      value = dat_split[[i]][2]
      
      metacycRow[[field]] = paste(metacycRow[[field]], value, sep='///')
    } else if(con_close[i]) {
      metacycRow = as.data.frame(metacycRow, stringsAsFactors=FALSE)
      metacyc = rbind.fill(metacyc, metacycRow)
      metacycRow = list()
    }
  }
    
  # Post-process
  for(i in names(metacyc)) {
    metacyc[[i]] = sub('^///', '', metacyc[[i]])
  }
  metacyc[is.na(metacyc)] = ''
   
  # Assemble reaction
  left = gsub('///', ' + ', metacyc$LEFT)
  right = gsub('///', ' + ', metacyc$RIGHT)
  ind_reverse = grep('RIGHT-TO-LEFT', metacyc$REACTION.DIRECTION)
  
  direction = sub('RIGHT-TO-LEFT', ' => ', metacyc$REACTION.DIRECTION)
  direction = sub('PHYSIOL-', '', direction)
  direction = sub('IRREVERSIBLE-', '', direction)
  direction = sub('LEFT-TO-RIGHT', ' => ', direction)
  direction = sub('REVERSIBLE', ' <=> ', direction)
  direction = sub('^$', ' <=> ', direction)
  
  reaction = matrix("", nrow(metacyc), 3)
  reaction[ind_reverse,] = cbind(right[ind_reverse], direction[ind_reverse], left[ind_reverse])
  
  ind_normal = (1:nrow(metacyc))[-ind_reverse]
  reaction[ind_normal,] = cbind(left[ind_normal], direction[ind_normal], right[ind_normal])

  reaction_assembled = apply(reaction, 1, paste, collapse='')
  
  reaction_assembled = gsub('CCO-IN', 'in', reaction_assembled)
  reaction_assembled = gsub('CCO-OUT', 'out', reaction_assembled)
  
  metacyc[['equation']] = reaction_assembled
  metacyc$NA. = NULL
  
  # Post-process DBLINK
  if('DBLINKS' %in% names(metacyc)) {
    metacyc_dblink = build.subtable(metacyc, 'UNIQUE.ID', 'DBLINKS', '///')
    regexp = '(\\()(.+)( ")(.+)(".*\\))'
    db = sub(regexp, '\\2', metacyc_dblink$DBLINKS)
    id = sub(regexp, '\\4', metacyc_dblink$DBLINKS)
  
    tmp_dblink = matrix("", nrow=nrow(metacyc), ncol=length(unique(db)))
    colnames(tmp_dblink) = unique(db)
    db = gsub('-', '.', db)
    tmp_dblink = data.frame(tmp_dblink, stringsAsFactors=F)
  
    cat('processing DBLINKS\n')
    for(i in 1:nrow(metacyc)) {
      index = metacyc_dblink$UNIQUE.ID == metacyc$UNIQUE.ID[i]
      db_index = db[index]
      id_index = id[index]
      if(length(index) > 0) {
        duplicated_entry = names(table(db[index])[table(db[index]) > 1])
        if(length(duplicated_entry) > 0) {
          for(j in duplicated_entry) {
            tmp_dblink[i, j] = paste(id[index][db[index] == j], collapse='///')
          }
          db_index = db[index][db[index] %in% duplicated_entry == FALSE]
          id_index = id[index][db[index] %in% duplicated_entry == FALSE]
        }
        tmp_dblink[i, db_index] = id_index
      }
    }
  
    tmp_dblink[is.na(tmp_dblink)] = ''
    metacyc2 = cbind(metacyc, tmp_dblink)
    metacyc2$DBLINK = NULL  
    metacyc2[is.na(metacyc2)] = ''
  } else {
    return(metacyc)
  }
  cat('Parsed reactions:', nrow(metacyc2), '\n')
  return(metacyc2)
}

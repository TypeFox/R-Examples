`createUniqueLongData` <- 
function(long.data) {

	YEAR <- ID <- NULL

	tmp.key <- key(long.data)
	tmp.last.year <- tail(sort(unique(long.data$YEAR)), 1)
	tmp.dups.index <- data.table(unique(long.data[duplicated(long.data)][, tmp.key, with=FALSE])[long.data, nomatch=0][,setdiff(tmp.key, "YEAR"), with=FALSE], key=setdiff(tmp.key, "YEAR"))
	setkeyv(long.data, setdiff(tmp.key, "YEAR"))
	tmp.unique.data <- long.data[!tmp.dups.index]
	tmp.past.dups.extended <- long.data[YEAR!=tmp.last.year][tmp.dups.index, allow.cartesian=TRUE, nomatch=0]
	tmp.current.dups <- data.table(long.data[unique(tmp.dups.index)][YEAR==tmp.last.year], key=setdiff(tmp.key, "YEAR"))
	tmp.all.dups.extended <- data.table(rbindlist(list(tmp.past.dups.extended, tmp.current.dups)), key=tmp.key)
	tmp.all.dups.extended[,ID:=paste(ID, "DUPS", tmp.all.dups.extended[,seq.int(.N), by=eval(tmp.key)][['V1']], sep="_")]
	return(data.table(rbindlist(list(tmp.unique.data, tmp.all.dups.extended)), key=tmp.key))
} ### END createUniqueLongData

dataIDs <- read.csv('dataIDs.csv')
dataNames <- read.csv('dataNames.csv')
dataFiles <- read.csv('dataFiles.csv')
levels(dataIDs$type) <- c('Example','Problem')
dataInfo <- cbind(dataNames,dataIDs)
save(dataInfo, file='../pkg/data/dataInfo.rda')


findData <- function(x, chapters=1:21, types=c('Example','Problem'), pattern='*',ignore.case=TRUE) {
	if (! missing(x) ) {
		if ( all( x %in% c('Example','Problem') ) ) { 
			types <- x 
		}
		if (is.character(x)) { pattern <- x }
		if (is.numeric(x)) { chapters <- x }
	}
	return(subset(dataInfo,
				chapter %in% chapters & 
				type %in% types &
				grepl(pattern,name,ignore.case=ignore.case) 
		)
	)
}

for (i in 1:nrow(dataInfo)) {
	read.csv( as.character(dataFiles$path)[i] ) -> foo
	name <- as.character(dataInfo$name[i])
	print(name)
	assign(noquote(name), foo)
}


Tobacco2 <- 
	data.frame(
		flower.length = c(
			rep(Tobacco$flower.length, times=Tobacco$f1.count),
			rep(Tobacco$flower.length, times=Tobacco$f2.count)),
		generation = rep(c('F1','F2'), times=c(sum(Tobacco$f1.count),sum(Tobacco$f2.count)))
	)

# Note: was error in original data (had bush winning popular vote against gore)
WillsDebates <- WillsDebates[,c(6,1,7,3,4,5)]
names(WillsDebates) <- c('year','winner','loser','winner.wills','loser.wills','diff.wills')

# add an extra column
SexualSelection$taxon.pair = toupper(letters[1:25])
Dioecy$taxon.pair <- 1:nrow(Dioecy)
BrookTrout$proportion.surviving <- with(BrookTrout, salmon.survived/salmon.released)

# some little tweaks

DayOfBirth$day <- with(DayOfBirth, ordered(day, levels=day))
Powerball$day <- with(Powerball, ordered(day, levels=day))


# save the data to .rda files
data_dir <- "data"
list <- c(as.character(dataInfo$name),'Tobacco2')
for (item in list) {
	try(save(list = item, file = file.path(data_dir, sprintf("%s.rda", item))))
}

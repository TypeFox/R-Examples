'gluesheets' <- function(PATH = NULL, out.file = 'fileout.txt'){

if(is.null(PATH))
  PATH <- paste(getwd(), '/', sep = '')

file.names <- dir(PATH)
n <- length(file.names)

allspp <- ''
sheet <- vector(mode = 'list', length = n)
for(i in 1:n){
  this <- read.table(paste(PATH, file.names[i], sep = ''),                      			header = TRUE, sep = '\t')
  cat('reading sheet', file.names[i], '\n')
  uniquifier <- duplicated(this[,1])
  uniquifier <- as.character(as.numeric(uniquifier))
  uniquifier[uniquifier == 0] <- ''
  rownames(this) <- paste(this[,1], uniquifier, sep = '')
  allspp <- c(allspp, rownames(this))
  this <- t(this[,3:(ncol(this)-2)])
  rownames(this) <- gsub('[ .]', '', rownames(this), perl = TRUE)
  colnames(this) <- gsub('[ .]', '', colnames(this), perl = TRUE)
  this <- data.frame(cbind(depth = rownames(this), this))
  sheet[[i]] <- data.frame(cbind(sheet = rep(i, nrow(this)), this))
}
allspp <- unique(allspp)[-1]

table.out <- as.data.frame(sheet[[1]])
for(i in 2:n){
  table.out <- merge(table.out, as.data.frame(sheet[[i]]), all = TRUE)
}

table.out <- apply(table.out, 2, as.character)

cat('this collection of data contains counts of', ncol(table.out), 'unique species at', nrow(table.out), 'depths\n')
cat(sum(is.na(table.out)), 'missing values replaced with zeros\n')

table.out[is.na(table.out)] <- 0

no <- as.numeric(table.out[,1]) - 1
na <- file.names[as.numeric(table.out[,1])]
na <- substr(na, 1, nchar(na) - 4)
dp <- table.out[,2]

rn <- paste('', no, na, dp, sep = '_')

table.out <- table.out[,-(1:2)]
table.out <- as.data.frame(table.out)

table.out <- apply(table.out, 2, as.numeric)
rownames(table.out) <- rn

sink(file = out.file)
cat('name ')
write.table(table.out, quote = FALSE)
sink()

}
"readGPDascii" <-
function(x){

everything <- scan(x, sep = '\n', what = '')
everything <- iconv(everything, "latin1", "ASCII", "?")
metadata <- gsub('^[^#]*', '', everything, perl = TRUE, useBytes = TRUE)
metadata <- metadata[metadata != '']
metadata <- substr(metadata, 2, nchar(metadata))
justdata <- gsub('#.*$', '', everything, perl = TRUE, useBytes = TRUE)
justdata <- justdata[justdata != '']

taxa <- as.numeric(unlist(strsplit(justdata[1], split = ' '))[1])
levels <- as.numeric(unlist(strsplit(justdata[1], split = ' '))[2])

cat('Number of taxa: ', taxa, '\n')
cat('Number of levels: ', levels, '\n')

short.tax.names <- substr(justdata[2:(taxa + 1)], 7, 14)
tax.cat <- as.factor(substr(justdata[2:(taxa + 1)], 16, 16))
tax.names <- substr(justdata[2:(taxa + 1)], 18,
                    max(nchar(justdata, type = 'width')))

justdata <- justdata[(taxa+2):length(justdata)]
data.labels <- justdata[grep(',', justdata, perl = TRUE)]
data.vector <- justdata[-grep(',', justdata, perl = TRUE)]

data.labels <- strsplit(data.labels, split = ',')
data.labels <- matrix(unlist(unlist(data.labels)), ncol = 3,
                      byrow = TRUE)
depths <- as.numeric(data.labels[,1])
sample.names <- data.labels[,2]
absolute.ages <- as.numeric(data.labels[,3])

count.vector <- paste(data.vector, collapse = ' ')
count.vector <- strsplit(count.vector, split = '[[:space:]]+',
                         perl = TRUE)[[1]]
count.matrix <- matrix(as.numeric(as.character(count.vector[-1])),
                       ncol = taxa, byrow = TRUE)
colnames(count.matrix) <- short.tax.names
rownames(count.matrix) <- depths

strat.col.out <- list(counts = count.matrix, depths = depths,
            sample.names = sample.names, absolute.ages = absolute.ages,
            taxa = tax.names, short.names = short.tax.names,
            tax.cat = tax.cat, metadata = metadata)
class(strat.col.out) <- 'strat.column'

return(strat.col.out)
} # End of function


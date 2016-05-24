## ----, echo=FALSE--------------------------------------------------------
knitr::knit_hooks$set(small.margin = function(before, options, envir) {
    if (before) par(mar = c(2,3,3,2))
})

## ----drug-type-----------------------------------------------------------
drug_class = 'PI' # Possible drug types are 'PI', 'NRTI', and 'NNRTI'. 

## ----raw-data------------------------------------------------------------
base_url = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')

gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')

## ----, results='asis', echo=FALSE----------------------------------------
knitr::kable(head(gene_df[,4:20]))
knitr::kable(head(tsm_df))

## ----cleaned-data--------------------------------------------------------
# Returns rows for which every column matches the given regular expression.
grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df))
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]

## ----design-matrix-------------------------------------------------------
# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[,4:(pos_start-1)]

## ----, results='asis', echo=FALSE----------------------------------------
knitr::kable(X[1:10,1:10])

## ----, results='asis', echo=FALSE----------------------------------------
knitr::kable(head(Y))

## ----small.margin=TRUE---------------------------------------------------
hist(Y[,1], breaks='FD')

## ----small.margin=TRUE---------------------------------------------------
hist(log(Y[,1]), breaks='FD')

## ----knockoff------------------------------------------------------------
library(knockoff)

knockoff_and_bhq <- function (X, y, fdr) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
    
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # Run the knockoff filter.
  result = knockoff.filter(X, y, fdr = fdr, knockoffs = 'equicorrelated')
  knockoff_selected = names(result$selected)
  
  # Run BHq.
  p = ncol(X)
  lm.fit = lm(y ~ X - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  cutoff = max(c(0, which(sort(p.values) <= fdr * (1:p) / p)))
  bhq_selected = names(which(p.values <= fdr * cutoff / p))
  
  list(Knockoff = knockoff_selected, BHq = bhq_selected)
}

fdr = 0.20
results = lapply(Y, function(y) knockoff_and_bhq(X, y, fdr))

## ------------------------------------------------------------------------
print(results[1])

## ----comparisons---------------------------------------------------------
get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

comparisons <- lapply(results, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})

## ------------------------------------------------------------------------
print(comparisons[1])

## ----small.margin=TRUE---------------------------------------------------
for (drug in names(comparisons)) {
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}


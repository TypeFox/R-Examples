### Nonspecific Filtering via Variance
# Input Matrix: cols=Samples, rows=Protein IDs
# Filtering via:  - IQR
#                 - overall Variance
#                 - noVar: no variance filter, only 1.filter is applied
# Cutoff Selection via: - shorth (shortest intervall containing 50% of the data) , default: NA
#                       - define quantile cutoff, e.g. 0.5 -> 50% data filtered of smallest Var or IQR
# limit = number of expected candidates in the data

#### Functions
#require(genefilter)

my.rowIQR <- function(mat) {          #adapted for matrix from genefilter:::rowIQRs 
rowQ(mat, ceiling(0.75 * ncol(mat))) - rowQ(mat, floor(0.25 * ncol(mat)))
}

varFilter <- function (mat, baittab, func= c("IQR", "overallVar", "noVar"), var.cutoff = NA, limit=0)
{
  func <- match.arg(func)
  
  Cpos <- match(baittab$V1[grep("C",baittab$V3)], colnames(mat) )     # columns in mat for "C"controls
  Tpos <- match(baittab$V1[grep("T",baittab$V3)], colnames(mat) )     # columns in mat for "T"samples
  ctrl.med <- apply(mat[,Cpos], 1, function(x){median(x)} )   # 1.filter: ctrl counts > bait counts
  bait.med <- apply(mat[,Tpos], 1, function(x){median(x)} )
  mat <- mat[-which(ctrl.med > bait.med), ]
cat("Biological filter reduced the number of proteins to: ", dim(mat)[1], "\n")
if (dim(mat)[1] < limit) {func <- "noVar"; cat("Argument limit is chosen too small, in order to hold it the variance filter cannot be applied\n") } 

  if (func == "IQR") {                                       #2. filter: small Var over all samples
        vars <- my.rowIQR(mat)     }
  else if (func == "overallVar")  {
       vars <- rowVars(mat)        }
  else if (func == "noVar")  {
       return(mat)        }
  else stop("Define filter function!\n")


    if (is.na(var.cutoff))  {                       # calculate shorth for cutoff
     var.cut <- shorth(vars, tie.action="min")
     if(var.cut <= median(vars))  selected <- which(vars > var.cut)   
     else stop("shorth calculation is not appropriate in this case, define a quantile for cutoff calculation!\n")  
     }
     
    else if (0 < var.cutoff && var.cutoff < 1)  {
     quant <- quantile(vars, probs = var.cutoff)   # cutoff corresponding to defined quantile
     selected <- which(vars > quant)        }
    else stop("Cutoff Quantile has to be between 0 and 1, alternative 'NA' for shortest interval calculation \n")

    mat.f <- mat[selected,]
    
    if (dim(mat.f)[1] > limit )  {
    cat("Variance filter reduced the number of proteins to: ", dim(mat.f)[1], "\n")
    return(mat.f)                 }
    else {
    cat("Chosen parameter setup filtered data below the threshold of the expected number of candidates \n")
    if (is.na(var.cutoff)) var.cutoff=0.4 
                                                # filter parameter setup failed to hold the limit, cutoff adjustment 
    while( dim(mat.f)[1] < limit)        {
    if (var.cutoff<=0.05) {cat("Decrease limit, variance filter cannot be applied \n"); return(mat) }
    else var.cutoff <- var.cutoff-0.05
    quant <- quantile(vars, probs = var.cutoff)   
    selected <- which(vars > quant) 
    mat.f <- mat[selected,]              
    }
    
    cat("Variance filter reduced the number of proteins to: ", dim(mat.f)[1], "\n")
    return(mat.f)
    }

}





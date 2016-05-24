#' Compare alternative multiple sequence alignments
#'
#' @param ref   The reference MSA (in fasta format)
#' @param com   The MSA to compare (in fasta format)
#'
#' @return Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal pairwise column alignment of two alternative MSAs of the same sequences, and summary statistics of the differences between them. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {reference_P}          {The numbered character matrix of the reference alignment}
#'  \item {comparison_Q}         {The numbered character matrix of the comparison alignment}
#'  \item {results_R}            {The results summary matrix (containing column averages of match, gapcon, merge, split, shift)}
#'  \item {similarity_S}         {The similarity matrix between the reference and comparison alignment columns}
#'  \item {dissimilarity_D}      {The dissimilarity matrix between the reference and comparison (containing match, gapcon, merge, split, shift)}
#'  \item {dissimilarity_simple} {The dissimilarity matrix with categories stacked into a single 2D matrix}
#'  \item {columnmatch}          {The column of the comparison alignment with the highest final match score}
#'  \item {cys}                  {The proportion of cysteines (relevant for cysteine rich proteins)}
#'  \item {reflen}               {The length of the reference alignment}
#'  \item {comlen}               {The length of the comparison alignment}
#'  \item {refcon}               {The consensus sequence of the reference alignment}
#'  \item {comcon}               {The consensus sequence of the comparison alignment}
#'  \item {score}                {The overall similarity score}
#' } 
#' 
#' @export
#' @examples
#' data(reference_alignment)
#' data(comparison_alignment)
#' PAC <- compare_alignments(reference_alignment,comparison_alignment)
#'
#' @note The `compare_alignments` compares two alternative multiple sequence alignments (MSAs) of the same sequences. The alternative alignments must contain the same sequences in the same order. The function classifies similarities and differences between the two MSAs. It produces the "pairwise alignment comparison" object required as the first step any other package functions. The function converts the MSAs into matrices of sequence characters labelled by their occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns have the highest similarty between the reference and comparison MSAs to generate a similarity matrix (excluding conserved gaps). From this matrix, the comparison alignment column with the similarity to each reference alignment column is used to calculate further statistics for dissimilarity matrix, summarised for each reference MSA column in the results matrix. Lastly, it calculates the overall similarity score between the two MSAs.
#'
compare_alignments <- function(ref,com){
  
  if (!is.data.frame(ref)){
    data.frame(seqinr::read.fasta(ref,set.attributes=FALSE)) -> ref
  }
  if (!is.data.frame(com)){
    data.frame(seqinr::read.fasta(com,set.attributes=FALSE)) -> com
  }
  if( !valid_alignments(ref,com)){
    stop("both alignments must contain the same sets of sequences in the same order")
  }
  
  
  ###########################################
  # Replacing letters with letter+occurance #
  ###########################################
  ref2  <- prepare_alignment_matrix(ref)
  com2  <- prepare_alignment_matrix(com)
  
  if (!is.data.frame(ref)){
    names <- row.names(as.matrix(seqinr::read.fasta(ref)))
  } else {
    names <- colnames(ref)
  }
  # Replacing "-" with NA in the test alignment means
  # that gaps don't count towards column matching score
  com[com=="-"]  <-NA
  com2[com2=="-"]<-NA
  
  ##################################
  # Alignment identity calculation #
  ##################################
  
  # res_list contains $results, $cat, $means 
  res_list = rcpp_align(ref2,com2)
  results  = res_list$results
  cat      = res_list$cat
  means    = res_list$means
  catnames = c("Match","Gapcon","Merge","Split","Shift")
  
  row.names(results)<-c("ColumnMatch",  # 1
                        "NonGap",       # 2
                        "Cys",          # 3
                        "RawMatch",     # 4
                        "Gapcon",       # 5
                        "Merge",        # 6
                        "Split",        # 7
                        "Shift",        # 8
                        "Match")        # 9

  # Format P and Q matrices for output
  ref3 <- t(ref2)
  com3 <- t(com2)
  rownames(ref3) <- names
  rownames(com3) <- names

  # Create dissimilarity (matrix D) from simplified dissimilarity (res_list$cat)
  dissimilarity_D <- array(dim      = c(ncol(ref), # rows
                                        nrow(ref), # columns
                                        5),        # stacks
                           dimnames = list(names,
                                           NULL,
                                           catnames))
                                           
  dissimilarity_D[,,1] <- 1*t(cat=="M") # "Match"
  dissimilarity_D[,,2] <- 1*t(cat=="g") # "Gapcon"
  dissimilarity_D[,,3] <- 1*t(cat=="m") # "Merge"
  dissimilarity_D[,,4] <- 1*t(cat=="s") # "Split"
  dissimilarity_D[,,5] <- 1*t(cat=="x") # "Shift"

  # Write category averages to results (R matrix)
  results_R           <- t(colMeans(dissimilarity_D))
  rownames(results_R) <- catnames

  # For each column of ref, which column of com is most similar
  columnmatch <- as.vector(res_list$results[1,]) 
 
  # Ref cysteine occurance
  cys <- as.vector((t(rowMeans(ref=="c"))))
 
  # Count alignment columns
  reflen <- nrow(ref)                              
  comlen <- nrow(com)
  
  # Alignment consensus sequences
  refcon <- seqinr::consensus(t(ref))
  comcon <- seqinr::consensus(t(com))

  # Final mean score
  score <- mean(cat=="M")/(1-mean(cat=="g"))
  
  # Create final object
  list(reference_P          = ref3,
       comparison_Q         = com3,
       results_R            = results_R,
       similarity_S         = means,
       dissimilarity_D      = dissimilarity_D,
       dissimilarity_simple = t(cat),
       columnmatch          = columnmatch,
       cys                  = cys,
       reflen               = reflen,
       comlen               = comlen,
       refcon               = refcon,
       comcon               = comcon,
       score                = score)
}


prepare_alignment_matrix <- function(commat){
  mat2 <- rcpp_prepare_alignment_matrix(as.matrix(commat))
  # Remove extra space and de-number gaps
  gsub(x = mat2, pattern = " ",     replacement = "")  -> mat2
  gsub(x = mat2, pattern = "[-].*", replacement = "-") -> mat2
  mat2
}


valid_alignments <- function(ref,com){
  checks = sapply(1:ncol(ref),function(i){
    r=as.character(ref[,i])
    a=as.character(com[,i])
    dg_ref = r[r!="-"]
    dg_com = a[a!="-"]
    all(dg_com==dg_ref)
  })
  all(checks)
}

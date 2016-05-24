score.model <- function(M, saved.scores, cross, addcov, intcov, threshold,
                        verbose = TRUE, ...)
{
  n.pheno <- ncol(M)
  model.score <- 0
  mymodel <- rep(NA,n.pheno)
  count.score <- 0
  update.scores <- list(code = NULL, pheno.col = NULL, bic = NULL)
  
  for(i in 1:n.pheno){
    pheno <- node.parents(M, i)
    mymodel[i] <- paste("(",paste(pheno$identifier,")",sep=""),sep="")

    ## Find saved BIC score if already computed.
    bic <- find.bic(pheno$code, i, update.scores, saved.scores, ...)

    if(is.na(bic)){
      run <- calc.bic(cross, pheno$code, i, pheno$parents,
                      addcov, intcov, threshold,
                      n.pheno, update.scores = update.scores,
                      saved.scores = saved.scores, ...)
      
      ## Update saved scores.
      update.scores$code <- c(update.scores$code, as.character(run$code))
      update.scores$pheno.col <- c(update.scores$pheno.col, run$pheno.col)
      update.scores$bic <- c(update.scores$bic, run$bic)

      ## Print scan info if verbose.
      count.score <- count.score + 1
      if(verbose) {
        if(count.score == 1) cat("\nscan ")
        cat(paste("(", paste(run[, "pheno.col"], collapse = ","), "|",
                  paste(pheno$parents, collapse = ","), ")", sep = ""), "")
      }
      bic <- run[i == run[, "pheno.col"], "bic"]
    }
    ## Accumulate model score.
    model.score <- model.score + bic
  }
    
  list(model.score = model.score,
       update.scores = update.scores,
       model.name = paste(mymodel,collapse=""))
}
###########################################################################################
agree.covs <- function(x,y) {
  out <- length(x) == length(y)
  if(out & length(x))
    out <- all(sort(x) == sort(y))
  out
}
###########################################################################################
find.bic <- function(code, pheno.col, update.scores = NULL, saved.scores = NULL, ...)
{
  bic <- rep(NA, length(code))
  
  if(!is.null(saved.scores)) {
    n.scores <- nrow(saved.scores)
    score.pointer <- match(code, dimnames(saved.scores)[[1]])
    bic <- saved.scores[score.pointer + n.scores * (pheno.col - 1)]
  }
  bic.na <- is.na(bic)
  if(any(bic.na) & !is.null(update.scores)) {
    wh <- match(paste(code[bic.na], pheno.col[bic.na], sep = "."),
                paste(update.scores$code, update.scores$pheno.col, sep = "."),
                nomatch = 0)
    if(!all(wh == 0))
      bic[bic.na][wh > 0] <- update.scores$bic[wh]
  }
  bic 
}
######################################################################
node.parents <- function(M, node)
{
  aux <- which(M[,node] == 1)
  if(length(aux) == 0){ 
    parents <- NULL
    identifier <- as.character(node)
  }
  else{
    parents <- aux
    identifier <- paste(node,paste(parents,collapse=","),sep="|")
  }
  code <- find.code(node, parents)
  list(parents=parents, identifier=identifier, code = code)
}

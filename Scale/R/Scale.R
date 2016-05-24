Scale <- function(data=data.frame(), orders=list(), 
                  orders_id=integer(), reverse=c(), 
                  items=character(), col_names=character()){
  # A class to represent scale data
  out <- list(data=data, orders=orders, 
              orders_id=orders_id, reverse=reverse, 
              items=items, col_names=col_names)
  class(out) <- "ScaleData"
  return(out)
}


PreProc <- function(sc){
  # operates on a ScaleData object
  if (class(sc)!= "ScaleData") stop("Input should be of class ScaleData.")
  if (!(any(class(sc$data) == "data.frame"))){ 
    stop("Data has to be a data.frame.\n
         See ?data.frame, ?cbind, ?matrix")
  }
  if (length(sc$orders) != 0){
    # This should handle various orderings of printed questionnaires
    if (length(sc$orders_id)==0) stop("You need to provide orders_id.\n
                                      This is an indicator of which order
                                      each participant received.")
    reordered <- lapply(seq_along(sc$orders), function(x){
      df <- sc$data[sc$orders_id==x,]
      df <- df[,order(sc$orders[[x]])]
      df
    })
    data <- data.frame(Reduce(rbind, reordered))
    
    
  } else{data <- sc$data}
  if (length(sc$reverse)!=0){
    data[,sc$reverse] <- (max(data, na.rm=T)+1) - data[,sc$reverse]
  } else {warning("Continue without reverse items.")}
  if (length(sc$col_names)!=0) colnames(data) <- sc$col_names
  if (length(sc$items)!=0){
      items <- sc$items
    }
  if (length(sc$items)==0){ 
    items <- NULL
    warning("Items' vector unspecified. out$items is NULL...")}
  return(list(data=data, items=items))
}

.item_set <- function(x, rcut=0.4){
  low_cor <- which(x$item.stats$r.cor < rcut) 
  a_drop <- which(x$alpha.drop$raw_alpha >= x$total$raw_alpha)
  list(low_cor=low_cor, a_drop=a_drop)
}

ItemAnalysis <- function(prep, method="spearman", fm="gls", nfactors=1, 
                         rcut=0.3, score_type = "z", exclude=c()){
  
  data <- if (length(exclude)!=0) prep$data[-exclude] else prep$data
  if (method=="polychoric") cor_matrix <- psych::polychoric(data)$rho
  if (method=="spearman") cor_matrix <- Hmisc::rcorr(as.matrix(data), 
                                                     type="spearman")$r
  al <- psych::alpha(cor_matrix)
  sug <- .item_set(al, rcut)
  rely <- list(alpha= al, k=ncol(data), 
               title=deparse(substitute(prep)), method=method, 
               suggest = sug)
  class(rely) <- "reliability"
  my_kmo <- psych::KMO(cor_matrix)
  my_bartlett <- psych::cortest.bartlett(cor_matrix, n=nrow(data))
  my_fa <- psych::fa(cor_matrix, fm=fm, nfactors=nfactors)
  scores <- apply(data, 1, function(x) sum(scale(x) * my_fa$loadings[,1]))
  if (score_type== "t") scores <- (scores*10)+50
  if (score_type == "sten") scores <- (scores*2)+5.5
  valid <- list(model=my_fa, method=fm,loadings=my_fa$Structure[,1], kmo=my_kmo, 
                bartlett=my_bartlett, scores=scores)
  class(valid) <- "validity"
  items <- prep$items
  out <- list(data=data, items=items, rely=rely, valid=valid)
  class(out) <- "ItemAnalysis"
  out
}

print.ItemAnalysis<- function(x,...){
  print(x$rely)
  print(x$valid)
}

print.validity <- function(x,...){
  mod <- x$model$fa
  cat("A", x$method, "factor analysis was conducted. Items were regressed to
      a single factor. Their loadings are the following:\n")
  print(sort(x$loadings))
}

print.reliability <- function(x,...){
  ## or Rmd file?
  cat("\nReliability Analysis of", x$title, "ScaleData object.",
      "\n\nA",x$method, "correlation matrix of", x$k,  
      "items was calculated and submitted to Reliability analysis.
      \nThe overall Cronbach's Alpha was", round(x$alpha$total$raw_alpha, 2), ".\n")
  if (length(x$suggest$low_cor)!=0){
    cat("Item(s) that exhibited low correlation with the rest of the scale were:\n",
        .comma_and(x$suggest$low_cor), ".\n")}
  if (length(x$suggest$a_drop)!=0){
    cat("Furthermore, deleting item(s)", .comma_and(x$suggest$a_drop),
        "may improve reliability.")}
}

.comma_and <- function(x){
  if (length(x)==1) return(x)
  paste(paste(x[1:(length(x)-1)], collapse=","), "and", x[length(x)])
}

ReportTable <- function(it, write_file=FALSE, sep=";"){
  item_means <- colMeans(it$data, na.rm=T)
  item_sds <- apply(it$data, 2, sd, na.rm=T)
  out1 <- data.frame(Item = colnames(it$data),
                     Corr = it$rely$alpha$item.stats[,"r.drop"],
                     Loadings = it$valid$loadings, 
                     Mean = item_means, 
                     SD = item_sds)
  out1 <- out1[order(out1$Loadings, decreasing=T),]
  colnames(out1) <- c("Item", "Corr. to scale", 
                      "Factor Loading", 
                      "Mean", "SD")
  if (write_file==TRUE){write.table(out1, "ItemAnalysisOutput.csv", 
                         sep=sep, row.names=F)}
  out1
}

GetScores <- function(it, write_file=FALSE, sep=";", 
                      scale_name="My_Scale"){
  my_scores <- data.frame(it$valid$scores)
  colnames(my_scores) <- scale_name
  if (write_file==TRUE){
    fname <- paste(scale_name, ".csv")
    write.table(my_scores, fname, sep=sep, row.names=F)}
  my_scores
}

ChooseBest <- function(it, n=5){
  loads <- it$valid$loadings
  if (n > length(loads)) stop("This n exceeds the number of available items!")
  names(head(sort(loads, T), n))
}

ShowItems <- function(it, n=5, 
                      write_file=FALSE, scale_name="MyItems"){
  if (is.null(it$items)) stop("You have not specified the items argument, 
                              when defining the ScaleData Object.
                              See ?ScaleData.")
  best <- ChooseBest(it, n)
  ind <- which(colnames(it$data) %in% best)
  out <- it$items[ind]
  if (write_file) {
    fn <- paste(scale_name, ".txt", collapse='')
    writeLines(out, fn)}
  list(best, out)
}

ReadItems <- function(filename, enc="UTF-8"){
  readLines(filename, encoding=enc)
}
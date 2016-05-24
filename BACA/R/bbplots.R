DAVIDsearch <- function(gene.lists,
                        david.user,
                        idType="AFFYMETRIX_3PRIME_IVT_ID", 
                        listType="Gene",
                        easeScore=1,
                        annotation="KEGG_PATHWAY",
                        species=NA)
{
  if(length(gene.lists) > 0 & !is.na(david.user)) {
    david.objs <- list()
    list.david.result <- list()
    for(i in 1:length(gene.lists)) {
      gene.list <- gene.lists[[i]]
      david.objs[[i]] <- DAVIDWebService$new(email=david.user)
      list.david.result[[i]] <- NA
      if(is.connected(david.objs[[i]])) {
        result <- addList(david.objs[[i]], gene.list, idType=idType, listName=names(gene.lists)[i], listType=listType)
        if(!is.na(species))  setCurrentSpecies(david.objs[[i]], species)
        cat("For the list ", names(gene.lists)[i], " you have:", "\n", sep="")
        cat(" - number of genes:",length(gene.list), "\n", sep="")
        cat(" - ", names(result)[1],":", result[[1]], "\n", sep="")
        cat(" - ", names(result)[2],":", length(result[[2]]), "\n", sep="")
        cat(" - species:", "\n")
        cat(getSpecieNames(david.objs[[i]]), "\n")
        setAnnotationCategories(david.objs[[i]], annotation)
        list.david.result[[i]] <- getFunctionalAnnotationChart(david.objs[[i]], threshold=easeScore, count=2)
        cat(" - number of ",annotation," terms found: ",dim(list.david.result[[i]])[1], "\n", sep="")
      }
      else
        warning(paste("Not able to connect with DAVID webservice by using:"), david.user)
    }
    return(list.david.result)
  }
  else stop("Provide a list of up-/down-regulate gene lists and a correct david user email.")
}


BBplot <- function(list.david.obj,
                   max.pval = 0.01,
                   min.ngenes = 5,
                   max.ngenes = 500,
                   adj.method = "Benjamini",
                   title = "BBplot",
                   name.com = "***",
                   labels = c("down", "up"),
                   colors = c("#009E73", "red"),
                   print.term = "full")
{
  if(length(list.david.obj) > 0 & all(unlist(lapply(list.david.obj, FUN = function(x) {inherits(x, "DAVIDFunctionalAnnotationChart")})))) {
    filt.dat <- lapply(list.david.obj, FUN = function(l) {
      if(adj.method == "") l <- l[which(l$PValue <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "Bonferroni") l <- l[which(l$Bonferroni <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "Benjamini") l <- l[which(l$Benjamini <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "FDR") l <- l[which(l$FDR <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      l <- l[which(l$PValue <= max.pval & l$Count >= min.ngenes),]
      if(print.term == "name") l$Term <- unlist(lapply(l$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
      else if(print.term == "description") l$Term <- unlist(lapply(l$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
      else if(print.term != "full") stop("Invalid input value for print.term.")
      l
    })
    list.paths <- unique(unlist(lapply(filt.dat, FUN=function(x) {unlist(x$Term)})))
    if(!is.null(list.paths)) {
      row <- rep(list.paths, length(filt.dat))
      col <- rep(1:length(filt.dat), each = length(list.paths))
      dn.up <- rep(rep(c("dn","up"), each = length(list.paths)),(length(filt.dat)/2))
      circle.size <- vector(mode = "numeric", length = (length(list.paths)*length(filt.dat)))
      significance <- vector(mode = "numeric", length = (length(list.paths)*length(filt.dat)))
      for(i in 1:length(filt.dat)) {
        sp <- ((i-1)*length(list.paths))
        for(j in 1:length(list.paths)){
          crt.path <- list.paths[j]
          ids <- which(filt.dat[[i]]$Term == crt.path)
          if(length(ids) > 0) {
            circle.size[(sp+j)] <- filt.dat[[i]]$Count[ids[1]]
            significance[(sp+j)] <- 1 - filt.dat[[i]]$PValue[ids[1]]
          }
          else {
            circle.size[(sp+j)] <- 0
            significance[(sp+j)] <- 0
          }
        }
      }
      bbplot <- NULL
      if((length(list.david.obj) %% 2) == 0 & length(colors) == 2) {
        grid.info <-  factor(rep(name.com, each = (length(list.paths)*2)), levels = name.com)
        dataset <- data.frame(row=row, col=col, circle.size=circle.size, significance=significance, dn.up=dn.up, exp=grid.info)
        bbplot <- ggplot(dataset, aes(y = factor(row),  x = factor(col))) +
          geom_point(data = subset(dataset, circle.size>0), aes(colour = dn.up, size = circle.size, alpha=significance)) +
          scale_colour_manual(breaks = c("dn", "up"), labels = labels, values = colors, name="Status gene") +
          scale_alpha(guide="none") +
          facet_grid(facets=.~exp, scales="free_x", space="free_x")  +
          scale_size(range = c(1, 10)) + 
          labs(title = title) +
          theme_bw() +
          theme(axis.text.x = element_blank(), plot.title = element_text(size = rel(1), colour = "blue")) +
          labs(x=NULL, y = NULL)
      }
      else {
        grid.info <- factor(rep(name.com, each = length(list.paths)), levels = name.com)
        dataset <- data.frame(row=row, col=col, circle.size=circle.size, significance=significance, exp=grid.info)
        ## build the BBplot
        bbplot <- ggplot(dataset, aes(y = factor(row),  x = factor(col))) +
          geom_point(data = subset(dataset, circle.size>0), aes(size = circle.size, alpha=significance), colour = colors[1]) +
          scale_alpha(guide="none") +
          facet_grid(facets=.~exp, scales="free_x", space="free_x")  +
          scale_size(range = c(1, 10)) + 
          labs(title = title) +
          theme_bw() +
          theme(axis.text.x = element_blank(), plot.title = element_text(size = rel(1), colour = "blue")) +
          labs(x=NULL, y = NULL)
      }
      return(bbplot)
    }
    else stop("The list of selected annoations is empty.")
  }
  else stop("Provide a list of DAVIDFunctionalAnnotationChart objects.")
}

Jplot <- function(david.obj.1,
                  david.obj.2,
                  max.pval = 0.01,
                  min.ngenes = 5,
                  title = "Jplot",
                  print.term = "full")
{
    if(inherits(david.obj.1, "DAVIDFunctionalAnnotationChart") & inherits(david.obj.2, "DAVIDFunctionalAnnotationChart")) {
        david.obj.1 <- david.obj.1[which(david.obj.1$PValue <= max.pval & david.obj.1$Count >= min.ngenes),]
        david.obj.2 <- david.obj.2[which(david.obj.2$PValue <= max.pval & david.obj.2$Count >= min.ngenes),]
        if(print.term == "name") {
          david.obj.1$Term <- unlist(lapply(david.obj.1$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
          david.obj.2$Term <- unlist(lapply(david.obj.2$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
        }
        else if(print.term == "description") {
            david.obj.1$Term <- unlist(lapply(david.obj.1$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
            david.obj.2$Term <- unlist(lapply(david.obj.2$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
        }
        else if(print.term != "full") 
          stop("Invalid input value for cod.term. Please, indicate a number between 0 and 2.")
        size.v <- length(david.obj.1$Term)*length(david.obj.2$Term)
        row <- rep(david.obj.1$Term, each=length(david.obj.2$Term))
        col <- rep(david.obj.2$Term, length(david.obj.1$Term))
        ij <- rep(0, length(david.obj.1$Term)*length(david.obj.2$Term))
        cnt <- 1
        for(i in seq_along(david.obj.1$Term)){
            for(j in seq_along(david.obj.2$Term)) {
                genes.i <- unlist(strsplit(david.obj.1$Genes[[i]], ", "))
                genes.j <- unlist(strsplit(david.obj.2$Genes[[j]], ", "))
                ij[cnt] <- length(intersect(genes.i,genes.j))/length(union(genes.i,genes.j))
                cnt <- cnt + 1
            }
        }
        data <- data.frame(row=factor(row), col=factor(col), ij=ij)
        data$ij<-cut(data$ij,
        breaks=c(0,0.01,0.25,0.5,0.75,0.99,1),
        include.lowest=TRUE,
        label=c("0%","10-25%","25-50%","50-75%","75-99%","1"))
        jplot <- ggplot(data, aes(x=row, y=col)) +
        geom_tile(aes(fill=ij),colour="black") +
        scale_fill_brewer(palette = "YlOrRd",name="Similarity score") +
        theme(axis.text.x=element_text(angle=-90, hjust=0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=20))
        return(jplot)
    }
    else stop("Provide two DAVIDFunctionalAnnotationChart objects.")
}


#' Function to visualise running enrichment score for a given sample and a gene set
#'
#' \code{visGSEA} is supposed to visualise running enrichment score for a given sample and a gene set. To help understand the underlying running enrichment score, the input gene scores are also displayed. Positions for members in the given gene set are color-coded in both displays (red line for the positive gene scores, and green line for the negative). 
#'
#' @param eTerm an object of class "eTerm"
#' @param which_sample which sample will be used. It can be index or sample names
#' @param which_term which term will be used. It can be index or term ID or term names 
#' @param weight type of score weigth. It can be "0" for unweighted (an equivalent to Kolmogorov-Smirnov, only considering the rank), "1" for weighted by input gene score (by default), and "2" for over-weighted, and so on
#' @param orientation the orientation of the plots. It can be either "vertical" (default) or "horizontal"
#' @param hit.linewidth the line width for the hits (ie genes in the gene set)
#' @param newpage logical to indicate whether to open a new page. By default, it sets to true for opening a new page
#' @return
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{dGSEA}}, \code{\link{dGSEAview}}
#' @include visGSEA.r
#' @examples
#' #visGSEA(eTerm, which_sample=1, which_term=1)

visGSEA <- function(eTerm, which_sample=1, which_term="GO:0006281", weight=1, orientation=c('vertical','horizontal'), hit.linewidth=0.5, newpage=T) 
{

    if (class(eTerm) != "eTerm" ){
        stop("The function must apply to a 'eTerm' object.\n")
    }
    
    orientation <- match.arg(orientation)
    
    data <- eTerm$data
    geneid <- rownames(data)
    nGene <- nrow(data)
    nSample <- ncol(data)
    sample_names <- colnames(data)
    
    gs <- eTerm$gs
    set_info <- eTerm$set_info
    nSet <- length(gs)

    ## determine which sample    
    if(which_sample==round(which_sample)){
        if(which_sample > nSample){
            which_sample <- nSample
        }else if(which_sample < 1){
            which_sample <- 1
        }
    }else{
        flag <- match(which_sample, sample_names)
        if(!is.na(flag)){
            which_sample <- flag
        }else{
            which_sample <- 1
        }
    }
    
    ## determine which gene set
    if(is.integer(which_term)){
        if(which_term > nSet){
            which_term <- nSet
        }else if(which_term < 1){
            which_term <- 1
        }
    }else{
        flag <- match(which_term, set_info$setID)
        if(!is.na(flag)){
            which_term <- flag
        }else{
            flag <- match(which_term, set_info$name)
            if(!is.na(flag)){
                which_term <- flag
            }else{
                which_term <- 1
            }
        }
    }
    
    
    rank.score <- data[,which_sample]
        
    ind <- order(rank.score, decreasing=T)
    rank.score.sorted <- rank.score[ind]
    geneid.sorted <- geneid[ind]
    geneid2ind <- ind

    nHit <- length(gs[[which_term]])
    nMiss <- nGene - nHit
            
    ## observed
    observed.point <- rep(-1/nMiss, nGene)
    flag <- match(gs[[which_term]], geneid.sorted)
    if(weight==0) {
        observed.point[flag] <- 1/nHit
    }else if(weight==1){
        hit_tmp <- abs(rank.score.sorted[flag])
        observed.point[flag] <- hit_tmp/sum(hit_tmp)
    }else{
        hit_tmp <- abs(rank.score.sorted[flag] ** weight)
        observed.point[flag] <- hit_tmp/sum(hit_tmp)      
    }
    RES <- cumsum(observed.point)
    max.RES <- max(RES)
    min.RES <- min(RES)
    es.observed <- signif(ifelse(max.RES>abs(min.RES), max.RES, min.RES), digits=5)
    es.position <- ifelse(max.RES>abs(min.RES), which.max(RES), which.min(RES))
    
    max.caxis <- ifelse(max.RES<0.4, 0.4, max.RES)
    min.caxis <- ifelse(min.RES>-0.4, -0.4, min.RES)
    caxis <- max(abs(max.caxis), abs(min.caxis))
    max.caxis <- caxis
    min.caxis <- -1*caxis
    
    
    ######################################################################################
    ## Visualisation
    if (newpage){
        grDevices::dev.new()
    }
    
    fontsize_axis <- 0.7
    
    if(orientation == "vertical"){

        ## For input data
        graphics::par(fig=c(0.05,1,0.4,1))
        plot(0,
            xlab="", xaxt="n",
            ylab="Input gene score",
            xlim=c(1, nGene), ylim=c(min(rank.score.sorted), max(rank.score.sorted)), 
            type="l", las=2, lwd=2, cex=1, cex.axis=fontsize_axis, col="black"
        )
    
        # value line 
        graphics::lines(1:nGene, rank.score.sorted, type="l", lwd=1, cex=1, col="black")
        # zero line
        graphics::lines(c(1,nGene), c(0,0), lwd=1, lty=1, cex=1, col="grey") 
        # left line
        graphics::lines(c(1,1), c(0,rank.score.sorted[1]), lwd=1, lty=1, cex=1, col="grey")
        # right line
        graphics::lines(c(nGene,nGene), c(0,rank.score.sorted[nGene]), lwd=1, lty=1, cex=1, col="grey")

        temp <- sapply(match(gs[[which_term]], geneid.sorted), function(x) {
            if(rank.score.sorted[x]>=0){
                graphics::lines(c(x, x), c(0, rank.score.sorted[x]), lwd=hit.linewidth, lty=1, cex=1, col="red")  # enrichment tags
            }else{
                graphics::lines(c(x, x), c(0, rank.score.sorted[x]), lwd=hit.linewidth, lty=1, cex=1, col="green")  # enrichment tags
            }
        })
    
        ## Running enrichment plot
        graphics::par(fig=c(0.05,1,0,0.6), new=TRUE)
        plot(1:nGene, RES,
            xlab=paste("Ranked genes (in a decreasing order) with ", length(gs[[which_term]]), " in geneset", sep="", collapse=""), xaxt="n",
            ylab="Running enrichment score",
            xlim=c(1, nGene), ylim=c(min.caxis, max.caxis), 
            type="l", las=2, lwd=2, cex=1, cex.axis=fontsize_axis, col="black"
        )
        graphics::axis(1, at=c(1, sum(rank.score.sorted>=0), nGene), cex.axis=fontsize_axis)
        graphics::lines(c(1, nGene), c(0, 0), lwd=1, lty=1, cex=1, col="grey") # zero RES line
    
        temp <- sapply(match(gs[[which_term]], geneid.sorted), function(x) {
            if(rank.score.sorted[x]>=0){
                graphics::lines(c(x, x), c(0, RES[x]), lwd=hit.linewidth, lty=1, cex=1, col="red")  # enrichment tags
            }else{
                graphics::lines(c(x, x), c(0, RES[x]), lwd=hit.linewidth, lty=1, cex=1, col="green")  # enrichment tags
            }
        })
    
        graphics::points(es.position, RES[es.position], lwd=2, col="blue") # at peak
    
        pvalue <- signif(eTerm$pvalue[which_term,which_sample],digits=2)
        adjp <- signif(eTerm$adjp[which_term,which_sample],digits=2)
    
        if(RES[es.position]<0){
    
            graphics::lines(c(es.position, es.position), c(RES[es.position], 0), lwd=2, lty=3, cex=1, col="blue") # peak line 
        
            leading <- sum(flag >= es.position) # number of leading genes
            leg.txt <- paste("Peak at rank=", es.position,
                             "\nLeading gene number=", leading,
                             "\nnormalized ES=", signif(eTerm$nes[which_term,which_sample],digits=3),
                             "\np-value=", ifelse(pvalue<0.01 & pvalue!=0, format(pvalue,scientific=T), pvalue),
                             "\nadjusted p-value=", ifelse(adjp<0.01 & adjp!=0, format(adjp,scientific=T), adjp),
                             sep="", collapse="")
            graphics::text(x=es.position, y=max.caxis, adj=c(0,1), col="blue", labels=leg.txt, cex=0.5)

        }else{
            graphics::lines(c(es.position, es.position), c(0, RES[es.position]), lwd=2, lty=3, cex=1, col="blue") # peak line
        
            leading <- sum(flag <= es.position) # number of leading genes
            leg.txt <- paste("Peak at rank=", es.position,
                             "\nLeading gene number=", leading,
                             "\nnormalized ES=", signif(eTerm$nes[which_term,which_sample],digits=3),
                             "\np-value=", ifelse(pvalue<0.01 & pvalue!=0, format(pvalue,scientific=T), pvalue),
                             "\nadjusted p-value=", ifelse(adjp<0.01 & adjp!=0, format(adjp,scientific=T), adjp),
                             sep="", collapse="")
            graphics::text(x=es.position, y=-1*max.caxis/10, adj=c(0,1), col="blue", labels=leg.txt, cex=0.5)        

        }
        
        graphics::mtext(paste("Sample:", sample_names[which_sample], "\nGeneset:", set_info$name[which_term]), side=3, outer=TRUE, line=-3, cex=1)
        
    }else if(orientation == "horizontal"){
    
        ## For input data
        graphics::par(fig=c(0.05,0.55,0,1))
        plot(0,
            ylab=paste("Ranked genes (in a decreasing order) with ", length(gs[[which_term]]), " in geneset", sep="", collapse=""), yaxt="n",
            xlab="Input gene score", 
            ylim=c(nGene, 1), xlim=c(min(rank.score.sorted), max(rank.score.sorted)), 
            type="l", lwd=2, cex=1, cex.axis=fontsize_axis, col="black"
        )
        graphics::axis(2, at=c(1, sum(rank.score.sorted>=0), nGene), las=2, cex.axis=fontsize_axis)
    
        if(0){
        for (x in seq(1, nGene, ceiling(nGene/500))) {
            if(rank.score.sorted[x]>=0){
                graphics::lines(c(0,rank.score.sorted[x]), c(x,x), lwd=1, cex=1, col="red")
            }else{
                graphics::lines(c(0,rank.score.sorted[x]), c(x,x), lwd=1, cex=1, col="green")
            }
        }
        }
    
        # value line 
        graphics::lines(rank.score.sorted, 1:nGene, type="l", lwd=1, cex=1, col="black") 
        # zero line
        graphics::lines(c(0,0), c(1,nGene), lwd=1, lty=1, cex=1, col="grey") 
        # top line
        graphics::lines(c(0,rank.score.sorted[1]), c(1,1), lwd=1, lty=1, cex=1, col="grey")
        # bottom line
        graphics::lines(c(0,rank.score.sorted[nGene]), c(nGene,nGene), lwd=1, lty=1, cex=1, col="grey")

        temp <- sapply(match(gs[[which_term]], geneid.sorted), function(x) {
            if(rank.score.sorted[x]>=0){
                graphics::lines(c(0, rank.score.sorted[x]), c(x, x), lwd=hit.linewidth, lty=1, cex=1, col="red")  # enrichment tags
            }else{
                graphics::lines(c(0, rank.score.sorted[x]), c(x, x), lwd=hit.linewidth, lty=1, cex=1, col="green")  # enrichment tags
            }
        })
    
        ## Running enrichment plot
        graphics::par(fig=c(0.4,1,0,1), new=TRUE)
        plot(RES, 1:nGene,
            ylab="", yaxt="n",
            xlab="Running enrichment score",
            ylim=c(nGene, 1), xlim=c(min.caxis, max.caxis), 
            type="l", lwd=2, cex=1, cex.axis=fontsize_axis, col="black"
            )
        graphics::lines(c(0, 0), c(1, nGene) , lwd=1, lty=1, cex=1, col="grey") # zero RES line
    
    
        temp <- sapply(match(gs[[which_term]], geneid.sorted), function(x) {
            if(rank.score.sorted[x]>=0){
                graphics::lines(c(0, RES[x]), c(x, x), lwd=hit.linewidth, lty=1, cex=1, col="red")  # enrichment tags
            }else{
                graphics::lines(c(0, RES[x]), c(x, x), lwd=hit.linewidth, lty=1, cex=1, col="green")  # enrichment tags
            }
        })
    
        graphics::points(RES[es.position], es.position, lwd=2, col="blue") # at peak
    
    
        pvalue <- signif(eTerm$pvalue[which_term,which_sample],digits=2)
        adjp <- signif(eTerm$adjp[which_term,which_sample],digits=2)
    
        if(RES[es.position]<0){
    
            graphics::lines(c(RES[es.position], 0), c(es.position, es.position), lwd=2, lty=3, cex=1, col="blue") # peak line 
        
            leading <- sum(flag >= es.position) # number of leading genes
            leg.txt <- paste("Peak at rank=", es.position,
                             "\nLeading gene number=", leading,
                             "\nnormalized ES=", signif(eTerm$nes[which_term,which_sample],digits=3),
                             "\np-value=", ifelse(pvalue<0.01 & pvalue!=0, format(pvalue,scientific=T), pvalue),
                             "\nadjusted p-value=", ifelse(adjp<0.01 & adjp!=0, format(adjp,scientific=T), adjp),
                             sep="", collapse="")
            graphics::text(x=max.caxis/10, y=es.position, adj=c(0,1), col="blue", labels=leg.txt, cex=0.5)

        }else{
            graphics::lines(c(0, RES[es.position]), c(es.position, es.position), lwd=2, lty=3, cex=1, col="blue") # peak line
        
            leading <- sum(flag <= es.position) # number of leading genes
            leg.txt <- paste("Peak at rank=", es.position,
                             "\nLeading gene number=", leading,
                             "\nnormalized ES=", signif(eTerm$nes[which_term,which_sample],digits=3),
                             "\np-value=", ifelse(pvalue<0.01 & pvalue!=0, format(pvalue,scientific=T), pvalue),
                             "\nadjusted p-value=", ifelse(adjp<0.01 & adjp!=0, format(adjp,scientific=T), adjp),
                             sep="", collapse="")
            graphics::text(x=min.caxis, y=es.position, adj=c(0,1), col="blue", labels=leg.txt, cex=0.5)        

        }
        
        graphics::mtext(paste("Sample:", sample_names[which_sample], "\nGeneset:", set_info$name[which_term]), side=3, outer=TRUE, line=-3, cex=1)
    }
    
    invisible()
}

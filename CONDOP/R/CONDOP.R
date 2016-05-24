#' Prepare data inputs for the main function \code{run.CONDOP()}.
#'
#' Load the annotation files and a list of count tables (or coverage vectors). 
#' Each count table is related to a specific experimental condition 
#' and it must contain two columns: fwd (coverage depth on the forward strand) 
#' and rev (coverage depth on the reverse strand). 
#' The annotations files are:
#'  
#'   - GFF-like file, it can be downloaded from the NCBI genomes ftp directory, ftp://ftp.ncbi.nih.gov/genomes.
#'   
#'   - DOOR-like file, it can be downloaded from http://csbl.bmb.uga.edu/DOOR/displayspecies.php.
#'   
#'   - FASTA-like file, it can be downloaded from www.ncbi.nlm.nih.gov.
#'  
#' @param gff.file A full local path indicating the GFF-like file to load <Gene annotations>. 
#' @param door.op.file A full local path indicating the DOOR-like file to load (DOOR-operon annotations).
#' @param fasta.file A full local path indicating the FASTA-like file to load or 
#'                   a character string representing the accession number of the genome sequence to download. 
#' @param list.cov.dat List of count tables.
#' @param remove.cov List of character values. Each charcater value corresponds to a specific type of annotated features. 
#'                                The coverage depth from those annotated feature will be removed. 
#'                                The default list contains "rRNA". The coverage depth of "rRNA" features will be removed. 
#' @param log2.expr Logical value indicating whether CONDOP will be using logged values of expression.
#'                  The expression values are compiled in RPKM values. Default logical value is TRUE.
#' @param sw Numeric value specifying the sliding window size. Default value is 100.
#' @param save.data.file Character string naming a file. The file will contain the input for the CONDOP main process.
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return A list of data inputs for the main process \code{run.CONDOP}.
#'  \item{genes.and.ops}{A merged dataframe containing information about genes/features and operons merged.} 
#'  \item{gseq}{A character vector representing the genome sequence of the target organism. }
#'  \item{igr.pos}{A dataframe containing information about intergenic regions (IRGs) - forward (+) strand.} 
#'  \item{igr.neg}{A dataframe containing information about intergenic regions (IRGs) - reverse (-) strand.} 
#'  \item{tl.cds}{A list of dataframes containing the expression levels of annotated coding sequences (CDS regions). One dataframe for each count table.}
#'  \item{tl.igr.pos}{A list of dataframes containing the expression levels of intergenic sequences (IGR regions) - forward (+) strand. One dataframe for each count table.}
#'  \item{tl.igr.neg}{A list of dataframes containing the expression levels of intergenic sequences (IGR regions) - reverse (-) strand. One dataframe for each count table.}
#'  \item{sid.points}{A list of dataframes containing information about boundaries of transcriptionally active regions.}
#'  \item{cut.lhe}{A list of numeric vectors indicating the cut-off values to distinguish low expressed RNA-seq data from high expression data on the forward and reverse strands. One dataframe for each count table.}
#'@examples
#' \dontrun{
#'     file_operon_annot <- system.file("extdata", "1944.opr", package="CONDOP")
#'     file_genome_seq   <- system.file("extdata", "EC-k12-MG1655.fasta", package="CONDOP")
#'     data(ct1)
#'     data.in <- pre.proc(file_genome_annot, file_operon_annot, "NC_000913", 
#'                         list.cov.dat = list(ct1 = ct1)) 
#' }
#' @note Use the \code{pre.proc} function before running \code{run.CONDOP}.
#'       You do not have to worry about how to make the input data structures for the the \code{run.CONDOP} function.
#' @author Vittorio Fortino
#' @importFrom seqinr read.fasta uco
#' @importFrom plyr ldply
#' @importFrom earth earth
#' @importFrom mclust Mclust
#' @importFrom S4Vectors metadata
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom stats cor.test ks.test mad median na.omit quantile sd
#' @importFrom utils URLdecode read.delim read.table write.table
#' @importClassesFrom IRanges IRanges
#' @import GenomicRanges
#' @export
pre.proc <- function(gff.file, door.op.file, fasta.file, list.cov.dat,
                     remove.cov = list("rRNA"), log2.expr = TRUE, sw = 100, 
                     save.data.file = NULL, verbose = TRUE) {
  if(!is.character(gff.file) & !is.character(door.op.file) & !is.character(fasta.file) & length(list.cov.dat) == 0) {
    stop("Input files are not specified.")
  }
  else {
      ## file.exists();
    if(!file.exists(gff.file))    stop("The GFF annotation file has not be specified.")
    if(!file.exists(door.op.file)) stop("The DOOR operon annotation file has not be specified.")
    if(!file.exists(fasta.file))	{
      if(!is.character(fasta.file))
        stop("Please specify the accession number or the FASTA file of the genome sequence has not be specified.")
      if(verbose) 
        cat("Searching and downloading the FASTA file.\n")
      gseq  <- get.NCBI.seq(fasta.file)
    }
    else{
      gseq  <- seqinr::read.fasta(fasta.file, seqtype="DNA")[[1]]
    }
    if(length(list.cov.dat) == 0)   stop("Coverage depth data has not be specified.")
    
    genes <- read.gff.annotations(gff.file, verbose=verbose)
    operons <- read.door.annotations(door.op.file, verbose=verbose)
  }
  igr.pos  <- get.intergenic.regions(genes, "+")
  igr.neg  <- get.intergenic.regions(genes, "-")
  genes.and.ops <- join.genes.and.operons(genes, operons) 
  not.ops <- names(table(genes.and.ops$operonID))[which(table(genes.and.ops$operonID) == 1)]
  id.not.ops <- which((genes.and.ops$operonID %in% not.ops) == TRUE)
  genes.and.ops$operonID[id.not.ops] <- NA
  if(verbose) {
    cat("Annotated feature types:\n")
    cat("      ",names(table(genes.and.ops$feature)))
    cat("\n\n")
  }
  if(verbose) cat("2) Compiling transcription levels...\n")	
  tl.cds <- list()
  tl.igr.pos <- list()
  tl.igr.neg <- list()
  for(i in 1:length(list.cov.dat)) {
    # if requested, remove depth-cov from a specific region
    if(!is.na(remove.cov) & length(remove.cov) > 0) {	
      for(j in 1:length(remove.cov)) {
        if(verbose) cat(paste("\nRemoving coverage depth from", remove.cov[[j]],"\n"))
        list.cov.dat[[i]] <- remove.cov.depth.from.aFeat(genes.and.ops, list.cov.dat[[i]], feature = remove.cov[[j]])
      }
    }
    tl.cds[[i]] <- comp.gene.transc.levels(genes.and.ops, list.cov.dat[[i]])  
    tl.igr.pos[[i]] <- comp.igr.transc.levels(genes.and.ops, igr.pos, list.cov.dat[[i]], tl.cds[[i]], "+");
    tl.igr.neg[[i]] <- comp.igr.transc.levels(genes.and.ops, igr.neg, list.cov.dat[[i]], tl.cds[[i]], "-");
    # set the 'expr' field
    if(log2.expr) {
      tl.cds[[i]]$expr <- tl.cds[[i]]$log2NormDepth
      tl.igr.pos[[i]]$expr <- tl.igr.pos[[i]]$log2NormDepth
      tl.igr.neg[[i]]$expr <- tl.igr.neg[[i]]$log2NormDepth
    }
    else {
      tl.cds[[i]]$expr <- tl.cds[[i]]$normDepth
      tl.igr.pos[[i]]$expr <- tl.igr.pos[[i]]$normDepth
      tl.igr.neg[[i]]$expr <- tl.igr.neg[[i]]$normDepth
    }
    if(verbose) {
      cat("\n  - quantiles in transcription level (CDSs +)...\n")
      cat(stats::quantile(tl.cds[[i]]$expr[tl.cds[[i]]$strand == "+"]))
      cat("\n  - quantiles in transcription level (CDSs -)...\n")
      cat(stats::quantile(tl.cds[[i]]$expr[tl.cds[[i]]$strand == "-"]))
      cat("\n  - quantiles in transcription level (IGRs +)...\n")
      cat(stats::quantile(tl.igr.pos[[i]]$expr))
      cat("\n  - quantiles in transcription level (IGRs -)...\n")
      cat(stats::quantile(tl.igr.neg[[i]]$expr))
      cat("\n")
    }
  }
  if(verbose) cat("\n3) Finding start/end transcription points...\n")	
  sid.points <- list()
  #if(parallel) {
  #  num.cores <- detectCores(all.tests = FALSE, logical = FALSE)
  #  registerDoMC(num.cores)
  #  sid.points <- foreach(i=1:length(list.cov.dat)) %dopar% {
  #    detect.sid.points(cd = list.cov.dat[[i]], sizeWindow = sw, verbose=FALSE)
  #  }
  #}
  #else{ 
    for(i in 1:length(list.cov.dat)) {
      if(verbose) cat(paste("   - transcriptome profiles:",i,"\n"))	
      sid.points[[i]] <- detect.sid.points(cd = list.cov.dat[[i]], sizeWindow = sw, verbose=FALSE)
    }
  #}
  if(verbose) cat("\n4) Calculating cut-off values...\n")
  cut.lhe <- list() # cut off values to distinguish between low and high expression
  for(i in 1:length(list.cov.dat)) {
    if(verbose) cat(paste("   - transcriptome profiles:",i,"\n"))	
    cut.lhe[[i]] <- qcut(list.cov.dat[[i]])
  }
  if(!is.null(save.data.file)) {
    if(verbose) cat("\n5) Saving data objects ...\n")	
    #save.data.tps(paste(save.tps.dir, "_", i, sep=""), sid.points[[i]])
    save(list = c("genes.and.ops","igr.pos","igr.neg",
                  "tl.cds","tl.igr.pos","tl.igr.neg",
                  "sid.points","cut.lhe"), file = save.data.file)
  }
  list(genes.and.ops = genes.and.ops, gseq = gseq, igr.pos = igr.pos, igr.neg = igr.neg, # annotations    
       tl.cds = tl.cds, tl.igr.pos = tl.igr.pos, tl.igr.neg = tl.igr.neg,                # transc. levels  
       sid.points = sid.points, cut.lhe = cut.lhe)                                       # transc. start/end points (+ cut off values)
}

#' Build condition-dependent operon maps.
#'
#' It develops an ensemble operon pair classifier that combines 
#' both genomic and transcriptomic features. 
#' The ensemble classifier consists of three machine-learning models 
#' that are trained on a small set of confirmed operon pairs (OPs) and 
#' non-operon pairs (NOPs). The set of OPs and NOPs is identified by crosschecking 
#' the DOOR annotation with consecutive, active coding-sequence and intergenic regions,
#' indicated with CDSs and IGR respectively. The trained ensemble classifier is used 
#' to predict the operon status of all the gene-pairs, including DOOR-based operon pairs, 
#' namely DOPs, and putative operon pairs (POPs). 
#' Finally, a linkage process is exploited to combine consecutive 
#' operon-pairs classified as OP, and to build the map of condition-dependent operons.
#'
#' @param data.in The output of the \code{pre.proc} function.
#' @param bkgExprCDS A threshold to be used for finding active coding-sequence regions. Default values is 0.1.
#' @param bkgExprIGR A threshold to be used for finding the active/transcribed intergenic regions. Default values is 0.25.
#' @param maxLenIGR Maximum length for the intergenic regions. Default values is 150.
#' @param win.start.trp Specify the maximum and the minimum distance from the beginning of a coding region. It is important to validate transcription start points. Defauls values are 100 (max) and 10 (min).
#' @param win.end.trp Specify the minimum and maximum distance from the end of a coding region. It is important to validate transcription end points. Defauls values are 10 (min) and 100 (max).
#' @param norm.type Character vector indicating the method to use for the normalization step. Default value is "n1".
#'                   n0 - without normalization; n1 - standardization ((x-mean)/sd); 
#'                   n2 - positional standardization ((x-median)/mad); n3 - unitization ((x-mean)/range);
#'                   n4 - unitization with zero minimum ((x-min)/range); n5 - normalization in range <-1,1> ((x-mean)/max(abs(x-mean))).
#' @param cl.run Number of runs of training/validation. Default values is 30.
#' @param nfolds Indicate the number of folds to be used for the cross-validation step. Default values is 5.
#' @param cons Indicate the minimum number of positive votes necessary to classify a gene pair as operon pair. Default values is 2.
#' @param find.ext To add putative operon pairs classified as OP to the condition-dependent operon map. Defaults to FALSE.
#' @param save.TAB.file Character string naming a file. The final condition operon map is saved in a tab-delimeted text file. Default values is NULL - the cond. operon map is not saved.
#' @param save.BED.file Character string naming a file. The final condition operon map is saved in a BED-like file. Default values is NULL - the cond. operon map is not saved.
#' @param return.all Logical value indicating if extra data must be provided in output.
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return List of data structures built by CONDOP.
#' If \code{return.all} is FALSE:
#'  \item{ndata}{A list of dataframes containing OPs and NOPs used for the traing/validation process. One for each count table.}
#'  \item{cls}{A list of OP classifiers for each count table.}
#'  \item{ev.cls}{A data.frame containing the evalaution result for the trained classifiers. One for each count table.}
#'  \item{pred.TS}{A list of dataframes containing the classification results on the training set. One for each count table.}
#'  \item{pred.POPs}{A list of dataframes containing the prediction results on the POPs. One for each count table.}
#'  \item{pred.DOPs}{A list of dataframes containg the prediction results on the DOPs. One for each count table.}
#'  \item{comap}{A list of condition-dependent operon maps (comaps). One for each count table.}
#'  \item{info}{A list of generic information on the confirmed DOOR based operons. One for each count table.}
#' If \code{return.all} is TRUE the \code{run.CONDOP()} function also provides..
#'  \item{osp}{A list of dataframes containing confirmed operon start points. One for each count table.}
#'  \item{oep}{A list of dataframes containing confirmed operon end points. One for each count table.}
#'  \item{cops}{A list of dataframes containing confirmed operons. One for each count table.}
#'  \item{OPs}{A list of dataframes containing OPs. One for each count table.}
#'  \item{NOPs}{A list of dataframes containing NOPs. One for each count table.}
#'  \item{POPs}{A list of dataframes containing POPs. One for each count table.}
#'  \item{DOPs}{A list of dataframes containing DOPs. One for each count table.}
#' @examples
#' \dontrun{
#'     file_operon_annot <- system.file("extdata", "1944.opr", package="CONDOP")
#'     file_genome_seq   <- system.file("extdata", "EC-k12-MG1655.fasta", package="CONDOP")
#'     data(ct1)
#'     data.in   <- pre.proc(file_genome_annot, file_operon_annot, "NC_000913",
#'                           list.cov.dat = list(ct1 = ct1)) 
#'     res.comap <- run.CONDOP(data.in = data.in, bkgExprCDS = 0.2, bkgExprIGR = 0.2, 
#'                             maxLenIGR = 150, find.ext = TRUE)                      
#' }
#' @author Vittorio Fortino
#' @importFrom randomForest randomForest
#' @importFrom randomForest importance
#' @import rminer
#' @export
run.CONDOP <- function(data.in,
                       bkgExprCDS = 0.1, bkgExprIGR = 0.25, maxLenIGR = 150, win.start.trp = c(100,10), win.end.trp = c(10,100),
                       norm.type = "n1", cl.run = 30, nfolds = 5, cons = 2, find.ext = FALSE, save.TAB.file = NULL, save.BED.file = NULL, 
                       return.all = FALSE, verbose = TRUE) {
  if(length(data.in) == 0) {
    stop("Missing input data.\n")
  }	
  osp  <- list()
  oep  <- list()
  cops <- list()
  OPs  <- list()
  NOPs <- list()
  POPs <- list()
  DOPs <- list()
  ndata <- list()	
  models <- list()
  evals  <- list()
  pred.TS  <- list()
  pred.POPs  <- list()
  pred.DOPs  <- list()
  comap <- list()
  info <- list()
  num.transc.profiles <- length(data.in$sid.points)
  for(i in 1:num.transc.profiles) {
    if(verbose) cat("TRANSCRIPTOME PROFILE #",i,"\n",sep="")	
    if(verbose) cat("1) Determine operon start- and end-points (OSPs and OEPs)...\n")	
    osp[[i]] <- get.operon.start.points(data.in$sid.points[[i]]$fwd.incs, 
                                        data.in$sid.points[[i]]$rev.incs, 
                                        data.in$genes.and.ops, data.in$igr.pos, data.in$igr.neg, 
                                        data.in$tl.cds[[i]], 
                                        borders = win.start.trp, 
                                        max.start.transc = data.in$cut.lhe[[i]],
                                        minExprCDS = bkgExprCDS,
                                        verbose = verbose) 							
    oep[[i]] <- get.operon.end.points(data.in$sid.points[[i]]$fwd.decs, 
                                      data.in$sid.points[[i]]$rev.decs, 
                                      data.in$genes.and.ops, data.in$igr.pos, data.in$igr.neg, 
                                      data.in$tl.cds[[i]], 
                                      borders = win.end.trp, 
                                      max.end.transc = data.in$cut.lhe[[i]],
                                      minExprCDS = bkgExprCDS,
                                      verbose = verbose) 		
    if(verbose) cat("\n2) Verify confirmed operons...\n")
    cops[[i]] <- compile.confirmed.operons(data.in$genes.and.ops, 
                                           data.in$tl.cds[[i]], 
                                           data.in$tl.igr.pos[[i]], data.in$tl.igr.neg[[i]], 
                                           osp[[i]]$POSSs, oep[[i]]$POESs, 
                                           minExprCDS = bkgExprCDS, minExprIGR = bkgExprIGR,
                                           max.start.transc = data.in$cut.lhe[[i]], 
                                           max.end.transc = data.in$cut.lhe[[i]],
                                           verbose = verbose)
    if(verbose) cat("\n3) Compile the features for OPs and NOPs...\n")
    OPs[[i]]  <-	select.ops(cops[[i]], data.in$genes.and.ops, data.in$tl.cds[[i]], 
                             data.in$tl.igr.pos[[i]], data.in$tl.igr.neg[[i]], data.in$gseq,
                             verbose = verbose)		   
    NOPs[[i]] <-	select.nops(data.in$genes.and.ops, OPs[[i]], data.in$tl.cds[[i]], 
                              data.in$tl.igr.pos[[i]], data.in$tl.igr.neg[[i]], 
                              osp[[i]]$POSSs, oep[[i]]$POESs, data.in$gseq, 
                              max.start.transc = data.in$cut.lhe[[i]], 
                              max.end.transc = data.in$cut.lhe[[i]],
                              verbose = verbose) 
    
    if(verbose) cat("\n4) Find gene pairs with an operon status to (re)define...\n")
    POPs[[i]] <- 	select.pops(OPs[[i]], NOPs[[i]], data.in$genes.and.ops, data.in$tl.cds[[i]], 
                              data.in$tl.igr.pos[[i]], data.in$tl.igr.neg[[i]], 
                              minExprCDS = bkgExprCDS, minExprIGR = bkgExprIGR, 
                              maxLenIGR = maxLenIGR, osp[[i]]$POSSs, oep[[i]]$POESs,
                              max.start.transc = data.in$cut.lhe[[i]], 
                              max.end.transc = data.in$cut.lhe[[i]], data.in$gseq,
                              verbose = verbose) 
    
    if(verbose) cat("\n5) Find gene pairs annotated in DOOR database that were not confirmed... \n")
    DOPs[[i]] <-	select.ops.indoor(data.in$genes.and.ops, data.in$tl.cds[[i]], 
                                   data.in$tl.igr.pos[[i]], data.in$tl.igr.neg[[i]], data.in$gseq,
                                   verbose = verbose)
    
    if(verbose) cat("\n6) Normalize the feature values for OPs, NOPs, POPs and DOPs...\n")
    ndata[[i]] <- pre.processing(OPs[[i]], NOPs[[i]], 
                                 POPs[[i]], DOPs[[i]], 
                                 type = norm.type,
                                 verbose = verbose)		
    
    if(verbose) cat("\n7) Tune and train the classification models...\n")
    models[[i]] <- tune.cls(data=ndata[[i]]$dat[,-ncol(ndata[[i]]$dat)], 
                            class=ndata[[i]]$dat$class,
                            verbose = verbose)
    
    if(verbose) cat("\n8) Validate the classification models...\n")
    evals[[i]] <- validate.cls(data=ndata[[i]]$dat[,-ncol(ndata[[i]]$dat)], 
                               class=ndata[[i]]$dat$class, models[[i]], 
                               runs = cl.run, kf = nfolds,
                               verbose = verbose)
    
    if(verbose) print(evals[[i]])
    
    #if(verbose) cat("\n7) Train the operon classifier...\n")
    #cls[[i]] <- train.RFs(data=ndata[[i]]$dat[,-ncol(ndata[[i]]$dat)],
    #        							   class=ndata[[i]]$dat$class, run = cl.run, 
    #        							   p = training, ntr = rf.ntr, mtree=rf.mtree)
    
    if(verbose) cat("\n9) Predict the operon status of POPs and DOPs\n")
    pred.TS[[i]]   <- pred.operon.status(models[[i]], ndata[[i]]$dat, 
                                         type="TS", cons = cons,
                                         verbose = verbose)
    pred.POPs[[i]] <- pred.operon.status(models[[i]], ndata[[i]]$dat.pops, 
                                         type="POPs", cons = cons,
                                         verbose = verbose)
    pred.DOPs[[i]] <- pred.operon.status(models[[i]], ndata[[i]]$dat.dops, 
                                         type="DOPs", cons = cons,
                                         verbose = verbose)
    pred.DOPs[[i]] <- data.frame(rownames(pred.DOPs[[i]]), op.st = pred.DOPs[[i]]$cons, 
                                 op.ref = DOPs[[i]]$refOp, row.names = 1)
    pred.POPs[[i]] <- data.frame(rownames(pred.POPs[[i]]), op.st = pred.POPs[[i]]$cons, 
                                 op.ref1 = POPs[[i]]$refOp1, op.ref2 = POPs[[i]]$refOp2, row.names=1) 
    if(verbose) cat("\n10) Build the condition dependent operon map\n")
    
    if(!is.null(save.BED.file))
      comap[[i]] <- getCondOperonMap(pred.POPs[[i]], pred.DOPs[[i]], data.in$genes.and.ops,
                                     osp[[i]]$POSSs, oep[[i]]$POESs,
                                     max.start.transc = data.in$cut.lhe[[i]],
                                     max.end.transc = data.in$cut.lhe[[i]],
                                     find.ext = find.ext,
                                     BED.file = paste(save.BED.file, i, "bed", sep="."),
                                     verbose = verbose)
    else if(!is.null(save.TAB.file))
      comap[[i]] <- getCondOperonMap(pred.POPs[[i]], pred.DOPs[[i]], data.in$genes.and.ops,
                                     osp[[i]]$POSSs, oep[[i]]$POESs,
                                     max.start.transc = data.in$cut.lhe[[i]],
                                     max.end.transc = data.in$cut.lhe[[i]],
                                     find.ext = find.ext,
                                     TAB.file = paste(save.TAB.file, i, "txt", sep="."),
                                     verbose = verbose)
    else 
      comap[[i]] <- getCondOperonMap(pred.POPs[[i]], pred.DOPs[[i]], data.in$genes.and.ops,
                                     osp[[i]]$POSSs, oep[[i]]$POESs,
                                     max.start.transc = data.in$cut.lhe[[i]],
                                     max.end.transc = data.in$cut.lhe[[i]],
                                     find.ext = find.ext,
                                     verbose = verbose)
    
    if(verbose) cat("\n11) Compile generic information about the condition-dependent operon predictions\n")
    
    info[[i]] <- get.info(data.in$genes.and.ops, evals[[i]], pred.POPs[[i]], pred.DOPs[[i]], comap[[i]], verbose = verbose)
    
    if(verbose) cat("\n")
  }
  #if(!is.null(save.data.file)) {
  #  if(verbose) cat("\nSaving data objects...\n")
  #  save(list = c("osp","oep","cops","OPs","NOPs","POPs","DOPs",
  #                "ndata","cls","pred.POPs","pred.DOPs","comap"), 
  #       file = save.data.file)
  #}
  if(!return.all)	
    return(list(ndata=ndata, cls=models, pred.TS=pred.TS, pred.POPs=pred.POPs, pred.DOPs=pred.DOPs, 
                comap = comap, info=info))
  else
    return(list(ndata=ndata, cls=models, ev.cls = evals, pred.TS=pred.TS, pred.POPs=pred.POPs, pred.DOPs=pred.DOPs, 
                osp = osp, oep = oep, cops = cops, OPs = OPs, NOPs = NOPs, POPs = POPs, DOPs = DOPs, comap=comap,
                info=info))
}

#' Find and get a genome sequence by specifyng the accession number. Read gene annotations in GFF format
#'
#' @param gff_file A full local path indicating the GFF file to load.
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return Annotation data table.
#' @keywords internal
#' @author Vittorio Fortino
#' @importFrom seqinr choosebank getSequence query closebank
#' @export
get.NCBI.seq <- function(accession,  list.dbs = NULL, verbose = TRUE) {
    if(length(list.dbs) == 0) list.dbs <- seqinr::choosebank()
    # first find which ACNUC database the accession is stored in:
    for (db in list.dbs){
        cat(c("Looking in", db, "...\n"))
        seqinr::choosebank(db)
        # check if the sequence is in ACNUC database 'db':
        resquery <- try(seqinr::query(".tmpquery", paste("AC=", accession)), silent = TRUE)
        if (!(inherits(resquery, "try-error"))) {
            query2 <- seqinr::query("query2", paste("AC=",accession,sep=""))
            # see if a sequence was retrieved:
            seq <- seqinr::getSequence(query2$req[[1]])
            closebank()
            return(tolower(seq))
        }
        closebank()
    }
    if(verbose) cat("\nERROR: accession",accession,"was not found")
    NULL
}

#' Read gene annotations in GFF format 
#'
#' Internal function to read GFF data file downloaded from the NCBI genomes ftp directory, 
#' ftp://ftp.ncbi.nih.gov/genomes.
#' @param gff_file A full local path indicating the GFF file to load.
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return Annotation data table.
#' @keywords internal
#' @author Vittorio Fortino
#' read.gff.annotations()
read.gff.annotations <- function(gff_file, verbose=TRUE, ...) {
  gff.annot <- GenomicRanges::as.data.frame(read.annot.from.gff(gff_file))
  if(verbose) cat("\nGene annoatiotions loaded.")
  data.frame(locusTag = gff.annot$locus, feature = as.character(gff.annot$feature),
             start = gff.annot$start, end = gff.annot$end,
             strand = as.character(gff.annot$strand), gid = as.character(gff.annot$gid), 
             gene = as.character(gff.annot$gene),
             row.names = 1, stringsAsFactors = FALSE)
}	

#' Read a GFF file from NCBI and return a GRanges object.
#'
#' @param gff.file a GFF file
#' @param locus.tags only return genes with locus tags. Defaults to TRUE.
#' @param nrows number of rows to read.
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return GRanges with 4 elementMetadata columns: locus, proetin id, gene id, feature, description and gene name. 
#'         If all rows are returned (locus.tags=FALSE), then score, phase and tags are included. 
#' @author Vittorio Fortino
#' read.annot.from.gff()
read.annot.from.gff <- function(gff.file, locus.tags = TRUE, nrows = -1, verbose=TRUE) {
  x <- utils::read.delim(gff.file, stringsAsFactors = FALSE, comment.char = "#", 
                  header = FALSE, nrows = nrows)
  colnames(x) <- c("seqid", "source", "feature", "start", "end", 
                   "score", "strand", "phase", "tags")
  n <- x$strand == "."
  if (any(n)) {
    x$strand[n] <- "+"
    print("WARNING: changing '.' to '+' in strand column")
  }
  seqid <- unique(x$seqid)
  if (all(grepl("\\.[0-9]$", seqid))) {
    #x$seqid <- genomes::strsplit2(x$seqid, ".", fixed = TRUE)
    x.split <- strsplit(as.character(x$seqid), split = ".", fixed = TRUE)
    x$seqid <- sapply(x.split, "[", n = 1)
  }
  x$id <- gsub("ID=([^;]*).*", "\\1", x$tags)
  if (!locus.tags) {
    gff <- GenomicRanges::GRanges(seqnames = x$seqid, ranges = IRanges::IRanges(x$start, x$end), 
                                  strand = x$strand, data.frame(x[, c(10, 3, 6, 8, 9)]))
  }
  else {
    #n <- x$tags %like% "*locus_tag=*"
    n <-  grepl("*locus_tag=*", x$tags)
    genes <- subset(x, n)
    y <- subset(x, !n)
    y$parent <- NA
    #n <- y$tags %like% "*Parent=*"
    n <- grepl("*Parent=*", y$tags)
    y$parent[n] <- gsub(".*Parent=([^;]*).*", "\\1", y$tags[n])
    y$gid <- ""
    #n <- y$tags %like% "*GeneID=*"
    n <- grepl("*GeneID:*", y$tags)
    y$gid[n] <- gsub(".*GeneID:([^;]*).*", "\\1", y$tags[n])
    y$pid <- ""
    #n <- y$tags %like% "*protein_id=*"
    n <- grepl("*protein_id=*", y$tags)
    if (sum(n) > 0) 
      y$pid[n] <- gsub(".*protein_id=([^;.]*).*", "\\1",  y$tags[n])
    y$product <- ""
    #n <- y$tags %like% "*product=*"
    n <- grepl("*product=*",y$tags)
    if (sum(n) > 0) 
      y$product[n] <- gsub(".*product=([^;]*).*", "\\1",y$tags[n])
    #n <- is.na(y$product) & y$tags %like% "*Note=*"
    n <- grepl("*Note=*", y$tags)
    n <- is.na(y$product) & n
    y$product[n] <- gsub(".*Note=([^;]*).*", "\\1", y$tags[n])
    n <- grep("%", y$product)
    if (length(n) > 0) 
      y$product[n] <- as.vector(sapply(y$product[n], utils::URLdecode))
    genes$locus <- gsub(".*;locus_tag=([^;]*).*", "\\1", genes$tags)
    n <- match(genes$id, y$parent)
    n2 <- !is.na(n)
    genes$feature[n2] <- y$feature[n[n2]]
    genes$description <- ""
    genes$description[n2] <- y$product[n[n2]]
    genes$gid[n2] <- y$gid[n[n2]]
    genes$pid[n2] <- y$pid[n[n2]]
    genes$feature[grepl("*pseudo=true*", genes$tags)] <- "pseudo"
    n <- which(genes$feature == "gene")
    if (length(n) > 0) {
      n2 <- match(paste(genes$start[n], genes$end[n]), 
                  paste(y$start, y$end))
      n3 <- !is.na(n2)
      genes$feature[n[n3]] <- y$feature[n2[n3]]
      genes$description[n[n3]] <- y$product[n2[n3]]
    }
    genes$feature[genes$feature == "transcript"] <- "miscRNA"
    genes$feature[genes$feature == "region"] <- "other"
    genes$feature[genes$feature == "gene"] <- "other"
    genes$gene <- ""
    n <- grep("gene=", genes$tag)
    if (length(n) > 0) 
      genes$gene[n] <- gsub(".*gene=([^;]*).*", "\\1", genes$tags[n])
    n <- grep("gene_synonym=", genes$tag)
    if (length(n) > 0) 
      genes$gene[n] <- paste(genes$gene[n], gsub(".*gene_synonym=([^;]*).*", 
                                                 "\\1", genes$tags[n]), sep = ",")
    n <- grep("%", genes$gene)
    if (length(n) > 0) 
      genes$gene[n] <- as.vector(sapply(genes$gene[n], utils::URLdecode))
    if (any(table(genes$locus) > 1)) {
      print("Warning: grouping coordinates for duplicate locus tags")
      genes2 <- by(genes, genes$locus, function(x) {
        data.frame(locus = x$locus[1], seqid = x$seqid[1], 
                   start = min(x$start), end = max(x$end), strand = x$strand[1], gid = x$gid[1],
                   pid = x$pid[1], feature = x$feature[1], description = x$description[1], 
                   gene = x$gene[1])
      })
      genes <- do.call("rbind", genes2)
    }
    if (any(diff(order(genes$start)) != 1)) {
      genes <- genes[order(genes$start), ]
    }
    gff <- GenomicRanges::GRanges(seqnames = genes$seqid, ranges = IRanges::IRanges(genes$start,genes$end), 
                                  strand = genes$strand, genes[, c("locus", "pid", "gid", "feature", "description", "gene")])
  }
  if (length(seqid) == 1) GenomeInfoDb::seqlengths(gff) <- max(x$end)
  S4Vectors::metadata(gff) <- list(source = unique(x$source), defline = seqid)
  gff
}

#' Read the operon(s) data file downloaded from http://csbl.bmb.uga.edu/DOOR/displayspecies.php
#'
#' Internal function to read the operon(s) data file downloaded from http://csbl.bmb.uga.edu/DOOR/displayspecies.php
#' @param opr_file A full local path indicating the .opr file to load. 
#' @param verbose Indicate whether information about the process should be reported. Defaults to TRUE.
#' @return Operon data table.
#' @keywords internal
#' @author Vittorio Fortino
#' read.door.annotations()
read.door.annotations <- function(opr_file, verbose=TRUE) {   		
  opr.info <- NULL;
  if(file.exists(opr_file)) {
    colNames <- c("operonID", "gi", "locusTag", "start", "end", 
                  "strand", "length", "cog", "description")
    colClasses <- c("character","character", "character", "integer", "integer",    
                    "character", "numeric", "character", "character")
    opr.info <- utils::read.table(opr_file, sep="\t", header = TRUE, row.names = 3, quote = "",
                           col.names = colNames, colClasses = colClasses, stringsAsFactors = FALSE)
  }
  if(verbose)	cat("\nDOOR operon annoatiotions loaded.\n")
  return(opr.info)
}

#' Join gene(s) and operon(s) annotations. 
#'
#' Internal function to merge the annotations with DOOR data.
#' @param genes An annotation data table.
#' @param ops An operon data table.
#' @return Operon data table.
#' @keywords internal
#' @author Vittorio Fortino
#' join.genes.and.operons()
join.genes.and.operons <- function(genes, ops) {     
  #operonID  <- vector(mode = "character", length = nrow(genes))
  #gi  <- vector(mode = "character", length = nrow(genes))
  #cog  <- vector(mode = "character", length = nrow(genes))
  for(i in 1:nrow(genes)) {
    id <- which(rownames(genes[i,]) == rownames(ops))
    if(length(id) > 0) {
      genes$operonID[i] <- ops$operonID[id]
    }
    else {
      genes$operonID[i] <- NA
    }
  }
  return(genes)
}

#' Build a data table containing generic information on intergenic regions.
#'
#' Internal function to build a data table containing information of the intergenic regions on a given strand.
#' @param genes An annotation data table.
#' @param str A given strand. Defaults to "+".
#' @keywords internal
#' @author Vittorio Fortino
#' get.intergenic.regions()
get.intergenic.regions <- function(genes, str = "+") {  
  genes <- genes[which(genes$strand == str),]
  numIGRs <- nrow(genes) - 1
  igr_lengths <- vector(mode = "numeric", length = numIGRs)
  startLocusTags  <- vector(mode = "character", length = numIGRs)
  endLocusTags  <- vector(mode = "character", length = numIGRs)
  start  <- vector(mode = "numeric", length = numIGRs)
  end    <- vector(mode = "numeric", length = numIGRs)
  #type   <- vector(mode = "character", length = numIGRs) 
  for(i in 1:numIGRs) {
    igr_lengths[i]  <- (genes$start[i+1] - genes$end[i]) - 1
    # default value
    type <- "distant"
    startLocusTags[i] <- rownames(genes)[i]
    endLocusTags[i] <- rownames(genes)[i+1]
    start[i] <- genes$end[i] + 1
    end[i]   <- genes$start[i+1] - 1
  }
  intergenic_regions <- data.frame(startLocusTag = startLocusTags, 
                                   endLocusTag = endLocusTags,
                                   length = igr_lengths, 
                                   start = start, 
                                   end = end,
                                   type = type,
                                   stringsAsFactors = FALSE)
  pos_overlapping_genes <- which(intergenic_regions$length < 0)
  pos_adjacent_genes <- which(intergenic_regions$length == 0)
  intergenic_regions$type[pos_overlapping_genes] = "O"           ## overlapping
  intergenic_regions$type[pos_adjacent_genes] = "A"              ## adjacent
  ## new categories for different types of distance
  q <- stats::quantile(intergenic_regions$length)	
  ##   - PE (Poorly Distant):  0 - q25%				
  ##   - ME (Moderately Distant): q25% - q50%						
  ##   - HE (Highly Distant): q50% - q75%		
  ##   - VE (Very Highly Distant): >= q75%		
  poorly_distant_igr <- which(intergenic_regions$length > 0 & intergenic_regions$length < q[2])
  intergenic_regions$type[poorly_distant_igr] = "PD"
  moderately_distant_igr <- which(intergenic_regions$length >= q[2] & intergenic_regions$length < q[3])
  intergenic_regions$type[moderately_distant_igr] = "MD"
  highly_distant_igr <- which(intergenic_regions$length >= q[3] & intergenic_regions$length < q[4])
  intergenic_regions$type[highly_distant_igr] = "HD"
  very_highly_distant_igr <- which(intergenic_regions$length >= q[4])
  intergenic_regions$type[very_highly_distant_igr] = "VD"
  intergenic_regions$type <- as.factor(intergenic_regions$type)
  return(intergenic_regions)
}

#' Remove the read coverage on a given feature (e.g. rRNA and tRNA).
#'
#' Internal function to remove the coverage depth from a given features.
#' @param genes An annotation data table.
#' @param cd A data table containing the coverage depth of an RNA-seq expression profile(s).
#' @param feature A given feature type. Defaults to "rRNA".
#' @keywords internal
#' @author Vittorio Fortino
#' remove.cov.depth.from.aFeat()
remove.cov.depth.from.aFeat <- function(genes, cd, feature = "rRNA") { 
  for(i in 1:nrow(genes)) {
    if(genes$feature[i] == feature) { 
      if(genes$strand[i] == "+") { # FWD strand
        cd$fwd[genes$start[i]:genes$end[i]] <- 0
      }
      else { # REV strand
        cd$rev[genes$start[i]:genes$end[i]] <- 0
      }
    }
  }
  return(cd)
}

#' Compile the transcription levels for the coding regions.
#'
#' Internal function to quantify the expression at CDS-level.
#' @param genes An annotation data table.
#' @param cd A data table representing the coverage depth for a given RNA-seq expression profile.
#' @keywords internal
#' @author Vittorio Fortino
#' comp.gene.transc.levels()
comp.gene.transc.levels <- function(genes, cd) {  
  num_genes <- nrow(genes)
  readsGene		   <- vector(mode = "numeric", length = num_genes)
  normDepth 	   <- vector(mode = "numeric", length = num_genes)
  log2NormDepth	 <- vector(mode = "numeric", length = num_genes)
  tot.cov <- sum(cd$fwd) + sum(cd$rev)
  for(i in 1:num_genes) {
    if(genes$strand[i] == "+") { # FWD strand
      readsGene[i] <- sum(cd$fwd[genes$start[i]:genes$end[i]])
    }
    else { # REV strand
      readsGene[i] <- sum(cd$rev[genes$start[i]:genes$end[i]])
    }
    normDepth[i] <- readsGene[i] / ((tot.cov * ((genes$end[i]-genes$start[i])+1))/ 10^9)   # RPKM normalization
    log2NormDepth[i] <- 0
    if(normDepth[i] > 1) log2NormDepth[i]  <- log2(normDepth[i])
  }
  return(data.frame(locusTag = rownames(genes), 
                    readsGene = readsGene,
                    normDepth = round(normDepth, digits=4),
                    log2NormDepth = round(log2NormDepth, digits=4),
                    strand = genes$strand,
                    row.names = 1,
                    stringsAsFactors = FALSE))					   
}

#' Compile the transcription levels for the intergenic regions.
#'
#' Internal function to quantify expression at the intergenic-level.
#' @param genes Data table containing gene annotations.
#' @param genes Data tbale containing intergenic region annotations. 
#' @param cd Data table containing the coverage depth of an RNA-seq expression profile(s).
#' @param transcCDSs Transcription levels for the coding regions.
#' @param str Charcater value indicating the strand. Defauls value is "+".
#' @keywords internal
#' @author Vittorio Fortino
#' comp.igr.transc.levels()
comp.igr.transc.levels  <- function(genes, igrs, cd, transcCDSs, str="+") {  
  normDepth      <- vector(mode = "numeric", length = nrow(igrs))
  log2NormDepth  <- vector(mode = "numeric", length = nrow(igrs))
  tot.cov <- sum(cd$fwd) + sum(cd$rev)
  for(i in 1:nrow(igrs)) {
    if(igrs$length[i] > 0) { ## distant genes
      if(str == '+') {
        readsGene    <- sum(cd$fwd[igrs$start[i]:igrs$end[i]])
        normDepth[i] <- readsGene / ((tot.cov * ((igrs$end[i]-igrs$start[i])+1))/ 10^9)   # RPKM normalization
      }
      else {
        readsGene    <- sum(cd$rev[igrs$start[i]:igrs$end[i]])
        normDepth[i] <- readsGene / ((tot.cov * ((igrs$end[i]-igrs$start[i])+1))/ 10^9)   # RPKM normalization
      }
    }
    else {
      normDepth[i] <- 0
    }
  }
  log2NormDepth  <- ifelse(normDepth > 1, log2(normDepth), 0)
  log2NormDepth  <- ifelse(!is.infinite(log2NormDepth),log2NormDepth,0)
  log2NormDepth  <- ifelse(!is.nan(log2NormDepth),log2NormDepth,0)
  igrs$normDepth     <- round(normDepth, digits = 4)
  igrs$log2NormDepth <- round(log2NormDepth, digits = 4)
  ## adjust the transcription level for the overlapping regions
  for(i in 1:nrow(igrs)) {
    igr <- igrs[i,]
    if(igr$type == "O") {
      idG1 <- which(rownames(genes) == igr$startLocusTag)
      idG2 <- which(rownames(genes) == igr$endLocusTag)
      lengthG1 <- (genes$end[idG1]-genes$start[idG1])+1
      lengthG2 <- (genes$end[idG2]-genes$start[idG2])+1
      # log2NormDepth
      exprG1 <- transcCDSs$log2NormDepth[which(rownames(transcCDSs) == igr$startLocusTag)]
      exprG2 <- transcCDSs$log2NormDepth[which(rownames(transcCDSs) == igr$endLocusTag)]
      diff <- abs(exprG1 - exprG2) 
      weight <- diff * (abs(igr$length) / min(lengthG1,lengthG2)) 
      igrs$log2NormDepth[i] <- round(min(exprG1, exprG2) + weight, digits = 4)
      # normDepth
      exprG1 <- transcCDSs$normDepth[which(rownames(transcCDSs) == igr$startLocusTag)]
      exprG2 <- transcCDSs$normDepth[which(rownames(transcCDSs) == igr$endLocusTag)]
      diff <- abs(exprG1 - exprG2) 
      weight <- diff * (abs(igr$length) / min(lengthG1,lengthG2)) 
      igrs$normDepth[i] <- round(min(exprG1, exprG2) + weight, digits = 4)
    }
  }
  return(igrs)
}

#' Determine cutoff values for a given RNA-seq expression profile.
#'
#' Internal function to estimate the cutoff values to distinguish low expressed RNA-seq data from high expression data.
#' @param data A coverage-depth table.
#' @keywords internal
#' @author Vittorio Fortino
#' qcut()
qcut <- function(data) {
  #set vector for cutoff values
  cutv <- rep(0,0)
  cutv.2 <- rep(0,0)
  for (i in 1:ncol(data)) {
    #specify array and remove 0 counts
    xx <- data[,i]
    #remove 0 counts
    xx <- xx[-which(xx==0)]
    xx <- stats::na.omit(xx)
    #take log2 of data
    log2xx <- log2(xx)
    dlog2 <- data.frame(LogC=log2xx)
    #vector to store Kolmogorov Smirnov distance statistics
    vv <- rep(0,0)
    #select start point
    start <- length(log2xx[log2xx==min(log2xx)])/length(log2xx)
    #set sequence
    s <- seq(round(start,2),0.5,by=0.005)
    #loop through cuts of the data to determine targeted K-S statistic
    for(q in s) {
      #select data greater than a quantile and run Mclust on that data to determine theoretical distribution
      d <- log2xx[which(log2xx>stats::quantile(log2xx,q,na.rm=T))]
      out <- mclust::Mclust(d,G=1)
      ks  <- suppressWarnings(stats::ks.test(d,"pnorm",out$parameter$mean, out$parameter$variance$sigmasq))
      vv  <- c(vv,ks$statistic)
    }
    #determine first left-most local minima
    out <- earth::earth(s,vv,thresh=0.005)
    #save suggested cut
    cutv <- c(cutv,min(out$cuts[out$cuts>0]))
  }
  #send results to outfile
  cutv
}

#' Statistical test to find a putative transcription start (or end) point.
#' 
#' Test the correlation coefficient between a segment of coverage depth and a vector of 100 integers modeling 
#' a simple shape of sharp increases (or decreases) in transcription: x =???[0..0,1..1] (or x =???[1..1,0..0]). 
#' 
#' @param x Segment of coverage depth.
#' @param sharpInc Default value is 1.
#' @param cutCorr Cutoff value for the correlation coefficient. Default values is 0.7.
#' @param cutPVal Cutoff value for the p-value. Default values is 0.0000001.
#' @param flag Indicate weather the correlation test is for a start- or end-point in transcription. Default values is 0 (start-point).
#' @keywords internal
#' @author Vittorio Fortino
#' test.corr()
test.corr <- function(x, sharpInc = 1, cutCorr = 0.7, cutPVal = 0.0000001, flag = 0) { 
  result <- c(test = 0, cc = NA, avg.cov.left = NA, avg.cov.right = NA, diff.min.max = NA)
  sw <- length(x)
  hw <- sw/2
  if(sum(x) > 0) {
    avgCovL <- mean(x[1:hw])
    avgCovR <- mean(x[(hw+1):sw])
    if(flag == 0) 
      testLogExpr <- log2((avgCovR+1)/(avgCovL+1)) >= sharpInc   # IncFWD / DecREV  (___|^^^)
    else 
      testLogExpr <- log2((avgCovL+1)/(avgCovR+1)) >= sharpInc   # IncREV / DecFWD  (^^^|___)
    
    if(testLogExpr == TRUE) {
      if(flag == 0)  
        test.corr <- stats::cor.test(x,rep(c(0,1),c(hw,hw))) # __________|^^^^^^^^^^
      else
        test.corr <- stats::cor.test(x,rep(c(1,0),c(hw,hw))) # ^^^^^^^^^^|__________
      if(!is.na(test.corr$estimate)) {
        if(test.corr$estimate > cutCorr & test.corr$p.value < cutPVal) {
          if(flag == 0) 
            diffMinMax <- abs(min(x[1:hw]) - max(x[(hw+1):sw]))
          else
            diffMinMax <- abs(max(x[1:hw]) - min(x[(hw+1):sw]))	
          result <- c(test = 1, cc = test.corr$estimate, avg.cov.left = avgCovL, avg.cov.right = avgCovR, diff.min.max = diffMinMax)
        }
      }
    }
  }	
  return(result)
}

#' Find strat/end transcription points.
#'
#' Internal function to identify the boundaries of transcriptionally active regions using a sliding window algorithm.
#' The sliding window algorithm uses fixed windows with a length of 100 nt that slides across the coverage-depth data table
#' and finds segments of coverage depth highly and statistically correlated with a vector of 100 integers modeling 
#' a simple shape of sharp increases (or decreases) in transcription: x =???[0..0,1..1] (or x =???[1..1,0..0])). 
#' 
#' With default values, segments having a positive correlation coefficient (exceeding 0.7) and a significant correlation test p-value (<10-7) are selected. 
#' are slected. The vector of the sliding window of 100 integers is a good trade-off between the accuracy of sharp increases/decreases 
#' in transcription and the computational costs of the procedure. P-value 10-7 allows reliable identification of sharp increases/decreases
#' in transcription.  
#' 
#' @param cd A data table containing the coverage depth of an RNA-seq expression profile(s).
#' @param sizeWindow An annotation data table.
#' @param verbose 
#' @keywords internal
#' @author Vittorio Fortino
#' detect.sid.points()
detect.sid.points <- function(cd, sizeWindow, verbose = TRUE) {  
  stopifnot(length(cd$fwd) >= sizeWindow) 
  result <- c(numeric(1),double(4))
  hw <- sizeWindow/2
  nm1 <- sizeWindow-1
  p <- sizeWindow/2
  sg <- length(cd$fwd)
  ## Inc FWD / Dec REV  (___|^^^)
  if(verbose) cat("Identifying ___|^^^ on FWD\n")
  incFWD  <- lapply((sizeWindow+1):sg, FUN = function(i) {test.corr(x=cd$fwd[(i-nm1):i],flag=0)})
  incFWD.pl  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(rep(0,p-i), cd$fwd[1:(i+p)]), flag=0)})
  incFWD.pr  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(cd$fwd[(sg-(p+i-1)):sg], rep(0,p-i)), flag=0)})
  dataIncFWD <- rbind(as.matrix(plyr::ldply(incFWD.pl)), as.matrix(plyr::ldply(incFWD)), as.matrix(plyr::ldply(incFWD.pr)))
  remove(list=c("incFWD","incFWD.pl","incFWD.pr"))
  if(verbose) cat("Identifying ___|^^^ on FWD\n")
  decREV     <- lapply((sizeWindow+1):sg, FUN = function(i) {test.corr(x=cd$rev[(i-nm1):i],flag=0)})
  decREV.pl  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(rep(0,p-i), cd$rev[1:(i+p)]), flag=0)})
  decREV.pr  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(cd$rev[(sg-(p+i-1)):sg], rep(0,p-i)), flag=0)})
  dataDecREV <- rbind(as.matrix(plyr::ldply(decREV.pl)), as.matrix(plyr::ldply(decREV)), as.matrix(plyr::ldply(decREV.pr)))
  remove(list=c("decREV","decREV.pl","decREV.pr"))
  ## Inc REV / Dec FWD  (^^^|___)
  if(verbose) cat("Identifying ^^^|___ on FWD\n")
  incREV     <- lapply((sizeWindow+1):sg, FUN = function(i) {test.corr(x=cd$rev[(i-nm1):i],flag=1)})
  incREV.pl  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(rep(0,p-i), cd$rev[1:(i+p)]), flag=1)})
  incREV.pr  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(cd$rev[(sg-(p+i-1)):sg], rep(0,p-i)), flag=1)})
  dataIncREV <- rbind(as.matrix(plyr::ldply(incREV.pl)), as.matrix(plyr::ldply(incREV)), as.matrix(plyr::ldply(incREV.pr)))
  remove(list=c("incREV","incREV.pl","incREV.pr"))
  if(verbose) cat("Identifying ^^^|___ on FWD\n")
  decFWD     <- lapply((sizeWindow+1):sg, FUN = function(i) {test.corr(x=cd$fwd[(i-nm1):i],flag=1)})
  decFWD.pl  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(rep(0,p-i), cd$fwd[1:(i+p)]), flag=1)})
  decFWD.pr  <- lapply(1:(sizeWindow/2), FUN = function(i) {test.corr(x = c(cd$fwd[(sg-(p+i-1)):sg], rep(0,p-i)), flag=1)})
  dataDecFWD <- rbind(as.matrix(plyr::ldply(decFWD.pl)), as.matrix(plyr::ldply(decFWD)), as.matrix(plyr::ldply(decFWD.pr)))
  remove(list=c("decFWD","decFWD.pl","decFWD.pr"))
  return(list(fwd.incs = as.data.frame(dataIncFWD), 
              fwd.decs = as.data.frame(dataDecFWD),
              rev.incs = as.data.frame(dataIncREV), 
              rev.decs = as.data.frame(dataDecREV)))
}

#' Determine operon start-points (OSPs).
#'
#' Internal function to estimate the gene-level expression values using the RPKM method.
#' @param fwd.sh.incs Data table containing information on the sharp increases in transcription found on the forward strand. See \code{detect.sid.points}.
#' @param rev.sh.incs Data table containing information on the sharp increases in transcription found on the reverse strand. See \code{detect.sid.points}.
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. 
#' @param igrs.p Data table containing generic information of the intergenic regions on the forward strand. See \code{get.intergenic.regions}.
#' @param igrs.n Data table containing generic information of the intergenic regions on the reverse strand. See \code{get.intergenic.regions}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param borders A vector.
#' @param max.start.transc Maximum acacepted start transcription level.
#' @param minExprCDS Minimum expression level for the coding sequence regions (CDSs). 
#' @param verbose 
#' @keywords internal
#' @author Vittorio Fortino
#' get.operon.start.points() 
get.operon.start.points <- function(fwd.sh.incs, rev.sh.incs, genes.and.operons, igrs.p, igrs.n, 
                                    transcCDSs, borders = c(100,10), max.start.transc = c(0.1,0.1),
                                    minExprCDS = 0.1, verbose = TRUE, ...) { 
  
  if(verbose) cat("\n   Identifying confirmed and putative operon start/end points\n")
  cnt.match <- 0
  name.cds  <- vector(mode = "character")
  start.tr  <- vector(mode = "integer")
  end.tr    <- vector(mode = "integer")
  expr.cds  <- vector(mode = "numeric")
  start.cds <- vector(mode = "integer")
  strand    <- vector(mode = "character")
  cc        <- vector(mode = "numeric")
  point     <- vector(mode = "numeric")
  operon.id <- vector(mode = "numeric")
  count <- 0
  for(i in 1:nrow(igrs.p)) {
    if(igrs.p$type[i] == "O") next
    start <- igrs.p$start[i] 
    end   <- igrs.p$end[i] 
    if((end-start) > borders[1])
      dat <- fwd.sh.incs[(end-borders[1]):(end+borders[2]),]
    else
      dat <- fwd.sh.incs[start:(end+borders[2]),]
    one_matches <- which(dat$test == 1)
    if(length(one_matches) > 0) {
      idc <- which(dat$avg.cov.left[one_matches] == min(dat$avg.cov.left[one_matches]))
      maxs <- one_matches[idc]
      idw <- which(dat$avg.cov.right[maxs] == max(dat$avg.cov.right[maxs]))
      match <- maxs[idw[1]] 
      x <- dat$avg.cov.left[match]
      y <- dat$avg.cov.right[match]
      test_log2FC <- FALSE
      test_Ygt1   <- FALSE
      if(x != 0) test_log2FC <- log2(y/x) > 1 # sharpSize
      else test_Ygt1 <- y > 1                 # to avoid transc regions with 0000011111
      test_start_expr <- log2(x) <= max.start.transc[1] 	
      test_min_expr   <- transcCDSs[igrs.p$endLocusTag[i],]$expr > minExprCDS
      id <- which(rownames(genes.and.operons) == igrs.p$endLocusTag[i])
      if(length(id) > 0) operon <- genes.and.operons$operonID[id]
      else operon <- NA
      ###
      if((test_log2FC == TRUE | test_Ygt1 == TRUE) & test_min_expr & !test_start_expr)	count <- count + 1
      if((test_log2FC == TRUE | test_Ygt1 == TRUE) & test_start_expr & test_min_expr) {
        cnt.match            <- cnt.match + 1
        name.cds[cnt.match]  <- igrs.p$endLocusTag[i]
        start.cds[cnt.match] <- end
        operon.id[cnt.match] <- operon
        strand[cnt.match]  	 <- "+"
        expr.cds[cnt.match]  <- transcCDSs[igrs.p$endLocusTag[i],]$expr
        start.tr[cnt.match]  <- x
        end.tr[cnt.match]    <- y
        cc[cnt.match]        <- dat$cc[match]
        if((end-start) > borders[1])
          point[cnt.match] <- ((end-borders[1]) + match) - 1
        else
          point[cnt.match] <- (start + match) - 1
      }
    } 
  }
  for(i in 1:nrow(igrs.n)) {
    if(igrs.n$type[i] == "O") next
    start <- igrs.n$start[i] 
    end   <- igrs.n$end[i] 
    if((end-start) > borders[1])
      dat   <- rev.sh.incs[(start-borders[2]):(start+borders[1]),]
    else
      dat   <- rev.sh.incs[start:end,]
    one_matches <- which(dat$test == 1)
    if(length(one_matches) > 0) {
      idc <- which(dat$avg.cov.right[one_matches] == min(dat$avg.cov.right[one_matches]))
      maxs <- one_matches[idc]
      idw <- which(dat$avg.cov.left[maxs] == max(dat$avg.cov.left[maxs]))
      match <- maxs[idw[1]]
      x <- dat$avg.cov.left[match]
      y <- dat$avg.cov.right[match]
      test_log2FC <- FALSE
      test_Xgt1 <- FALSE
      test_max_start_transc <- FALSE
      if(y != 0) test_log2FC <- log2(x/y) > 1
      else test_Xgt1 <- x > 1
      test_start_expr <- log2(y) <= max.start.transc[2]
      test_min_expr   <- transcCDSs[igrs.n$startLocusTag[i],]$expr > minExprCDS
      id <- which(rownames(genes.and.operons) == igrs.n$startLocusTag[i])
      if(length(id) > 0) operon <- genes.and.operons$operonID[id]
      else operon <- NA
      if((test_log2FC == TRUE | test_Xgt1 == TRUE) & test_min_expr & !test_start_expr)	count <- count + 1
      if((test_log2FC == TRUE | test_Xgt1 == TRUE) & test_start_expr & test_min_expr) {
        cnt.match            <- cnt.match + 1
        name.cds[cnt.match]  <- igrs.n$startLocusTag[i]
        start.cds[cnt.match] <- start
        operon.id[cnt.match] <- operon
        strand[cnt.match]    <- "-"
        expr.cds[cnt.match]  <- transcCDSs[igrs.n$startLocusTag[i],]$expr
        end.tr[cnt.match]    <- x
        start.tr[cnt.match]  <- y
        cc[cnt.match]  	     <- dat$cc[match]
        if((end-start) > borders[1])
          point[cnt.match] <- (((start-borders[2]) + match) - 1)
        else
          point[cnt.match] <- ((start + match) - 1)
      } 
    } 
  } 
  pOSSs <- data.frame(name.cds = name.cds, lp = start.tr, hp = end.tr, 
                      expr.cds = expr.cds, start.cds = start.cds, 
                      strand = strand, cc = cc, point = point, 
                      operon.id = operon.id, stringsAsFactors = FALSE)
  
  if(verbose) cat("    - Result:  ", nrow(pOSSs),
                  " identified operon start points(pos=", 
                  length(which(pOSSs$strand == "+")),
                  ", neg=", length(which(pOSSs$strand == "-")), ") \n", sep="")
  
  # To select the best pOSS for each operon
  unique.cds <- unique(pOSSs$name.cds)
  if(length(unique.cds) < length(pOSSs$name.cds)) {
    if(verbose) cat("    - Eliminating duplicated rows  \n")
    best.poss.cds <- vector("integer", length = length(unique.cds))
    for(i in 1: length(unique.cds)) {
      list_poss <- which(pOSSs$name.cds == unique.cds[i])
      strand <- pOSSs$strand[list_poss[1]]
      id_min <- which(pOSSs[list_poss,]$lp == min(pOSSs[list_poss,]$lp))
      if(length(id_min) > 1) {
        id_max <- which(pOSSs[list_poss[id_min],]$hp == max(pOSSs[list_poss[id_min],]$hp))
        best.poss.cds[i] <- list_poss[id_min[id_max[1]]]
      }
      else {
        best.poss.cds[i] <- list_poss[id_min[1]]
      }
    }
    pOSSs <- pOSSs[best.poss.cds,]
    rownames(pOSSs) <- pOSSs[,1]
    pOSSs <- pOSSs[,-1]
  }
  rownames(pOSSs) <- pOSSs[,1]
  pOSSs <- pOSSs[,-1]
  if(verbose) cat("    - Number of unique POSSs = ", nrow(pOSSs), "\n", sep="")
  return(list(POSSs = pOSSs[!is.na(pOSSs$operon.id),], newPOSSs = pOSSs[is.na(pOSSs$operon.id),]))
}

#' Determine operon end-points (OEPs).
#'
#' Internal function to estimate the gene-level expression values using the RPKM method.
#' @param fwd.sh.incs Data table containing information on the sharp decreases in transcription found on the forward strand. See \code{detect.sid.points}.
#' @param rev.sh.incs Data table containing information on the sharp decreases in transcription found on the reverse strand. See \code{detect.sid.points}.
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. 
#' @param igrs.p Data table containing generic information of the intergenic regions on the forward strand. See \code{get.intergenic.regions}.
#' @param igrs.n Data table containing generic information of the intergenic regions on the reverse strand. See \code{get.intergenic.regions}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param borders A numeric vector. 
#' @param max.start.transc Maximum log2.
#' @param minExprCDS Minimum expression level for the coding sequence regions (CDSs). Default values is 0.1.
#' @param verbose 
#' @keywords internal
#' @author Vittorio Fortino
#' get.operon.start.points()
get.operon.end.points <- function(fwd.sh.decs, rev.sh.decs, genes.and.operons, igrs.p, igrs.n, 
                                  transcCDSs, borders = c(10,100), max.end.transc = c(0.1,0.1),
                                  minExprCDS = 0.1, verbose = TRUE,  ...) { 
  if(verbose) cat("\n   Identifying confirmed and putative operon start/end points\n")
  cnt.match <- 0
  name.cds  <- vector(mode = "character")
  start.tr  <- vector(mode = "integer")
  end.tr    <- vector(mode = "integer")
  expr.cds  <- vector(mode = "numeric")
  end.cds   <- vector(mode = "integer")
  strand    <- vector(mode = "character")
  cc        <- vector(mode = "numeric")
  point     <- vector(mode = "numeric")
  operon.id <- vector(mode = "numeric")
  count <- 0
  for(i in 1:nrow(igrs.p)) {
    if(igrs.p$type[i] == "O") next
    start <- igrs.p$start[i] 
    end   <- igrs.p$end[i] 
    if((end-start) > borders[2])
      dat   <- fwd.sh.decs[(start-borders[1]):(start+borders[2]),]
    else
      dat   <- fwd.sh.decs[start:end,]
    one_matches <- which(dat$test == 1) 
    if(length(one_matches) > 0) {
      idc <- which(dat$avg.cov.left[one_matches] == max(dat$avg.cov.left[one_matches]))
      maxs <- one_matches[idc]
      idw <- which(dat$avg.cov.right[maxs] == min(dat$avg.cov.right[maxs]))
      match <- maxs[idw[1]] 
      x <- dat$avg.cov.left[match]
      y <- dat$avg.cov.right[match]
      test_log2FC <- FALSE
      test_Xgt1 <- FALSE
      if(y != 0) test_log2FC <- log2(x/y) > 1
      else test_Xgt1 <- x > 1
      test_end_expr <- log2(y) <= max.end.transc[1] 	
      test_min_expr   <- transcCDSs[igrs.p$startLocusTag[i],]$expr > minExprCDS
      id <- which(rownames(genes.and.operons) == igrs.p$startLocusTag[i])
      if(length(id) > 0) operon <- genes.and.operons$operonID[id]
      else operon <- NA
      ###
      if((test_log2FC == TRUE | test_Xgt1 == TRUE) & !test_end_expr & test_min_expr)	count <- count + 1
      if((test_log2FC == TRUE | test_Xgt1 == TRUE) & test_end_expr & test_min_expr) {
        cnt.match            <- cnt.match + 1
        name.cds[cnt.match]  <- igrs.p$startLocusTag[i]
        end.cds[cnt.match]   <- start
        operon.id[cnt.match] <- operon
        strand[cnt.match]    <- "+"
        expr.cds[cnt.match]  <- transcCDSs[igrs.p$startLocusTag[i],]$expr
        start.tr[cnt.match]  <- x
        end.tr[cnt.match]    <- y
        cc[cnt.match]        <- dat$cc[match]
        if((end-start) > borders[1])
          point[cnt.match] <- ((start-borders[1]) + match) - 1
        else
          point[cnt.match] <- (start + match) - 1
      } 
    }
  }
  for(i in 1:nrow(igrs.n)) {
    if(igrs.n$type[i] == "O") next
    start <- igrs.n$start[i] 
    end   <- igrs.n$end[i] 
    if((end-start) > borders[2])
      dat   <- rev.sh.decs[(end-borders[2]):(end+borders[1]),]
    else
      dat   <- rev.sh.decs[start:end,]
    one_matches <- which(dat$test == 1)
    if(length(one_matches) > 0) {
      idc <- which(dat$avg.cov.right[one_matches] == max(dat$avg.cov.right[one_matches]))
      maxs <- one_matches[idc]
      idw <- which(dat$avg.cov.left[maxs] == min(dat$avg.cov.left[maxs]))
      match <- maxs[idw[1]]
      x <- dat$avg.cov.left[match]
      y <- dat$avg.cov.right[match]
      test_log2FC <- FALSE
      test_Ygt1 <- FALSE
      if(x != 0) test_log2FC <- log2(y/x) > 1
      else test_Ygt1 <- y > 1	
      test_end_expr <- log2(x) <= max.end.transc[2]
      test_min_expr   <- transcCDSs[igrs.n$endLocusTag[i],]$expr > minExprCDS
      id <- which(rownames(genes.and.operons) == igrs.n$endLocusTag[i])
      if(length(id) > 0) operon <- genes.and.operons$operonID[id]
      else operon <- NA
      if((test_log2FC == TRUE | test_Ygt1 == TRUE) & test_min_expr & !test_end_expr)	count <- count + 1
      if((test_log2FC == TRUE | test_Ygt1 == TRUE) & test_min_expr & test_end_expr) {
        cnt.match            <- cnt.match + 1
        name.cds[cnt.match]  <- igrs.n$endLocusTag[i]
        end.cds[cnt.match]   <- end
        operon.id[cnt.match] <- operon
        strand[cnt.match]    <- "-"
        expr.cds[cnt.match]  <- transcCDSs[igrs.n$endLocusTag[i],]$expr
        end.tr[cnt.match]    <- x
        start.tr[cnt.match]  <- y
        cc[cnt.match]        <- dat$cc[match]
        if((end-start) > borders[1])
          point[cnt.match] <- (((end-borders[2]) + match) - 1)
        else
          point[cnt.match] <- ((start + match) - 1)
      }
    } 
  } 
  pOESs <- data.frame(name.cds = name.cds, hp = start.tr, lp = end.tr, 
                      expr.cds = expr.cds, end.cds = end.cds, 
                      strand = strand, cc = cc, point = point, 
                      operon.id = operon.id, stringsAsFactors = FALSE)
  
  if(verbose) cat("    - Result:  ", nrow(pOESs),
                  " identified operon end points(pos=", 
                  length(which(pOESs$strand == "+")),
                  ", neg=", length(which(pOESs$strand == "-")), ") \n", sep="")
  
  # To select the best pOES for each operon
  unique.cds <- unique(pOESs$name.cds)
  if(length(unique.cds) < length(pOESs$name.cds)) {
    if(verbose) cat("    - Eliminating duplicated rows  \n")
    best.poes.cds <- vector("integer", length = length(unique.cds))
    for(i in 1: length(unique.cds)) {
      list_poes <- which(pOESs$name.cds == unique.cds[i])
      strand <- pOESs$strand[list_poes[1]]
      id_min <- which(pOESs[list_poes,]$lp == min(pOESs[list_poes,]$lp))
      if(length(id_min) > 1) {
        id_max <- which(pOESs[list_poes[id_min],]$hp == max(pOESs[list_poes[id_min],]$hp))
        best.poes.cds[i] <- list_poes[id_min[id_max[1]]]
      }
      else {
        best.poes.cds[i] <- list_poes[id_min[1]]
      }
    }
    pOESs <- pOESs[best.poes.cds,]
    rownames(pOESs) <- pOESs[,1]
    pOESs <- pOESs[,-1]
  }
  rownames(pOESs) <- pOESs[,1]
  pOESs <- pOESs[,-1]
  if(verbose) cat("    - Number of unique POESs = ", nrow(pOESs), "\n", sep="")
  return(list(POESs = pOESs[!is.na(pOESs$operon.id),], newPOESs = pOESs[is.na(pOESs$operon.id),]))
}

#' Compile a set of coinfirmed operons. 
#' 
#' Test the correlation coefficient beween a segment of coverage depth and a vector of 100 integers modeling 
#' a simple shape of sharp increases (or decreases) in transcription: x =[0..0,1..1] (or x =[1..1,0..0]). 
#' 
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param transcIGRs.pos Transcription levels for the intergenic regions (forard strand). See \code{comp.igr.transc.levels}.
#' @param transcIGRs.neg Transcription levels for the intergenic regions (reverse strand). See \code{comp.igr.transc.levels}.
#' @param POSSs Data table representing a set of putative operon start-points.
#' @param POESs Data table representing a set of putative operon end-points.
#' @param minExprCDS Minimum expression level for the coding sequence regions (CDSs). Default values is 0.1.
#' @param minExprIGR Minimum expression level for the intergenic regions (IGRs). Default values is 0.25.
#' @param max.start.transc Cutoff values for the start transcription points. Default values is 0.1.
#' @param max.end.transc Cutoff values for the end transcription points. Default values is 0.1.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' compile.confirmed.operons()
compile.confirmed.operons <- function(genes.and.operons, transcCDSs, transcIGRs.pos, transcIGRs.neg, 
                                      POSSs, POESs, minExprCDS = 0.1, minExprIGR = 0.25, 
                                      max.start.transc = c(0.1,0.1), max.end.transc = c(0.1,0.1), 
                                      verbose = TRUE, ...) {
  cnt.OPs <- 0
  operon.id  <- vector(mode = "character") 
  operon.ref <- vector(mode = "character") 
  start.cds  <- vector(mode = "character") 
  strand     <- vector(mode = "character") 
  expr       <- vector(mode = "character")
  length     <- vector(mode = "integer")
  idgs       <- vector(mode = "integer") 
  i <- 1
  while(i <= nrow(genes.and.operons)) {
    a.gene <- genes.and.operons[i,]
    if(!is.na(a.gene$operonID) & a.gene$strand == '+') {	
      opRef <- a.gene$operonID
      gene.s <- genes.and.operons[which(genes.and.operons$operonID == opRef),]
      numGenesInOp <- nrow(gene.s) 
      exprVals     <- vector(mode = "numeric", length = numGenesInOp)
      possTags     <- vector(mode = "integer", length = numGenesInOp)
      poesTags     <- vector(mode = "integer", length = numGenesInOp)
      test_IGR     <- FALSE
      test_len_IGR <- FALSE
      test_CDS     <- FALSE
      test_X       <- FALSE
      test_Y       <- FALSE
      startNewOp   <- -1
      endNewOp     <- -1
      for(j in 1:numGenesInOp) {
        possTags[j] <- length(which(rownames(POSSs) == rownames(gene.s)[j]))
        poesTags[j] <- length(which(rownames(POESs) == rownames(gene.s)[j]))
        exprVals[j] <- transcCDSs$expr[which(rownames(transcCDSs) == rownames(gene.s)[j])]
        # 5 tests to check the end of the linkage process 
        if(j != numGenesInOp) { 
          id_IGR   <- which(transcIGRs.pos$startLocusTag == rownames(gene.s)[j])
          test_IGR <- transcIGRs.pos$expr[id_IGR] < minExprIGR
          id_CDS   <- which(rownames(transcCDSs) == rownames(gene.s)[j+1])
          test_CDS <- transcCDSs$expr[id_CDS] < minExprCDS
        }
        if(possTags[j] == 1 & startNewOp != -1) {
          id_poss <- which(rownames(POSSs) == rownames(gene.s)[j])
          test_X <- log2(POSSs$lp[id_poss]) <= max.start.transc[1]
        }
        if(poesTags[j] == 1 & startNewOp != -1) {
          id_poes <- which(rownames(POESs) == rownames(gene.s)[j])
          test_Y <- log2(POESs$lp[id_poes]) <= max.end.transc[1]
        }
        test_startOp <- (possTags[j] == 1 & startNewOp == -1)
        if(test_startOp == TRUE) { 
          startNewOp <- j
          endNewOp <- j
        } 
        if(startNewOp != -1 & (test_Y == TRUE | test_X == TRUE | test_IGR == TRUE | test_CDS == TRUE | j == numGenesInOp)) { 
          if(endNewOp > startNewOp & startNewOp != -1) { 
            cnt.OPs <- cnt.OPs + 1
            operon.id[cnt.OPs]  <- paste("Operon#",cnt.OPs,sep='')
            operon.ref[cnt.OPs] <- opRef
            start.cds[cnt.OPs]  <- rownames(gene.s)[startNewOp]
            strand[cnt.OPs]     <- a.gene$strand 
            length[cnt.OPs]     <- endNewOp - startNewOp + 1
            expr[cnt.OPs]       <- paste(exprVals[startNewOp:endNewOp],collapse='-')
            idgs[cnt.OPs]       <- paste(rownames(gene.s[startNewOp:endNewOp,]),collapse='-')
            startNewOp <- -1
            endNewOp   <- -1
          }
          else {  
            startNewOp <- -1
            endNewOp   <- -1
          }
        }
        else if(startNewOp != -1) endNewOp <- endNewOp + 1
      }
      i <- (i + numGenesInOp) - 1
    } 
    if(!is.na(a.gene$operonID) == TRUE & a.gene$strand == '-') {	
      opRef <- a.gene$operonID
      gene.s <- genes.and.operons[which(genes.and.operons$operonID == opRef),]
      numGenesInOp <- nrow(gene.s) 
      exprVals     <- vector(mode = "numeric", length = numGenesInOp)
      possTags     <- vector(mode = "integer", length = numGenesInOp)
      poesTags     <- vector(mode = "integer", length = numGenesInOp)
      test_IGR     <- FALSE
      test_len_IGR <- FALSE
      test_CDS     <- FALSE
      test_X       <- FALSE
      test_Y       <- FALSE
      startNewOp   <- -1
      endNewOp     <- -1
      for(j in numGenesInOp:1) {
        possTags[j] <- length(which(rownames(POSSs) == rownames(gene.s)[j]))
        poesTags[j] <- length(which(rownames(POESs) == rownames(gene.s)[j]))
        exprVals[j] <- transcCDSs$expr[which(rownames(transcCDSs) == rownames(gene.s)[j])]
        if(j != 1) { 
          id_IGR   <- which(transcIGRs.neg$endLocusTag == rownames(gene.s)[j])
          test_IGR <- transcIGRs.neg$expr[id_IGR] < minExprIGR
          id_CDS <- which(rownames(transcCDSs) == rownames(gene.s)[j-1])
          test_CDS <- transcCDSs$expr[id_CDS] < minExprCDS
        }
        if(possTags[j] == 1 & startNewOp != -1) {
          id_poss <- which(rownames(POSSs) == rownames(gene.s)[j])
          test_Y <- log2(POSSs$lp[id_poss]) <= max.start.transc[2]
        }
        if(poesTags[j] == 1 & startNewOp != -1) {
          id_poes <- which(rownames(POESs) == rownames(gene.s)[j])
          test_X <- log2(POESs$lp[id_poes]) <= max.end.transc[2]
        }
        test_startOp <- (possTags[j] == 1 & startNewOp == -1)
        if(test_startOp == TRUE) { 
          startNewOp <- j
          endNewOp <- j
        }
        if(startNewOp != -1 & (test_Y == TRUE | test_X == TRUE | test_IGR == TRUE | test_CDS == TRUE | j == 1)) { 
          if(startNewOp > endNewOp & startNewOp != -1) { 
            cnt.OPs <- cnt.OPs + 1
            operon.id[cnt.OPs]  <- paste("Operon#",cnt.OPs,sep='')
            operon.ref[cnt.OPs] <- opRef
            start.cds[cnt.OPs]  <- rownames(gene.s)[startNewOp]
            strand[cnt.OPs]     <- a.gene$strand 
            length[cnt.OPs]     <- startNewOp - endNewOp + 1
            expr[cnt.OPs]       <- paste(exprVals[endNewOp:startNewOp],collapse='-')
            idgs[cnt.OPs]       <- paste(rownames(gene.s[endNewOp:startNewOp,]),collapse='-')
            startNewOp <- -1
            endNewOp <- -1
          }
          else {
            startNewOp <- -1
            endNewOp   <- -1
          }
        }
        else { 
          if(startNewOp != -1)  endNewOp <- endNewOp - 1
        }
      }
      i <- (i + numGenesInOp) - 1
    }  
    i <- i + 1
  }	
  cops <- data.frame(opID = operon.id, operon.ref = operon.ref, start.cds = start.cds,
                     strand = strand, length = length, idgs = idgs, expr = expr, 
                     row.names = 1, stringsAsFactors = FALSE)
  id.dups <- which(duplicated(cops$idgs) == TRUE)
  if(length(id.dups) > 0)	cops <- cops[-id.dups,]
  if(verbose) cat("     - Number of confirmed operons:", nrow(cops), "\n")
  return(cops)
}

#' Define a set of OPs which is used to train the operon classifier.
#' 
#' Build a data table containing confirmed operon pairs (OPs) and the feature values. 
#' 
#' @param listOperons List of confirmed operons. See \code{compile.confirmed.operons}.
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param transcIGRs.pos Transcription levels for the intergenic regions (forard strand). See \code{comp.igr.transc.levels}.
#' @param transcIGRs.neg Transcription levels for the intergenic regions (reverse strand). See \code{comp.igr.transc.levels}.
#' @param wseq Sequence genome.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' select.ops()
select.ops <- function(listOperons, genes.and.operons, 
                       transcCDSs, transcIGRs.pos, transcIGRs.neg, wseq, verbose=TRUE) {
  cnt.OPs <- 1
  opRef		<- vector(mode = "character")
  nameG1      <- vector(mode = "character")
  nameG2      <- vector(mode = "character")
  diffExpr    <- vector(mode = "numeric")    
  lengthIGR   <- vector(mode = "integer")
  cuScore     <- vector(mode = "numeric") 
  exprIGR     <- vector(mode = "numeric")
  class       <- vector(mode = "character")
  for(j in 1:length(listOperons[,1])) {
    operon <- listOperons[j,]
    gene.s <- unique(strsplit(operon$idgs, "-")[[1]])
    if(operon$strand == '+') {
      for(i in 1: (length(gene.s)-1)) {
        gene1 <- genes.and.operons[gene.s[i],]
        gene2 <- genes.and.operons[gene.s[i+1],]
        geneSeq1 <- wseq[gene1$start:gene1$end]
        geneSeq2 <- wseq[gene2$start:gene2$end]
        rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
        rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
        cuScore[cnt.OPs] <- sum(rscu1 * rscu2, na.rm = TRUE)
        id_IGR <- which(transcIGRs.pos$startLocusTag == rownames(gene1))
        id_CDS <- which(rownames(transcCDSs) == rownames(gene1))
        lengthIGR[cnt.OPs] <- transcIGRs.pos$length[id_IGR]
        exprIGR[cnt.OPs]   <- transcIGRs.pos$expr[id_IGR]
        diffExpr[cnt.OPs]  <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]), digits=4)
        opRef[cnt.OPs]	   <- operon$operon.ref
        nameG1[cnt.OPs]	  <- rownames(gene1)
        nameG2[cnt.OPs]	  <- rownames(gene2)
        class[cnt.OPs]	  <- "OP"
        cnt.OPs <- cnt.OPs + 1
      }
    }
    else { 
      for(i in (length(gene.s)-1):1) {
        gene1 <- genes.and.operons[gene.s[i],]
        gene2 <- genes.and.operons[gene.s[i+1],]
        geneSeq1 <- wseq[gene1$start:gene1$end]
        geneSeq2 <- wseq[gene2$start:gene2$end]
        rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
        rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
        cuScore[cnt.OPs] <- sum(rscu1 * rscu2, na.rm = TRUE)
        id_IGR  <- which(transcIGRs.neg$startLocusTag == rownames(gene1))
        id_CDS  <- which(rownames(transcCDSs) == rownames(gene1))
        lengthIGR[cnt.OPs] <- transcIGRs.neg$length[id_IGR]
        exprIGR[cnt.OPs]   <- transcIGRs.neg$expr[id_IGR]
        diffExpr[cnt.OPs]  <- abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1])
        opRef[cnt.OPs]	  <- operon$operon.ref
        nameG1[cnt.OPs]	  <- rownames(gene1)
        nameG2[cnt.OPs]	  <- rownames(gene2)
        class[cnt.OPs]	  <- "OP"
        cnt.OPs <- cnt.OPs + 1
      }
    }
  }
  if(verbose) cat("    - nameG1:", length(nameG1), ", nameG2:", length(nameG2),
                  ", diffExpr:", length(diffExpr), ", lengthIGR:", length(lengthIGR), 
                  ", cuScore:", length(cuScore), ", opRef:", length(opRef), "\n", sep="")
  
  return(data.frame(G1 = nameG1, G2 = nameG2, 
                    diffExpr = diffExpr, exprIGR = exprIGR, lengthIGR = lengthIGR, cuScore = cuScore,
                    opRef = opRef, class = class, stringsAsFactors = FALSE))
}

#' Define a set of NOPs which is used to train the operon classifier.
#' 
#' Build a data table containing confirmed non-operon pairs (NOPs) and the feature values. 
#' 
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param OPs Data table including the confirmed operon pairs (OPs). See \code{select.ops }.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param transcIGRs.pos Transcription levels for the intergenic regions (forard strand). See \code{comp.igr.transc.levels}.
#' @param transcIGRs.neg Transcription levels for the intergenic regions (reverse strand). See \code{comp.igr.transc.levels}.
#' @param POSSs Data table representing a set of putative operon start-points.
#' @param POESs Data table representing a set of putative operon end-points.
#' @param wseq Sequence genome.
#' @param max.start.transc Cutoff values for the start transcription points. Default values is 0.1.
#' @param max.end.transc Cutoff values for the end transcription points. Default values is 0.1.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' select.nops()
select.nops <- function(genes.and.operons, OPs, transcCDSs, 
                        transcIGRs.pos, transcIGRs.neg, POSSs, POESs, wseq, 
                        max.start.transc = c(0.1,0.1), max.end.transc = c(0.1,0.1), 
                        verbose = TRUE, ...) {
  cnt.NOPs <- 1
  nameG1      <- vector(mode = "character")
  nameG2      <- vector(mode = "character")
  refOp1		<- vector(mode = "character")
  refOp2		<- vector(mode = "character")
  diffExpr    <- vector(mode = "numeric")    
  lengthIGR   <- vector(mode = "integer")
  cuScore     <- vector(mode = "numeric")  
  exprIGR     <- vector(mode = "numeric")
  class       <- vector(mode = "character")
  for(i in 1:nrow(transcIGRs.pos)) {
    gene1 <- transcIGRs.pos$startLocusTag[i]
    gene2 <- transcIGRs.pos$endLocusTag[i]
    g1 <- genes.and.operons[which(rownames(genes.and.operons) == gene1),]
    g2 <- genes.and.operons[which(rownames(genes.and.operons) == gene2),]	
    id1 <- which(OPs$G1 == gene1)
    id2 <- which(OPs$G2 == gene2)
    if(length(id1) > 0 & length(id2) > 0)
      if(id1 == id2) next
    test <- FALSE
    if(length(which(rownames(POSSs) == gene2)) > 0) {
      id_poss <- which(rownames(POSSs) == gene2)
      if(log2(POSSs$lp[id_poss]) <= max.start.transc[1]) test <- TRUE
    }
    if(length(which(rownames(POESs) == gene1)) > 0) {
      id_poes <- which(rownames(POESs) == gene1)
      if(log2(POESs$lp[id_poes]) <= max.end.transc[1]) test <- TRUE
    }
    if(test == TRUE) {
      geneSeq1 = wseq[g1$start:g1$end]
      geneSeq2 = wseq[g2$start:g2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt.NOPs] <- sum(rscu1 * rscu2, na.rm = TRUE)
      id_IGR   <- which(transcIGRs.pos$startLocusTag == gene1)
      lengthIGR[cnt.NOPs] <- transcIGRs.pos$length[id_IGR]
      exprIGR[cnt.NOPs]   <- transcIGRs.pos$expr[id_IGR]
      id_CDS  <- which(rownames(transcCDSs) == gene1)
      diffExpr[cnt.NOPs] <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt.NOPs]   <- gene1
      nameG2[cnt.NOPs]   <- gene2
      refOp1[cnt.NOPs]   <- g1$operonID
      refOp2[cnt.NOPs]   <- g2$operonID
      class[cnt.NOPs]	   <- "NOP"
      cnt.NOPs <- cnt.NOPs + 1
    }
  }
  for(i in 1:nrow(transcIGRs.neg)) {
    gene1 <- transcIGRs.neg$startLocusTag[i]
    gene2 <- transcIGRs.neg$endLocusTag[i]
    g1 <- genes.and.operons[which(rownames(genes.and.operons) == gene1),]
    g2 <- genes.and.operons[which(rownames(genes.and.operons) == gene2),]
    id1 <- which(OPs$G1 == gene1)
    id2 <- which(OPs$G2 == gene2)
    if(length(id1) > 0 & length(id2) > 0)
      if(id1 == id2) next
    test <- FALSE
    if(length(which(rownames(POSSs) == gene1)) > 0) {
      id_poss <- which(rownames(POSSs) == gene1)
      if(log2(POSSs$lp[id_poss]) <= max.start.transc[2]) test <- TRUE
    }
    if(length(which(rownames(POESs) == gene2)) > 0) {
      id_poes <- which(rownames(POESs) == gene2)
      if(log2(POESs$lp[id_poes]) <= max.end.transc[2]) test <- TRUE
    }
    if(test == TRUE) {
      geneSeq1 = wseq[g1$start:g1$end]
      geneSeq2 = wseq[g2$start:g2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt.NOPs] <- sum(rscu1 * rscu2, na.rm = TRUE)
      id_IGR   <- which(transcIGRs.neg$startLocusTag == gene1)
      lengthIGR[cnt.NOPs] <- transcIGRs.neg$length[id_IGR]
      exprIGR[cnt.NOPs]   <- transcIGRs.neg$expr[id_IGR]
      id_CDS  <- which(rownames(transcCDSs) == gene1)
      diffExpr[cnt.NOPs] <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt.NOPs]   <- gene1
      nameG2[cnt.NOPs]   <- gene2
      refOp1[cnt.NOPs]   <- g1$operonID
      refOp2[cnt.NOPs]   <- g2$operonID
      class[cnt.NOPs]	   <- "NOP"
      cnt.NOPs <- cnt.NOPs + 1
    }
  }	
  if(verbose)	cat("    - nameG1:", length(nameG1), ", nameG2:", length(nameG2),
                  ", diffExpr:", length(diffExpr), ", lengthIGR:", length(lengthIGR),
                  ", cuScore:",length(cuScore),"\n", sep="")
  
  return(data.frame(G1 = nameG1, 
                    G2 = nameG2, 
                    diffExpr = diffExpr,  
                    exprIGR = exprIGR,  
                    lengthIGR = lengthIGR,
                    cuScore = cuScore,
                    refOp1 = refOp1,
                    refOp2 = refOp2,
                    class = class,
                    stringsAsFactors = FALSE))
}

#' Define a set of gene pairs POPs with an "operon status" to (re)define represent used to train the operon classifier.
#' 
#' Internal function to build a data table containing gene pairs and the correspoding transcriptomic/genomic feature values. 
#' The operon classifier trained/validated on the OPs and NOPs is used to predict the operon status of these gene pairs.
#' 
#' @param OPs Data table including the confirmed operon pairs (OPs). See \code{select.ops}.
#' @param OPs Data table including the confirmed non-operon pairs (NOPs). See \code{select.nops}.
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param transcIGRs.pos Transcription levels for the intergenic regions (forard strand). See \code{comp.igr.transc.levels}.
#' @param transcIGRs.neg Transcription levels for the intergenic regions (reverse strand). See \code{comp.igr.transc.levels}.
#' @param minExprCDS Minimum expression level for the coding sequence regions (CDSs). Default values is 0.1.
#' @param minExprIGR Minimum expression level for the intergenic regions (IGRs). Default values is 0.25.
#' @param maxLenIGR Maximum length for the intergenic regions (IGRs). Default values is 150.
#' @param POSSs Data table representing a set of putative operon start-points.
#' @param POESs Data table representing a set of putative operon end-points.
#' @param max.start.transc Cutoff values for the start transcription points. Default values is 0.1.
#' @param max.end.transc Cutoff values for the end transcription points. Default values is 0.1.
#' @param wseq Sequence genome.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' select.pops()
select.pops  <- function(OPs, NOPs, genes.and.operons, transcCDSs, transcIGRs.pos, transcIGRs.neg, 
                         minExprCDS = 0.1, minExprIGR = 0.25, maxLenIGR = 150, 
                         POSSs, POESs, max.start.transc, max.end.transc, wseq, verbose = TRUE, ...) {
  cnt <- 1
  nameG1    <- vector(mode = "character")
  nameG2    <- vector(mode = "character")
  refOp1    <- vector(mode = "character")
  refOp2	  <- vector(mode = "character")
  diffExpr  <- vector(mode = "numeric")    
  lengthIGR <- vector(mode = "integer")
  cuScore   <- vector(mode = "numeric") 
  exprIGR   <- vector(mode = "numeric")
  class     <- vector(mode = "character")
  for(i in 1:nrow(transcIGRs.pos)) {
    g1 <- transcIGRs.pos$startLocusTag[i]
    g2 <- transcIGRs.pos$endLocusTag[i]
    if(length(which(rownames(genes.and.operons) == g1)) == 0 | length(which(rownames(genes.and.operons) == g2)) == 0) next;
    gene1 <- genes.and.operons[which(rownames(genes.and.operons) == g1),]
    gene2 <- genes.and.operons[which(rownames(genes.and.operons) == g2),]
    testOpDOOR <- (is.na(gene1$operonID) & !is.na(gene2$operonID)) | (!is.na(gene1$operonID) & is.na(gene2$operonID)) | (is.na(gene1$operonID) & is.na(gene2$operonID))
    test_transcActivity_g1 <- transcCDSs[rownames(gene1),]$expr > minExprCDS
    test_transcActivity_g2 <- transcCDSs[rownames(gene2),]$expr > minExprCDS
    test_transcIGR <- transcIGRs.pos$expr[i] > minExprIGR
    test_lenIGR  <- transcIGRs.pos$length[i] < maxLenIGR
    test_confOP  <- !(rownames(gene1) %in% OPs$G1) & !(rownames(gene2) %in% OPs$G2)
    test_confNOP <- !(rownames(gene1) %in% NOPs$G1) & !(rownames(gene2) %in% NOPs$G2)
    test_break_point <- FALSE
    test_feat <- gene1$feature == "CDS" & gene2$feature == "CDS"
    if(length(which(rownames(POSSs) == gene2)) > 0) {
      id_poss <- which(rownames(POSSs) == gene2)
      if(log2(POSSs$lp[id_poss]) <= max.start.transc[1]) 
        test_break_point <- TRUE
    }
    if(length(which(rownames(POESs) == gene1)) > 0) {
      id_poes <- which(rownames(POESs) == gene1)
      if(log2(POESs$lp[id_poes]) <= max.end.transc[1]) 
        test_break_point <- TRUE
    }
    # print(paste(test_confOP, test_confNOP, testOpDOOR, test_transcActivity_g1, test_transcActivity_g2, test_transcIGR))
    if(test_confOP & test_confNOP & testOpDOOR & test_transcActivity_g1 & test_transcActivity_g2 & 
         test_transcIGR & test_lenIGR & !test_break_point & test_feat) {
      geneSeq1 = wseq[gene1$start:gene1$end]
      geneSeq2 = wseq[gene2$start:gene2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt] <- sum(rscu1 * rscu2, na.rm = TRUE)
      lengthIGR[cnt] <- transcIGRs.pos$length[i]
      exprIGR[cnt]   <- transcIGRs.pos$expr[i]
      id_CDS  <- which(rownames(transcCDSs) == rownames(gene1))
      diffExpr[cnt]  <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt]	   <- rownames(gene1)
      nameG2[cnt]	   <- rownames(gene2)
      refOp1[cnt]	   <- gene1$operonID
      refOp2[cnt]    <- gene2$operonID
      class[cnt]	   <-  "Und"
      cnt <- cnt + 1
    }
  }	
  for(i in 1:nrow(transcIGRs.neg)) {
    g1 <- transcIGRs.neg$startLocusTag[i]
    g2 <- transcIGRs.neg$endLocusTag[i]
    if(length(which(rownames(genes.and.operons) == g1)) == 0 | length(which(rownames(genes.and.operons) == g2)) == 0) next;
    gene1 <- genes.and.operons[which(rownames(genes.and.operons) == g1),]
    gene2 <- genes.and.operons[which(rownames(genes.and.operons) == g2),]
    testOpDOOR <- (is.na(gene1$operonID) & !is.na(gene2$operonID)) | (!is.na(gene1$operonID) & is.na(gene2$operonID)) | (is.na(gene1$operonID) & is.na(gene2$operonID))
    test_transcActivity_g1 <- transcCDSs[rownames(gene1),]$expr > minExprCDS
    test_transcActivity_g2 <- transcCDSs[rownames(gene2),]$expr > minExprCDS
    test_transcIGR <- transcIGRs.neg$expr[i] > minExprIGR
    test_lenIGR  <- transcIGRs.neg$length[i] < maxLenIGR
    test_confOP  <- !(rownames(gene1) %in% OPs$G1) & !(rownames(gene2) %in% OPs$G2)
    test_confNOP <- !(rownames(gene1) %in% NOPs$G1) & !(rownames(gene2) %in% NOPs$G2)
    test_break_point <- FALSE
    test_feat <- gene1$feature == "CDS" & gene2$feature == "CDS"
    if(length(which(rownames(POSSs) == gene1)) > 0) {
      id_poss <- which(rownames(POSSs) == gene1)
      if(log2(POSSs$lp[id_poss]) <= max.start.transc[2]) 
        test_break_point <- TRUE
    }
    if(length(which(rownames(POESs) == gene2)) > 0) {
      id_poes <- which(rownames(POESs) == gene2)
      if(log2(POESs$lp[id_poes]) <= max.end.transc[2]) 
        test_break_point <- TRUE
    }
    
    # print(paste(test_confOP, test_confNOP, testOpDOOR, test_transcActivity_g1, test_transcActivity_g2, test_transcIGR))
    if(test_confOP & test_confNOP & testOpDOOR & test_transcActivity_g1 & test_transcActivity_g2 & 
         test_transcIGR & test_lenIGR & !test_break_point & test_feat) {
      ##
      geneSeq1 = wseq[gene1$start:gene1$end]
      geneSeq2 = wseq[gene2$start:gene2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt] <- sum(rscu1 * rscu2, na.rm = TRUE)
      id_IGR   <- which(transcIGRs.neg$startLocusTag == rownames(gene1))
      lengthIGR[cnt] <- transcIGRs.neg$length[id_IGR]
      exprIGR[cnt]   <- transcIGRs.neg$expr[id_IGR]
      id_CDS  <- which(rownames(transcCDSs) == rownames(gene1))
      diffExpr[cnt]  <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt]	   <- rownames(gene1)
      nameG2[cnt]	   <- rownames(gene2)
      refOp1[cnt]	   <- gene1$operonID
      refOp2[cnt]    <- gene2$operonID
      class[cnt]	   <-  "Und"
      cnt <- cnt + 1	
    }
  }	
  if(verbose) cat("    - nameG1:",length(nameG1),", nameG2:",length(nameG2),
                  ", diffExpr:",length(diffExpr),", lengthIGR:",length(lengthIGR),
                  ", cuScore:",length(cuScore),"\n", sep="")
  return(data.frame(G1 = nameG1, 
                    G2 = nameG2, 
                    diffExpr = diffExpr,  
                    exprIGR = exprIGR,  
                    lengthIGR = lengthIGR,
                    cuScore = cuScore,
                    refOp1 = refOp1,
                    refOp2 = refOp2,
                    class = class,
                    stringsAsFactors = FALSE))
}

#' Define the set of OPs annotated in DOOR. 
#' 
#' Build a data table containing the transcriptomic/genomic feature values of operon pairs annotated in DOOR.
#' The operon classifier trained/validated on the OPs and NOPs is used to predict the operon status of these gene pairs.
#' 
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param transcCDSs Transcription levels for the coding regions. See \code{comp.gene.transc.levels}.
#' @param transcIGRs.pos Transcription levels for the intergenic regions (forard strand). See \code{comp.igr.transc.levels}.
#' @param transcIGRs.neg Transcription levels for the intergenic regions (reverse strand). See \code{comp.igr.transc.levels}.
#' @param wseq Sequence genome.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' select.ops.indoor()
select.ops.indoor  <- function(genes.and.operons, transcCDSs, transcIGRs.pos, transcIGRs.neg, wseq, 
                               verbose = TRUE, ...) {
  cnt <- 1
  nameG1    <- vector(mode = "character")
  nameG2    <- vector(mode = "character")
  refOp     <- vector(mode = "character")
  diffExpr  <- vector(mode = "numeric")    
  lengthIGR <- vector(mode = "integer")
  cuScore   <- vector(mode = "numeric") 
  exprIGR   <- vector(mode = "numeric")
  class     <- vector(mode = "character")
  for(i in 1:nrow(transcIGRs.pos)) {
    g1 <- transcIGRs.pos$startLocusTag[i]
    g2 <- transcIGRs.pos$endLocusTag[i]
    gene1 <- genes.and.operons[which(rownames(genes.and.operons) == g1),]
    gene2 <- genes.and.operons[which(rownames(genes.and.operons) == g2),]
    # print(paste(test_confOP, test_confNOP, testOpDOOR, test_transcActivity_g1, test_transcActivity_g2, test_transcIGR))
    if(!is.na(gene1$operonID) & !is.na(gene2$operonID) & (gene1$operonID == gene2$operonID)) {
      geneSeq1 = wseq[gene1$start:gene1$end]
      geneSeq2 = wseq[gene2$start:gene2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt] <- sum(rscu1 * rscu2, na.rm = TRUE)
      lengthIGR[cnt] <- transcIGRs.pos$length[i]
      exprIGR[cnt]   <- transcIGRs.pos$expr[i]
      id_CDS  <- which(rownames(transcCDSs) == rownames(gene1))
      diffExpr[cnt]  <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt]	   <- rownames(gene1)
      nameG2[cnt]	   <- rownames(gene2)
      refOp[cnt]	   <- gene1$operonID
      class[cnt]	   <- "Door"
      cnt <- cnt + 1
    }
  }	
  for(i in 1:nrow(transcIGRs.neg)) {
    g1 <- transcIGRs.neg$startLocusTag[i]
    g2 <- transcIGRs.neg$endLocusTag[i]
    gene1 <- genes.and.operons[which(rownames(genes.and.operons) == g1),]
    gene2 <- genes.and.operons[which(rownames(genes.and.operons) == g2),]
    if(!is.na(gene1$operonID) & !is.na(gene2$operonID) & (gene1$operonID == gene2$operonID)) {
      geneSeq1 = wseq[gene1$start:gene1$end]
      geneSeq2 = wseq[gene2$start:gene2$end]
      rscu1 <- seqinr::uco(geneSeq1, index = "rscu")
      rscu2 <- seqinr::uco(geneSeq2, index = "rscu")
      cuScore[cnt] <- sum(rscu1 * rscu2, na.rm = TRUE)
      id_IGR   <- which(transcIGRs.neg$startLocusTag == rownames(gene1))
      lengthIGR[cnt] <- transcIGRs.neg$length[id_IGR]
      exprIGR[cnt]   <- transcIGRs.neg$expr[id_IGR]
      id_CDS  <- which(rownames(transcCDSs) == rownames(gene1))
      diffExpr[cnt]  <- round(abs(transcCDSs$expr[id_CDS] - transcCDSs$expr[id_CDS+1]),digits = 4)
      nameG1[cnt]	   <- rownames(gene1)
      nameG2[cnt]	   <- rownames(gene2)
      refOp[cnt]	   <- gene1$operonID
      class[cnt]	   <-  "Door"
      cnt <- cnt + 1	
    }
  }	
  if(verbose)	cat("    - nameG1:",length(nameG1),", nameG2:",length(nameG2),
                  ", diffExpr:",length(diffExpr),", lengthIGR:", length(lengthIGR),
                  ", cuScore:",length(cuScore),"\n", sep="")
  
  return(data.frame(G1 = nameG1, 
                    G2 = nameG2, 
                    diffExpr = diffExpr,  
                    exprIGR = exprIGR,  
                    lengthIGR = lengthIGR,
                    cuScore = cuScore,
                    refOp = refOp,
                    class = class,
                    stringsAsFactors = FALSE))
}

#' Normalize the data before the classification step. 
#' 
#' n0 - without normalization
#' n1 - standardization ((x-mean)/sd)
#' n2 - positional standardization ((x-median)/mad)
#' n3 - unitization ((x-mean)/range)
#' n4 - unitization with zero minimum ((x-min)/range)
#' n5 - normalization in range <-1,1> ((x-mean)/max(abs(x-mean)))
#' 
#' @param OPs Data table including the confirmed operon pairs (OPs). See \code{select.ops}.
#' @param OPs Data table including the confirmed non-operon pairs (NOPs). See \code{select.nops}.
#' @param POPs Data table including gene pairs with an operon status to (re)define (POPs). See \code{select.pops}.
#' @param DOPs Data table including the operon pairs annotated in DOOR database (DOPs). See \code{select.ops.indoor}.
#' @param type Charcater vector indicating the method to use for the normalization step. Default value is "n0".
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' pre.processing()
pre.processing <- function(OPs, NOPs, POPs, DOPs, type = "n0", verbose = TRUE, ...) {
  x <- rbind(OPs[,c(3:6)], NOPs[,c(3:6)], POPs[,c(3:6)], DOPs[,c(3:6)])
  if(verbose) cat("    - OPs#", nrow(OPs),", NOPs#", nrow(NOPs), ", POPs#", nrow(POPs),
                  ", DOPs#", nrow(DOPs), "\n", sep="")		
  ## normalization by column
  nx <- NULL
  for (i in 1:ncol(x)) 
    nx <- switch(type, n0 = cbind(nx, (x[,i])), 
                 n1 = cbind(nx, (x[,i] - mean(x[,i]))/(stats::sd(x[,i]))), 
                 n2 = cbind(nx, (x[,i] - stats::median(x[,i]))/(stats::mad(x[,i]))), 
                 n3 = cbind(nx, (x[,i] - mean(x[,i]))/(max(x[,i]) - min(x[,i]))), 
                 n4 = cbind(nx, (x[, i] - min(x[,i]))/(max(x[,i]) - min(x[,i]))), 
                 n5 = cbind(nx, (x[, i] - mean(x[,i]))/(max(abs(x[,i] - mean(x[,i]))))))
  
  nOPs  <- nx[1:nrow(OPs),]
  nNOPs <- nx[(nrow(OPs)+1):(nrow(OPs)+nrow(NOPs)),]
  nPOPs <- nx[(nrow(OPs)+nrow(NOPs)+1):(nrow(OPs)+nrow(NOPs)+nrow(POPs)),]
  nDOPs <- nx[(nrow(OPs)+nrow(NOPs)+nrow(POPs)+1):(nrow(OPs)+nrow(NOPs)+nrow(POPs)+nrow(DOPs)),]
  #if(verbose){
  #	cat("\n OPs#",nrow(nOPs),", NOPs#", nrow(nNOPs),", POPs#",nrow(nPOPs),", DOPs#",nrow(nDOPs),"\n", sep="")
  #}
  for(i in 1:3) OPs[,(i+2)] <- nOPs[,i]
  for(i in 1:3) NOPs[,(i+2)] <- nNOPs[,i]
  for(i in 1:3) POPs[,(i+2)] <- nPOPs[,i]
  for(i in 1:3) DOPs[,(i+2)] <- nDOPs[,i]
  dat <- as.data.frame(rbind(OPs[,c(3:6,ncol(OPs))],NOPs[,c(3:6,ncol(NOPs))]))
  rownames(dat) <- rbind(cbind(paste(OPs[,1], OPs[,2], sep="-")), 
                         cbind(paste(NOPs[,1], NOPs[,2], sep="-")))
  dat$class <- as.factor(dat$class)
  dat.pops  <- as.data.frame(POPs[,c(3:6)])
  rownames(dat.pops) <- paste(POPs[,1], POPs[,2], sep="-")
  dat.dops  <- as.data.frame(DOPs[,c(3:6)])
  rownames(dat.dops) <- paste(DOPs[,1], DOPs[,2], sep="-")
  #if(plot)	plot.feat.vals(dat)
  list(dat = dat, dat.pops = dat.pops, dat.dops = dat.dops)
} 

#' Tune and build the classification models. 
#' 
#' Internal function to tune and train all the classifiers to be used in distinguishing operon pairs (OPs) from non-operon pairs (NOPs)
#' on a given RNA-seq expression profile.
#' 
#' @param data Training/test data. See \code{select.ops} and \code{select.nops}.
#' @param class Vector of the class labels.
#' @param nf Number of folds for the cross-validation and the automatic selection of the model parameters.
#' @keywords internal
#' @author Vittorio Fortino
tune.cls <- function(data, class, nf=3, verbose = TRUE, ...) {
  data$class <- class
  ## Random Forest 
  rf.search <- list(smethod="grid",search=list(mtry=c(2:4)), convex=0,  metric="ACC", method = c("kfold", nf, 12345))
  tr.models <- list()
  tr.models[[1]] <- rminer::fit(class~., data, model="mlpe", task = "class", search="heuristic5", mpar=c(NA,NA,"kfold",nf,"ACC")) ## MLPE
  tr.models[[2]] <- rminer::fit(class~., data, model="randomforest", task = "class", search = rf.search) ## Random Forest
  tr.models[[3]] <- rminer::fit(class~., data, model="svm", task = "class", search="heuristic5", mpar=c(NA,NA,"kfold",nf,"ACC")) ## SVM
  names(tr.models) <- c("mlpe","rf","svm")
  if(verbose) {
    cat("    - Automatic parameter setting: \n", sep="")
    cat("       * MLPE-model: ", paste(paste(names(tr.models[[1]]@mpar), tr.models[[1]]@mpar, sep="="), collapse=" "), "\n", sep="")
    cat("       * RF-model: ", paste(paste(names(tr.models[[2]]@mpar), tr.models[[2]]@mpar, sep="="), collapse=" "), "\n", sep="")
    cat("       * SVM-model: ", paste(paste(names(tr.models[[3]]@mpar), tr.models[[3]]@mpar, sep="="), collapse=" "), "\n", sep="")
  }
  #list(mlpe = tr.models[[1]]@mpar, rf = tr.models[[2]]@mpar, svm = tr.models[[3]]@mpar)
  tr.models
}

#' Validate the classification models. 
#' 
#' Internal function to validate all the classifiers to be used in distinguishing operon pairs (OPs) from non-operon pairs (NOPs)
#' on a given RNA-seq expression profile.
#' 
#' @param data Training/test data. See \code{select.ops} and \code{select.nops}.
#' @param class Vector of the class labels.
#' @param runs Number of bootstraps to be used.
#' @param kf Number of folds for the cross-validation.
#' @keywords internal
#' @author Vittorio Fortino
validate.cls <- function(data, class, models, runs = 5, kf = 5, verbose = TRUE) {
  data$class <- class
  vl.models <- list()
  vl.models[[1]] <- rminer::mining(class~., data, task = "class", Runs=runs, method=c("kfold",kf), model="mlpe", 
                                   mpar = c(models[[1]]@mpar$nr, models[[1]]@mpar$maxit)) 
  vl.models[[2]] <- rminer::mining(class~., data, task = "class", Runs=runs, method=c("kfold",kf), model="randomforest", 
                                   mpar = c(models[[2]]@mpar$mtry))
  vl.models[[3]] <- rminer::mining(class~., data, task = "class", Runs=runs, method=c("kfold",kf), model="svm", 
                                   mpar = c(models[[3]]@mpar$C, models[[3]]@mpar$epsilon))
  acc <- rep(0, 3)
  sd.acc <- rep(0, 3)
  gm <- rep(0, 3)
  sd.gm <- rep(0, 3)
  err <- rep(0, 3)
  sd.err <- rep(0, 3)
  fsc <- rep(0, 3)
  sd.fsc <- rep(0, 3)
  names(acc) <- c("mlpe","rf","svm")
  names(sd.acc) <- c("mlpe","rf","svm") 
  names(gm) <- c("mlpe","rf","svm")
  names(sd.gm) <- c("mlpe","rf","svm")
  names(err) <- c("mlpe","rf","svm")
  names(sd.err) <- c("mlpe","rf","svm")
  names(fsc) <- c("mlpe","rf","svm")
  names(sd.fsc) <- c("mlpe","rf","svm")
  for(i in 1:length(vl.models)) {
    conf.list <- rminer::mmetric(vl.models[[i]], metric="CONF")
    metrics.list <- lapply(conf.list, function(cm) {
      metrics <- sum.conf.matrix(cm$conf)
      metrics
    })
    acc[i]    <- mean(unlist(lapply(metrics.list, function(m) m$acc)))
    sd.acc[i] <- stats::sd(unlist(lapply(metrics.list, function(m) m$acc)))
    gm[i]     <- mean(unlist(lapply(metrics.list, function(m) m$gmean)))
    sd.gm[i]  <- stats::sd(unlist(lapply(metrics.list, function(m) m$gmean)))
    err[i]    <- mean(unlist(lapply(metrics.list, function(m) m$err)))
    sd.err[i] <- stats::sd(unlist(lapply(metrics.list, function(m) m$err)))
    fsc[i]    <- mean(unlist(lapply(metrics.list, function(m) mean(m[[1]]$Fscore))))
    sd.fsc[i] <- stats::sd(unlist(lapply(metrics.list, function(m) mean(m[[1]]$Fscore))))
  }
  return(data.frame(acc = acc, sd.acc = sd.acc, 
                    gmean = gm, sd.gmean = sd.gm,
                    fscore = fsc, sd.fscore = sd.fsc,
                    err = err, sd.err = sd.err))
}

#' Provide classification metrics. 
#' 
#' Internal function to compile common classification metrics on a give confusion matrix. 
#' 
#' @param cm confusion matrix. 
#' @keywords internal
#' @author Vittorio Fortino
sum.conf.matrix <- function(cm, ...) {
  ## Number of groups
  Ngp <- nrow(cm)
  ## Total : TP + TN + FP + FN
  Tot <- sum(cm)
  ## TP : True positive item : All items on diagonal
  TP <- diag(cm)
  ## TP + TN : sum of diagonal = All correct identification
  TP_TN <- sum(TP)
  ## TP + FP : sum of columns : Automatic classification
  TP_FP <- colSums(cm)
  ## TP + FN : sum of rows : Manual classification
  TP_FN <- rowSums(cm)
  ## FP : False positive items
  FP <- TP_FP - TP   
  ## FN : False negative item
  FN <- TP_FN - TP
  ## TN : True Negative = Total - TP - FP - FN
  TN <- rep(Tot, Ngp) - TP - FP - FN
  ## The 8 basic ratios
  ## Recall = TP / (TP + FN) = 1 - FNR
  Recall <- TP / (TP_FN)
  ## Specificity = TN / (TN + FP) = 1 - FPR
  Specificity <- TN / (TN + FP)
  ## Precision = TP / (TP + FP) = 1 - FDR
  Precision <- TP / (TP_FP)
  ## F-score = F-measure = F1 score = Harmonic mean of Precision and recall
  Fscore <- 2 * ((Precision * Recall) / (Precision + Recall))
  Fscore[is.nan(Fscore)] <- 0
  Gmean <- prod(Recall)^(1/Ngp)
  ## General statistics
  Accuracy <- ((TP + TN) / (TP + TN + FP + FN))
  Error <- 1 - Accuracy
  ## Create a data frame with all results
  res <- list(data.frame(Fscore = Fscore, Recall = Recall, Precision = Precision, Specificity = Specificity), 
              acc = Accuracy, err = Error, gmean = Gmean)
  return(res)
}

#' Train and validate the operon classifier and evaluate the feature importance. 
#' 
#' Internal function to train an operon classifier to distinguish operon pairs (OPs) from non-operon pairs (NOPs)
#' on a given RNA-seq expression profile.
#' 
#' @param data Training/test data. See \code{select.ops} and \code{select.nops}.
#' @param class Vector of the class labels.
#' @param run   Number of runs of training/validation. Default values is 30.
#' @param p     Percentage of samples to be used for the training. Default values is 0.9.
#' @param ntr   Number of trees to use in the operon classifier. Default values is 1000.
#' @param mtree Number of variables randomly sampled as candidates at each split in each tree. Default values is 4. 
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
train.RFs <- function(data, class, run = 30, p = 0.9, ntr = 1000, mtree=4, verbose = TRUE, ...) {
  rank.aggreg <- rep(0, ncol(data))
  names(rank.aggreg) <- colnames(data)
  bootSam <- round(table(class)*p) 
  print(bootSam)
  if(verbose) cat("    - Training data: ", bootSam[1],  "(NOP) - ",  bootSam[2], "(OP)\n", sep="")
  nclass  <- names(table(class))
  id.class <- list()
  for(i in 1:length(nclass)) {
    id.class[[i]] <- which(nclass[i] == class)
  }
  accs <- rep(0, run)
  for(i in 1:run) {
    list.ids.train <- list()
    list.ids.test  <- list()
    for(j in 1:length(nclass)) {
      list.ids.train[[j]] <- sample(id.class[[j]],bootSam[j])
      list.ids.test[[j]]  <- id.class[[j]][which((id.class[[j]] %in% list.ids.train[[j]]) == FALSE)]
    }
    rf <- randomForest::randomForest(x=data[unlist(list.ids.train),], y=class[unlist(list.ids.train)], 
                       xtest=data[unlist(list.ids.test),], ytest=class[unlist(list.ids.test)], 
                       ntree=ntr, mtry=mtree, importance=TRUE)
    accs[i] <- mean(rf$test$predicted == class[unlist(list.ids.test)])
    rf.imp <- randomForest::importance(rf, type=1)
    sort.variables <- rownames(rf.imp)[sort(rf.imp, decreasing = TRUE, index.return = TRUE)[[2]]]
    for(j in 1:dim(rf.imp)[1]) {
      rank.aggreg[sort.variables[j]] <- rank.aggreg[sort.variables[j]] + (dim(rf.imp)[1] - j)
    }
  }
  rank.aggreg <- rank.aggreg[sort(rank.aggreg, decreasing = TRUE, index.return = TRUE)[[2]]]
  if(verbose) {
    cat("      * mean ACC: ", mean(accs), ", sd ACC: ", stats::sd(accs),"\n", sep="")
    cat("      * feature importance: \n", sep="")
    for(i in 1:length(rank.aggreg))	{
      cat("        #",names(rank.aggreg)[i],":",rank.aggreg[i],"\n", sep="")
    }
  }
  rf <- randomForest(x=data, y=class, ntree=ntr, mtry=mtree)              
  return(list(model=rf, accs=accs, rf=rank.aggreg))
}

#' Predict operon status of gene pairs (e.g., POPs, DOPs,...).
#' 
#' Internal function that uses the trained operon classifier to predict the operon status of POPs and DOPs.
#' 
#' @param cls List of classifiers. See \code{train.rfs}.
#' @param gene.pairs Gene pairs. 
#' @param type Character vector indicating the type of gene pairs (POPs or DOPs). See \code{select.pops} and \code{select.ops.indoor}.
#' @param cons Numeric values indicating the minimum number of postive votes to classify a gen pair as operon pai. See \code{select.pops} and \code{select.ops.indoor}.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
pred.operon.status <- function(list.cls, gene.pairs, type="", cons = 1, verbose=TRUE, ...) {
  pred.ops.mlpe <- predict(list.cls$mlpe, newdata=gene.pairs) 
  pred.ops.rf   <- predict(list.cls$rf,   newdata=gene.pairs) 
  pred.ops.svm  <- predict(list.cls$svm,  newdata=gene.pairs) 
  pred.ops <- data.frame(p.mlpe = pred.ops.mlpe, p.rf = pred.ops.rf, p.svm = pred.ops.svm)
  rownames(pred.ops) <- rownames(gene.pairs)
  consensus <- rep("NOP", length=nrow(pred.ops))
  for(i in 1:nrow(pred.ops)) {
    if(length(which(pred.ops[i,] == "OP")) >= cons) consensus[i] <- "OP"
  }
  pred.ops$cons <- consensus
  if(verbose) 
    cat("    - Number of ", type, " predicted as OPs:", length(which(consensus == "OP")), ".\n", sep="")
  #print(head(pred.ops,50))
  #
  #pred.ops <- predict(cls, newdata=gene.pairs) 
  #print(head(pred.ops))
  #num.pred.ops <- length(which(cbind(pred.ops)[,1] == 2)) 
  #if(verbose) 
  #	cat("    - Number of ", type, " predicted as OPs:", 
  #			num.pred.ops, "(%", (num.pred.ops/length(pred.ops)) * 100, ").\n")
  return(pred.ops)
}

#' Build the condition-dependent operon map for a given RNA-seq expression profile.
#' 
#' @param pred.POPs Data table containing the POPs predicted as OPs.
#' @param pred.DOPs Data table containing the DOPs predicted as OPs.
#' @param genes.and.operons Data table merging gene(s) and operon(s) annotations. See \code{join.genes.and.operons}.
#' @param POSSs Data table representing a set of putative operon start-points.
#' @param POESs Data table representing a set of putative operon end-points.
#' @param max.start.transc Cutoff values for the start transcription points. Default values is 0.1.
#' @param max.end.transc Cutoff values for the end transcription points. Default values is 0.1.
#' @param find.ext Add gene pairs predicted as OPs to confirmed operons. 
#' @param BED.file Save the condition-dependent operon map in a BED-like file.
#' @param TAB.file Save the condition-dependent operon map in a tab-delimited text file file. 
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' getCondOperonMap()
getCondOperonMap <- function(pred.POPs, pred.DOPs, genes.and.operons, POSSs, POESs, 
                             find.ext = FALSE, BED.file = NA, TAB.file = NA,
                             verbose = TRUE, ...) {
  first     <- vector(mode = "character") 
  strand    <- vector(mode = "character") 
  start     <- vector(mode = "numeric")
  end       <- vector(mode = "numeric") 
  length    <- vector(mode = "numeric")
  length.door  <- vector(mode = "numeric")	
  listGenes <- vector(mode = "character") 
  door.op   <- vector(mode = "character")
  type      <- vector(mode = "character")
  op.ref    <- unique(genes.and.operons$operonID)
  ids.remove <- c()
  flag.st.op <- FALSE
  cnt <- 0
  for(i in 1:length(op.ref)) {
    ids <- which(genes.and.operons$operonID == op.ref[i])
    genes <- genes.and.operons[ids,]
    if(nrow(genes) > 0) {
      for(j in 1:(nrow(genes)-1)) {
        if((ids[j] + 1) != ids[j+1]) {
          #print(paste("How many genes between: ", ids[j+1]-ids[j]))
          g1 <- rownames(genes)[ids[j]]
          g2 <- rownames(genes)[(ids[j]+1)] 
          g3 <- rownames(genes)[ids[j+1]]
          t1 <- pred.POPs[paste(g1, g2, sep="-"),]$op.st == "OP" & !flag.st.op
          t2 <- pred.POPs[paste(g2, g3, sep="-"),]$op.st == "OP" & !flag.st.op
          if(t1 & t2)	{
            listGenes[cnt] <- paste(listGenes[cnt], g2, g3, sep="-")
            if(verbose) cat("     A POP has been used to link an operon.\n")
          }
          next
        }
        g1 <- rownames(genes)[j]
        g2 <- rownames(genes)[j+1]
        if(pred.DOPs[paste(g1, g2, sep="-"),]$op.st == "OP" & !flag.st.op) {
          flag.st.op <- TRUE
          cnt <- cnt + 1
          type[cnt] <- "COP"
          first[cnt]    <- g1
          start[cnt]    <- genes[j,]$start
          end[cnt]      <- genes[j+1,]$end
          length[cnt]   <- 2
          length.door[cnt] <- length(ids)
          strand[cnt]   <- genes[j,]$strand
          door.op[cnt]  <- op.ref[i]
          listGenes[cnt] <- paste(g1, g2, sep="-")
        }
        else if(pred.DOPs[paste(g1, g2, sep="-"),]$op.st == "OP" & flag.st.op) {
          listGenes[cnt] <- paste(listGenes[cnt], g2, sep="-")
          end[cnt] <- genes[j+1,]$end
          length[cnt] <- length[cnt] + 1
        }
        else if(pred.DOPs[paste(g1, g2, sep="-"),]$op.st == "NOP" & flag.st.op) { 
          flag.st.op <- FALSE
        }
      }
    }	
    flag.st.op <- FALSE
  } # end for
  cnt.pop <- 0
  if(find.ext) {
    e.cop <- rep(0, nrow(pred.POPs)) 
    cop.e <- rep(0, nrow(pred.POPs)) 
    cons.ops <- rep(0, nrow(pred.POPs)) 
    cop.e[which(pred.POPs$op.st == "OP" & !is.na(pred.POPs$op.ref1) & is.na(pred.POPs$op.ref2))] <- 1
    e.cop[which(pred.POPs$op.st == "OP" & is.na(pred.POPs$op.ref1) & !is.na(pred.POPs$op.ref2))] <- 1
    cons.ops[which(pred.POPs$op.st == "OP" & is.na(pred.POPs$op.ref1) & is.na(pred.POPs$op.ref2))] <- 1
    prev.gene <- ""
    flag.new.op <- FALSE
    list.cons.ops <- list() 
    mask <- rep(0, nrow(pred.POPs))  
    for(i in 1:nrow(pred.POPs)) {
      pop <- pred.POPs[i,]
      g1 <- strsplit(rownames(pop), "-")[[1]][1]
      g2 <- strsplit(rownames(pop), "-")[[1]][2]
      if(prev.gene == g1 & pop$op.st == "OP") {
        if(!flag.new.ops & pred.POPs$op.st[i-1] == "OP") {
          list.cons.ops[[length(list.cons.ops)+1]] <- c(i-1,i)
          flag.new.ops <- TRUE
          mask[c(i-1,i)] <- 1
        }
        else if(flag.new.ops) {
          list.cons.ops[[length(list.cons.ops)]] <- c(list.cons.ops[[length(list.cons.ops)]],i)
          mask[i] <- 1
        }
      }
      else if(prev.gene != g1 | pop$op.st == "NOP"){
        flag.new.ops <- FALSE
      }
      prev.gene <- g2
    }
    ## finding single operon pair
    for(i in 1:nrow(pred.POPs)) {
      pop <- pred.POPs[i,]
      g1 <- strsplit(rownames(pop), "-")[[1]][1]
      g2 <- strsplit(rownames(pop), "-")[[1]][2]
      if(mask[i] == 0){ 
        flag.new.ops <- FALSE
        if(cop.e[i] == 1) {
          id <- which(door.op == pop$op.ref1)
          if(length(id) > 0) {
            id <- id[length(id)] 
            length[id] <- length[id] + 1
            end[id]    <- genes.and.operons[g2,]$end
            listGenes[id] <- paste(listGenes[id], g2, sep="-")
            type[id] <- "COP-E"  
            id.ext.pop <- id
          }
        } 
        else if(e.cop[i] == 1) {
          id <- which(door.op == pop$op.ref2)
          if(length(id) > 0) { 
            id <- id[1]
            length[id] <-  length[id] + 1
            start[id] <- genes.and.operons[g1,]$start
            listGenes[id] <- paste(g1, listGenes[id], sep="-")
            type[id] <- "E-COP"
            first[id] <- g1
            id.ext.pop <- id
          }
        }
        else if(cons.ops[i] == 1) {
          cnt <- cnt + 1
          cnt.pop <- cnt.pop + 1
          type[cnt] <- "POP"
          first[cnt]    <- g1
          start[cnt]    <- genes.and.operons[g1,]$start
          end[cnt]      <- genes.and.operons[g2,]$end
          length[cnt]   <- 2
          length.door[cnt] <- 0
          strand[cnt]   <- genes.and.operons[g1,]$strand
          door.op[cnt]  <- paste("pop", cnt.pop, sep="-")
          listGenes[cnt] <- paste(g1, g2, sep="-")
        }
      }
    }
    for(ops in list.cons.ops) {
      pops <- pred.POPs[ops,]
      # print(pops)
      # exract the genes
      genes <- strsplit(rownames(pops)[1], "-")[[1]]
      for(i in 2:nrow(pops)) {
        genes <- c(genes, strsplit(rownames(pops)[i], "-")[[1]][2])
      }
      # print(genes)
      pp <- as.numeric(as.character(pops$op.ref1[1]))
      ss <- as.numeric(as.character(pops$op.ref2[nrow(pops)]))
      # print(paste("GGG", pp, ss))
      # new putative operon
      if(all(is.na(pops$op.ref1)) & all(is.na(pops$op.ref2))) { 
        cnt <- cnt + 1
        cnt.pop <- cnt.pop + 1
        type[cnt] <- "POP"
        first[cnt]    <- genes[1]
        start[cnt]    <- genes.and.operons[genes[1],]$start
        end[cnt]      <- genes.and.operons[genes[length(genes)],]$end
        length[cnt]   <- nrow(pops)
        length.door[cnt] <- 0
        strand[cnt]   <- genes.and.operons[genes[1],]$strand
        door.op[cnt]  <- paste("pop", cnt.pop, sep="-")
        listGenes[cnt] <- paste(genes, collapse="-")
      }
      else { # there are operon pair annotated in door
        if(!is.na(pops$op.ref1[1]) & all(is.na(pops$op.ref2))) {
          id <- which(door.op == pops$op.ref1[1])
          if(length(id) > 0) {
            id <- id[length(id)] 
            length[id] <- length[id] + (length(genes)-1)
            end[id]    <- genes.and.operons[genes[length(genes)],]$end
            listGenes[id] <- paste(listGenes[id], paste(genes[-1], collapse="-"), sep="-")
            type[id] <- "COP-E"  
          }
        }
        else if(all(is.na(pops$op.ref1)) & !is.na(pops$op.ref2[nrow(pops)])) {
          id <- which(door.op == pops$op.ref2[nrow(pops)])
          if(length(id) > 0) {
            id <- id[1] 
            length[id] <- length[id] + (length(genes)-1)
            start[id]  <- genes.and.operons[genes[1],]$start
            listGenes[id] <- paste(paste(genes[-length(genes)], collapse="-"), listGenes[id], sep="-")
            type[id] <- "E-COP"  
          }
        }
        else if(!is.na(pops$op.ref1[1]) & !is.na(pops$op.ref2[nrow(pops)])) {
          id.1 <- which(door.op == pops$op.ref1[1])
          id.2 <- which(door.op == pops$op.ref2[nrow(pops)])
          # unify a door operon / two consecutive door operons
          # print(id.1)
          # print(id.2)
          # print("-----------------------------------------")
          if(as.character(pops$op.ref1[1]) == as.character(pops$op.ref2[nrow(pops)])) { 
            id.1 <- id.1[length(id.1)]
            end[id.1]    <- genes.and.operons[genes[length(genes)],]$end
            length[id.1] <- length[id.1] + length(genes) - 1
            listGenes[id.1] <- paste(listGenes[id.1], paste(genes[-1], collapse="-"), sep="-")
          }
          else if((pp+1) == ss) {
            if(length(id.1) == 0 & length(id.2) > 0) {
                id.2 <- id.2[1]
                length[id.2] <- length[id.2] + length(genes) - 1 
                start[id.2]    <- genes.and.operons[genes[1],]$start
                listGenes[id.2] <- paste(paste(genes[-length(genes)], collapse="-"), listGenes[id.2], sep="-")
                type[id.2] <- "E-COP" 
                #door.op[id.2] <- paste(pops$op.ref1[1], door.op[id.2], sep="-")
            }
            else if(length(id.1) > 0 & length(id.2) == 0) {
              id.1 <- id.1[length(id.1)]
              length[id.1] <- length[id.1] + (length(genes)-1)
              end[id.1]    <- genes.and.operons[genes[length(genes)],]$end
              listGenes[id.1] <- paste(listGenes[id.1], paste(genes[-1], collapse="-"), sep="-")
              type[id.1] <- "COP-E"  
              #door.op[id.1] <- paste(door.op[id.1], pops$op.ref2[length(genes)], sep="-")
            }
            else if(length(id.1) > 0 & length(id.2) > 0) {
              id.1 <- id.1[length(id.1)]
              id.2 <- id.2[1]
              length[id.1] <- length[id.1] + length[id.2] + (length(genes)-2)
              length.door[id.1] <- length.door[id.1] + length.door[id.2] 
              start[id.1]  <- start[id.1]
              end[id.1]    <- end[id.2]
              listGenes[id.1] <- paste(listGenes[id.1], paste(genes[-c(1,length(genes))], collapse="-"), 
                                       listGenes[id.2], sep="-")
              type[id.1] <- "U-COP"  
              door.op[id.1] <- paste(door.op[id.1], door.op[id.2], sep="-")
              ## remove the id
              ids.remove <- c(ids.remove, id.2)
            }
          }
        }
      }
    }
  }
  if(!is.na(BED.file)) {
    df <- data.frame(chr = rep("chr1",length(start)), start, end, name = paste("Operon",door.op,sep="."), 
                     length, strand, type)
    if(length(ids.remove) > 0) df <- df[-ids.remove,]
    utils::write.table(df, file = paste("COP.", BED.file, sep=""), quote = FALSE, sep = "\t", 
                       row.names = FALSE, col.names = FALSE)
    #utils::write.table(df[which(df$type == "POP"), -ncol(df)], file = paste("POP.", BED.file, sep=""), 
    #            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  df <- data.frame(first = first, strand = strand, start = start, end = end, len = length, len.door = length.door,
                   type = type, list.genes = listGenes, door.op = door.op, stringsAsFactors = FALSE)
  if(length(ids.remove) > 0) df <- df[-ids.remove,]
  if(!is.na(TAB.file)) {
    utils::write.table(df, file = paste("COP.", TAB.file, sep=""), quote = FALSE, sep = "\t", 
                       row.names = FALSE, col.names = TRUE)
  }
  return(df)
} 

#' Get generic info about the door operons that have been confirmed 
#' by using the ensemble operon classifier.
#' 
#' @param genes.and.ops Data table merging gene(s) and operon(s) annotations.
#' @param eval Data table reporting accuracy evaluations.
#' @param pred.POPs Data table containing the POPs predicted as OPs.
#' @param pred.DOPs Data table containing the DOPs predicted as OPs.
#' @param comap Data table representing the condition-dependent operon map.
#' @param verbose Default logical value is TRUE.
#' @keywords internal
#' @author Vittorio Fortino
#' get.info
get.info <- function(genes.and.ops, eval, pred.POPs, pred.DOPs, comap, verbose = TRUE) {
  id.OPs <- 0
  id.NOPs <- 0
  for(i in 1:(nrow(genes.and.ops)-1)) {
    t0 <- !is.na(genes.and.ops$operonID[i]) & !is.na(genes.and.ops$operonID[i+1])
    t1 <- genes.and.ops$operonID[i] == genes.and.ops$operonID[i+1]
    if(t0 & t1)  id.OPs <- id.OPs + 1
    else id.NOPs <- id.NOPs + 1
  }
  preds <- table(pred.DOPs$op.st)
  pc.door.ops <- (preds[2]/id.OPs) * 100
  pc.door.ops.to.nops <- (preds[1]/id.OPs) * 100
  new.preds <- table(pred.POPs$op.st)
  pc.door.nops.to.ops <- (new.preds[2]/id.NOPs) * 100
  if(is.na(pc.door.nops.to.ops)) pc.door.nops.to.ops <- 0
  #pc.door.nops <- (new.preds[1]/id.NOPs) * 100
  ndoor <- length(table(genes.and.ops$operonID))
  conf.complete.door.ops <- (sum(ifelse(comap$len == comap$len.door, 1, 0)) / ndoor) * 100
  conf.split.door.ops <- (length(which(duplicated(comap$door.op) == T & comap$type == "COP")) / ndoor) * 100
  conf.ext.door.ops <- (length(which(comap$type == "E-COP" | comap$type == "COP-E")) / ndoor) * 100
  conf.mer.door.ops <- (length(which(comap$type == "U-COP")) / ndoor) * 100
  new.ops <- (length(which(comap$type == "POP")) / nrow(comap)) * 100
  info <- c(mean(eval$acc), pc.door.ops, pc.door.ops.to.nops, pc.door.nops.to.ops, 
            conf.complete.door.ops,conf.split.door.ops, conf.ext.door.ops, conf.mer.door.ops, new.ops)
  names(info) <-c("over_acc","door_ops", "door_ops_to_nops","door_nops_to_ops",
                  "complete_door_operons", "spl_door_operons","ext_door_operons",
                  "door_operons_merged_into_one","putative_operons")
  info
}

# Gene usage analysis

#' @include Alakazam.R
NULL

#### Classes ####


#### Methods ####


#### Calculation functions ####

# Numerically order V(D)J gene assignments
#
# @param   genes   a character vector containing segment calls
# @return  a character vector of unique gene names ordered by chain, family and gene
orderGenes <- function(genes) {
    # TOPO:  replace this with tigger::sortAlleles()
    # Define sorting regular expression
    lex_regex <- "(IG[HLK]|TR[ABGD])[VDJ]"
    num1_regex <- "(?<=[HLKABGD][VDJ])[0-9]+"
    num2_regex <- "[0-9]+"
    sub_regex <- "(IG[HLK]|TR[ABGD])[VDJ][0-9]+[-/a-zA-Z]*"
    
    # Extract lexical and numerical sorting values from gene names
    gene_set <- unique(genes)
    lex <- stri_extract_first_regex(gene_set, lex_regex)
    sub <- stri_replace_first_regex(gene_set, sub_regex, "")
    num1 <- as.numeric(stri_extract_first_regex(gene_set, num1_regex))
    num2 <- as.numeric(stri_extract_first_regex(sub, num2_regex))
    
    # Sort on chain, family and gene
    gene_set <- gene_set[order(lex, num1, num2)]
    
    return(gene_set)
}


#' Tabulates V(D)J allele, gene or family usage.
#' 
#' Determines the count and relative abundance of V(D)J alleles, genes or families within
#' groups.
#'
#' @param    data    data.frame with Change-O style columns containing clonal assignments.
#' @param    gene    column containing allele assignments. Only the first allele in the
#'                   column will be considered.
#' @param    groups  columns containing grouping variables. If \code{NULL} do not group.
#' @param    copy    name of the \code{data} column containing copy numbers for each 
#'                   sequence. If this value is specified, then total copy abundance
#'                   is determined by the sum of copy numbers within each gene.
#' @param    mode    one of \code{c("gene", "family", "allele")} defining
#'                   the degree of specificity regarding allele calls. Determines whether 
#'                   to return counts for genes, families or alleles.
#' 
#' @return   A data.frame summarizing family, gene or allele counts and frequencies with
#'           columns:
#'           \itemize{
#'             \item \code{GENE}:        name of the family, gene or allele
#'             \item \code{SEQ_COUNT}:   total number of sequences for the gene.
#'             \item \code{SEQ_FREQ}:    frequency of the gene as a fraction of the total
#'                                       number of sequences within each grouping.
#'             \item \code{COPY_COUNT}:  sum of the copy counts in the \code{copy} column.
#'                                       for each gene. Only present if the \code{copy} 
#'                                       argument is specified.
#'             \item \code{COPY_FREQ}:   frequency of the gene as a fraction of the total
#'                                       copy number within each group. Only present if 
#'                                       the \code{copy} argument is specified.
#'           }
#'           Additional columns defined by the \code{groups} argument will also be present.
#'
#' @examples
#' # Load example data
#' file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Without copy numbers
#' genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="family")
#' genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="gene")
#' genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="allele")
#'
#' # With copy numbers and multiple groups
#' genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
#'                     copy="DUPCOUNT", mode="family")
#' genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
#'                     copy="DUPCOUNT", mode="gene")
#' genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
#'                     copy="DUPCOUNT", mode="allele")
#' 
#'@export
countGenes <- function(data, gene, groups=NULL, copy=NULL, 
                       mode=c("gene", "allele", "family")) {
    #groups=NULL
    #groups="PRCONS"
    #gene="V_CALL"
    #mode="gene"
    
    # Check input
    mode <- match.arg(mode)
    check <- checkColumns(data, c(gene, groups, copy))
    if (check != TRUE) { stop(check) }
    
    # Extract gene, allele or family assignments
    gene_func <- switch(mode,
                        allele=getAllele,
                        gene=getGene,
                        family=getFamily)
    data[[gene]] <- gene_func(data[[gene]], first=TRUE)
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize(SEQ_COUNT=n()) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(SEQ_COUNT)") %>%
            dplyr::rename_(.dots=c("GENE"=gene))
    } else {
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize_(SEQ_COUNT=interp(~length(x), x=as.name(gene)),
                              COPY_COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT")),
                           COPY_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COPY_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(COPY_COUNT)") %>%
            dplyr::rename_(.dots=c("GENE"=gene))
    }
    
    # Order genes
    gene_tab$GENE <- factor(gene_tab$GENE, levels=orderGenes(gene_tab$GENE))
    
    return(gene_tab)
}
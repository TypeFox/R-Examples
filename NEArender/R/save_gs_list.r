#' Create a TAB-delimited text file from AGS or FGS
#'
#' Each line in this file represents one gene/protein from an AGS/FGS and is accompanied with respective AGS/FGS ID. This format can be used e.g. as input at web site EviNet \url{https://www.evinet.org/}

#' @param gs.list a list created with \code{\link{samples2ags}}, \code{\link{mutations2ags}}, \code{\link{as_genes_fgs}}, or \code{\link{import.gs}}.
#' @param File output file name.

#' @seealso \code{\link{samples2ags}}, \code{\link{mutations2ags}}, \code{\link{as_genes_fgs}}, \code{\link{import.gs}}
#' @references \url{http://www.biomedcentral.com/1471-2105/13/226}
#' @references \url{https://www.evinet.org/}

#' @examples
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath);
#' fgs.genes <- as_genes_fgs(net);
#' save_gs_list(fgs.genes, File = "single_gene_ags.groups.tsv");
#' @importFrom utils write.table
#' @export


save_gs_list <- function(gs.list, File = "gs.list.groups") {
  t1 <- NULL;
  for (gs in names(gs.list)) {
    t1 <- rbind(t1, cbind(gs.list[[gs]], gs));
  }
  write.table(t1, file=File, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE)
}

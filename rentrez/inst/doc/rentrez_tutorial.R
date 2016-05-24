## ---- count_recs, echo=FALSE---------------------------------------------
library(rentrez)
count_recs <- function(db, denom) {
    nrecs <-  rentrez::entrez_db_summary(db)["Count"]
    round(as.integer(nrecs)/denom, 1)
}

## ---- dbs----------------------------------------------------------------
entrez_dbs()

## ---- cdd----------------------------------------------------------------
entrez_db_summary("cdd")

## ---- sra_eg-------------------------------------------------------------
entrez_db_searchable("sra")

## ----eg_search-----------------------------------------------------------
r_search <- entrez_search(db="pubmed", term="R Language")

## ----print_search--------------------------------------------------------
r_search

## ----search_ids----------------------------------------------------------
r_search$ids

## ----searchids_2---------------------------------------------------------
another_r_search <- entrez_search(db="pubmed", term="R Language", retmax=40)
another_r_search

## ---- Tt-----------------------------------------------------------------
entrez_search(db="sra",
              term="Tetrahymena thermophila[ORGN]",
              retmax=0)

## ---- Tt2----------------------------------------------------------------
entrez_search(db="sra",
              term="Tetrahymena thermophila[ORGN] AND 2013:2015[PDAT]",
              retmax=0)

## ---- Tt3----------------------------------------------------------------
entrez_search(db="sra",
              term="(Tetrahymena thermophila[ORGN] OR Tetrahymena borealis[ORGN]) AND 2013:2015[PDAT]",
              retmax=0)

## ---- sra_searchable-----------------------------------------------------
entrez_db_searchable("sra")

## ---- mesh---------------------------------------------------------------
entrez_search(db   = "pubmed",
              term = "(vivax malaria[MeSH]) AND (folic acid antagonists[MeSH])")

## ---- connectome, fig.width=5, fig.height=4, fig.align='center'----------
search_year <- function(year, term){
    query <- paste(term, "AND (", year, "[PDAT])")
    entrez_search(db="pubmed", term=query, retmax=0)$count
}

year <- 2008:2014
papers <- sapply(year, search_year, term="Connectome", USE.NAMES=FALSE)

plot(year, papers, type='b', main="The Rise of the Connectome")

## ----elink0--------------------------------------------------------------
all_the_links <- entrez_link(dbfrom='gene', id=351, db='all')
all_the_links

## ----elink_link----------------------------------------------------------
all_the_links$links

## ---- elink_pmc----------------------------------------------------------
all_the_links$links$gene_pmc[1:10]

## ---- elink_omim---------------------------------------------------------
all_the_links$links$gene_clinvar


## ---- elink1-------------------------------------------------------------
nuc_links <- entrez_link(dbfrom='gene', id=351, db='nuccore')
nuc_links
nuc_links$links

## ---- elinik_refseqs-----------------------------------------------------
nuc_links$links$gene_nuccore_refseqrna

## ---- outlinks-----------------------------------------------------------
paper_links <- entrez_link(dbfrom="pubmed", id=25500142, cmd="llinks")
paper_links

## ---- urls---------------------------------------------------------------
paper_links$linkouts

## ----just_urls-----------------------------------------------------------
linkout_urls(paper_links)

## ---- multi_default------------------------------------------------------
all_links_together  <- entrez_link(db="protein", dbfrom="gene", id=c("93100", "223646"))
all_links_together
all_links_together$links$gene_protein

## ---- multi_byid---------------------------------------------------------
all_links_sep  <- entrez_link(db="protein", dbfrom="gene", id=c("93100", "223646"), by_id=TRUE)
all_links_sep
lapply(all_links_sep, function(x) x$links$gene_protein)

## ---- Summ_1-------------------------------------------------------------
taxize_summ <- entrez_summary(db="pubmed", id=24555091)
taxize_summ

## ---- Summ_2-------------------------------------------------------------
taxize_summ$articleids

## ---- Summ_3-------------------------------------------------------------
taxize_summ$pmcrefcount

## ---- multi_summ---------------------------------------------------------
vivax_search <- entrez_search(db = "pubmed",
                              term = "(vivax malaria[MeSH]) AND (folic acid antagonists[MeSH])")
multi_summs <- entrez_summary(db="pubmed", id=vivax_search$ids)

## ---- multi_summ2--------------------------------------------------------
extract_from_esummary(multi_summs, "fulljournalname")

## ---- multi_summ3--------------------------------------------------------
date_and_cite <- extract_from_esummary(multi_summs, c("pubdate", "pmcrefcount",  "title"))
knitr::kable(head(t(date_and_cite)), row.names=FALSE)

## ---- transcript_ids-----------------------------------------------------
gene_ids <- c(351, 11647)
linked_seq_ids <- entrez_link(dbfrom="gene", id=gene_ids, db="nuccore")
linked_transripts <- linked_seq_ids$links$gene_nuccore_refseqrna
head(linked_transripts)

## ----fetch_fasta---------------------------------------------------------
all_recs <- entrez_fetch(db="nuccore", id=linked_transripts, rettype="fasta")
class(all_recs)
nchar(all_recs)

## ---- peak---------------------------------------------------------------
cat(strwrap(substr(all_recs, 1, 500)), sep="\n")

## ---- Tt_tax-------------------------------------------------------------
Tt <- entrez_search(db="taxonomy", term="(Tetrahymena thermophila[ORGN]) AND Species[RANK]")
tax_rec <- entrez_fetch(db="taxonomy", id=Tt$ids, rettype="xml", parsed=TRUE)
class(tax_rec)

## ---- Tt_list------------------------------------------------------------
tax_list <- XML::xmlToList(tax_rec)
tax_list$Taxon$GeneticCode

## ---- Tt_path------------------------------------------------------------
tt_lineage <- tax_rec["//LineageEx/Taxon/ScientificName"]
tt_lineage[1:4]

## ---- Tt_apply-----------------------------------------------------------
XML::xpathSApply(tax_rec, "//LineageEx/Taxon/ScientificName", XML::xmlValue)

## ---- asthma-------------------------------------------------------------
upload <- entrez_post(db="omim", id=600807)
upload

## ---- snail_search-------------------------------------------------------
entrez_search(db="nuccore", term="COI[Gene] AND Gastropoda[ORGN]")

## ---- snail_history------------------------------------------------------
snail_coi <- entrez_search(db="nuccore", term="COI[Gene] AND Gastropoda[ORGN]", use_history=TRUE)
snail_coi
snail_coi$web_history

## ---- asthma_links-------------------------------------------------------
asthma_clinvar <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", id=600807)
asthma_clinvar$web_histories

## ---- asthma_links_upload------------------------------------------------
asthma_variants <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", web_history=upload)
asthma_variants

## ---- links--------------------------------------------------------------
snp_links <- entrez_link(dbfrom="clinvar", db="snp", 
                         web_history=asthma_variants$web_histories$omim_clinvar,
                         cmd="neighbor_history")
snp_summ <- entrez_summary(db="snp", web_history=snp_links$web_histories$clinvar_snp)
knitr::kable(extract_from_esummary(snp_summ, c("chr", "fxn_class", "global_maf")))


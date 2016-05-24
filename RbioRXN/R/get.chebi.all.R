get.chebi.all <-
function() {
  url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl'
  url_chebi_compound = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz'
  url_chebi_chemical_data = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv'
  tmpdest = tempfile(pattern = "chebi")
  tmpdest_chebi_compound = tempfile(pattern='chebi_compound')
  tmpdest_chebi_chemical_data = tempfile(pattern='chemical_data')
  download.file(url, destfile=tmpdest)
  download.file(url_chebi_compound, destfile=tmpdest_chebi_compound)
  download.file(url_chebi_chemical_data, destfile = tmpdest_chebi_chemical_data)
  owl = readLines(tmpdest)
  chebi_compound = read.delim(tmpdest_chebi_compound, stringsAsFactors=F, quote="", comment.char="")
  chebi_chemical_data = read.delim(tmpdest_chebi_chemical_data, stringsAsFactors=F, quote="", comment.char="")
  chebi = .parse.chebi(owl, chebi_compound, chebi_chemical_data)
  return(chebi)
}

# chebi.candidates
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

## Returns the possible ChEBI names based on metabolite synonyms
chebi.candidates <- function(metabolite) {
  # Load ChEBI data
  data("chebi", envir = environment())
  chebi <- chebi
  # Search in metabolite synonyms the metabolite name
  chebi$name[grep(metabolite, chebi$synonyms, ignore.case = TRUE)]
}

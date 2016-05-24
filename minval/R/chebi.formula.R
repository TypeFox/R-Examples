# chebi.formula
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

## Returns the molecular formula associated to a ChEBI metabolite name
chebi.formula <- function(metabolite){
  # Load ChEBI data
  data("chebi", envir = environment())
  chebi <- chebi
  # Search in ChEBI database for a molecular formula based in metabolite name
  .safe.index(chebi[chebi$name%in%tolower(metabolite),4],1)
}

man.annot <-
function()
{
  write("\t\t#################################", file="")
  write("\t\t#       Annot Quick Start       #", file="")
  write("\t\t#################################", file="")
  write("\n\t# The lists of IDs must be in text tabulated format and placed alone in a folder.", file="")
  write("\n\t# The annotation file must be in text tabulated format and placed outside of the precedent folder.", file="")
  write("\n\t# The common identifiers must be placed in the first column in the lists and in the annotations file.", file="")
  write("\n\t# - Prototype: annot(annot_file=\"\", data_dir=\"\")", file="")
  write("\t# - annot_file: path for the annotations file.", file="")
  write("\t# - data_dir: path of the folder where are placed the lists.", file="")
  write("\t# - res_path: path of the resulting annotated lists.", file="")
  write("\t# - display [TRUE/FALSE]: if informations could be displayed during the process.", file="")
}


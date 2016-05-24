# @ http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#load
#solar> help file-pedigree                                                       
#
# The pedigree file consists of one record for each individual in the data
# set.  Each record must include the following fields:
#
#    ego ID, father ID, mother ID, sex
#
# In addition, a family ID is required when ego IDs are not unique across
# the entire data set.  If the data set contains genetically identical
# individuals, an MZ-twin ID must be present (as described below).  If an
# analysis of household effects is planned, a household ID can be included
# (also described below).
#
# The default field names are ID, FA, MO, SEX, FAMID, MZTWIN, and HHID.
#

fields <- c("id", "ID", "ids",
  "famid", "FAMID", "famidity",
  "mo", "MO", "mother", "MOTHER", "MOtrait", "motherland",
  "fa", "FA", "father", "FATHER", "fatherland",
  "sex", "SEX", "sexo")

### ID
# pass: id, ID
# filter: ids
grep("^id$|^ID$", fields, value = TRUE)

### FAMID
grep("^famid$|^FAMID$", fields, value = TRUE)

### MO
grep("^mo$|^MO$|^mother$|^MOTHER$", fields, value = TRUE)

### FA
grep("^fa$|^FA$|^father$|^FATHER$", fields, value = TRUE)

### SEX
grep("^sex$|^SEX$", fields, value = TRUE)


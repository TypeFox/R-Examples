# function to convert cancer type to ICD code

# cancers from IARC
#        All Sites, Oral Cavity and Pharynx, Oesophagus. Stomach, Colon, Rectum and Anus, 
#        Liver, Gallbladder, Pancreas, Larynx, Lung, Bone, Melanoma of skin, 
#        Prostate (Males only), Testis (Males only), Breast (Females only), Cervix uteri (Females only), 
#        Corpus uteri (Females only), Ovary and other uterine adnexa (Females only), Kidney, 
#        Bladder, Eye, Brain and Central Nervous System, Thyroid, Non-Hodgkin Lymphoma, 
#        Hodgkin Lymphoma, Multiple myeloma, Leukaemia.

# sites from SEER data (RESPIR*)
#  [1] "9999" "C220" "C300" "C301" "C310" "C311" "C312" "C313" "C318" "C319" "C320" "C321"
# [13] "C322" "C323" "C328" "C329" "C33 " "C340" "C341" "C342" "C343" "C348" "C349" "C37 "
# [25] "C381" "C382" "C383" "C384" "C388" "C390" "C398" "C399" "C419" "C449" "C509" "C709"
# [37] "C73 " "D020" "D021" "D022" "D023" "D038"

#returns a regular expression pattern
# D-sites are ignored, only malignant neoplasm considered
siteLookup <- function(cancer) {
  cancer <- gsub(" ","",tolower(cancer), fixed=TRUE)
  if (cancer=="all sites")      return("")
  if (cancer=="liver")          return("C22[[:digit:]]")
  if (cancer=="nasal and ear")  return("C30[[:digit:]]")
  if (cancer=="sinus")          return("C31[[:digit:]]")
  if (cancer=="larynx")         return("C32[[:digit:]]")
  if (cancer=="trachea")        return("C33")
  if (cancer=="lung")           return("C34[[:digit:]]")
  if (cancer=="thymus")         return("C37")
  if (cancer=="heart")          return("C38[[:digit:]]")
  if (cancer=="other respiratory") return("C39[[:digit:]]")
  if (cancer=="bone")           return("C41[[:digit:]]")
  if (cancer=="breast")         return("C509")
  if (cancer=="meninges")       return("C709")
  if (cancer=="thyroid")        return("C73")
  
  return("Cancer not found, listCancers() for list of supported cancers.")
}

# List cancers types that are supported by this package
# Can specify your own if you know the specific ICD code(s) (as a regular expression)
listCancers <- function() {
  clist = c("All Sites","Liver","Nasal and Ear","Sinus","Larynx","Trachea","Lung","Thymus","Heart","Other Respiratory","Bone","Breast","Meninges","Thyroid")
  return(sort(clist))
}



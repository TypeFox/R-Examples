format_detector<-function(file){
  

  if(grepl(".txt",file)){ return("TXT")}
  else if((grepl(".xls",file))  & !(grepl(".xlsx",file)) ){return("XLS")}
  else if(grepl(".xlsx",file)){return("XLSX")}
  else if(grepl(".sav",file)){return("SPSS_.sav")}
  else if(grepl(".mdb",file)){return("MICROSOFT_OFFICE_ACCES_.mdb")}
  else if(grepl(".dta",file)){return("STATA")}
  else if(grepl(".csv",file)){return("CSV")}
  else if(grepl(".Rda",file)){return("R_.Rda")}
  else if(grepl(".sas7bdat",file)){return("SAS_.sas7bdat")}
  else if(grepl(".xport",file)){return("SAS_.xport")}
  else { return("Specify format")}
  
}

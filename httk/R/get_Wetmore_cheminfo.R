# This function displays the information specified in "info=" for all chemicals from Wetmore et al. (2012) and (2013)
get_wetmore_cheminfo <- function(info="CAS",species="Human")
{
  Wetmore.data <- Wetmore.data
  valid.info <- c("Compound","CAS","MW","Raw.Experimental.Percentage.Unbound","Entered.Experimental.Percentage.Unbound","Fub","source_PPB","Renal_Clearance","Met_Stab","Met_Stab_entered" ,"r2","p.val","Concentration..uM.","Css_lower_5th_perc.mg.L.","Css_median_perc.mg.L.","Css_upper_95th_perc.mg.L.","Css_lower_5th_perc.uM.","Css_median_perc.uM.","Css_upper_95th_perc.uM.","Species")

  if (any(!(info %in% valid.info))) stop(paste("Data on",info[!(info %in% valid.info)],"not available. Valid options are:",paste(valid.info,collapse=" ")))
  
  if (!(toupper(species) %in% toupper(unique(Wetmore.data[,"Species"])))) stop(paste("Species",species,"not found.  Available data for:",unique(Wetmore.data[,"Species"])))
  
  return(unique(Wetmore.data[toupper(Wetmore.data[,"Species"])==toupper(species),info]))
}
# Method from Kilford et al. (2008) for fraction of unbound chemical in the 
#  hepatocyte intrinsic clearance assay
# value for Vr (ratio of ceel volume to incubation volume) is taken from 
#  Wetmore (unpublished)
calc_fu_hep <- function(Pow,Vr=0.005,pH=7.4,pKa_Donor=NA,pKa_Accept=NA) 
{
  if (!is_base(pH=pH,pKa_Donor=pKa_Donor,pKa_Accept=pKa_Accept))
  {
    logPD <- log10(calc_dow(Pow,pH=pH,pKa_Donor=pKa_Donor,pKa_Accept=pKa_Accept)) 
  } else logPD <- log10(Pow)
  
  fu_hep <- 1/(1+ 125*Vr*10^(0.072*logPD^2+0.067*logPD-1.126))
  if (fu_hep <0 | fu_hep>1) fu_hep <- 1
  
  return(fu_hep)
}

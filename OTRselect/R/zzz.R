.onAttach <- function( ... )
{
  msg <- paste("\n",
"OTRselect was developed in support of IMPACT, a comprehensive research\n", 
"program that aims to improve the health and longevity of people by\n",
"improving the clinical trial process. To learn more about our \n",
"research and available software visit www.impact.unc.edu. \n\n",
               sep = "")

  packageStartupMessage(msg)
}

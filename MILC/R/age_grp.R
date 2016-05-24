age_grp <-
function(age)
 {
  age_char <- c("<40", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", ">85")
  age_num  <- c(40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 120) 
  if(is.na(age)){agegrp=NA} else  if (age>max(age_num)) {agegrp <- age_char[length(age_char)]} else 
	  {ind_age  <- min(which(age_num>age)) ; agegrp <- age_char[ind_age]}
  return(agegrp)
 }

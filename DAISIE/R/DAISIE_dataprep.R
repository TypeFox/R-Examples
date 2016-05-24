DAISIE_dataprep = function(datatable,island_age,M,number_clade_types = 1,list_type2_clades = NA, prop_type2_pool = "proportional",epss = 1E-5)
{
	number_colonisations = nrow(datatable)
	datalist = list()

	if(number_clade_types == 1)
	{
     datalist[[1]] = list(island_age=island_age,not_present=(M - number_colonisations))
	}
	if(number_clade_types == 2)
	{
 	   number_type2_colonisations = length(list_type2_clades)
   	 number_type1_colonisations = number_colonisations - number_type2_colonisations
		 if (prop_type2_pool == "proportional")
		 {
       	not_present_type1 = roundn((M/number_colonisations) * number_type1_colonisations) - number_type1_colonisations
		    not_present_type2 = roundn((M/number_colonisations) * number_type2_colonisations) - number_type2_colonisations
		 } else {
		    not_present_type1 = roundn(M * (1 - prop_type2_pool)) - number_type1_colonisations
    		not_present_type2 = roundn(M * prop_type2_pool) - number_type2_colonisations
		 }
	   datalist[[1]] = list(island_age = island_age, not_present_type1 = not_present_type1, not_present_type2 = not_present_type2)
	}
	for (i in 1:nrow(datatable))
	{
  	 datalist[[i + 1]] = list(colonist_name = as.character(datatable[i,"Clade_name"]),branching_times = NA,stac = NA,missing_species = datatable[i,"Missing_species"], type1or2 = 1)
	   the_brts = rev(sort(as.numeric(unlist(strsplit(as.character(datatable[i,"Branching_times"]),split = ",")))))
		 if(length(the_brts) == 1) 
		 {
		    datalist[[i + 1]]$branching_times = c(island_age,min(the_brts,island_age - epss))
	   }
  	 if(length(the_brts) > 1)  
		 {
		    datalist[[i + 1]]$branching_times = c(island_age,the_brts)
	   }
		 if(datatable[i,"Status"] == "Non_endemic_MaxAge")
 		 {
		    datalist[[i + 1]]$stac = 1
	   }
	   if(datatable[i,"Status"] == "Endemic")
		 {
		    datalist[[i + 1]]$stac = 2
		 }
		 if(datatable[i,"Status"] == "Endemic&Non_endemic")
		 {
		    datalist[[i + 1]]$stac = 3
		 }
	   if(datatable[i,"Status"] == "Non_endemic")
		 {
		    datalist[[i + 1]]$stac = 4
		 }
		 if(number_clade_types == 2)
		 {
		    if(length(which(list_type2_clades == datatable[i,"Clade_name"])) > 0)
        {
           datalist[[i + 1]]$type1or2 = 2
        }
		 }
	}
  return(datalist)
}
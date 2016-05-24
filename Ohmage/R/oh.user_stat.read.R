oh.user_stat.read <- function(campaign_urn="urn:ohmage:special:all", user="urn:ohmage:special:all", ... ){
	xhr <- oh.call("/user_stats/read", campaign_urn=campaign_urn, user=user, ...);		
	return(xhr);
}

test.mobilize <- function(serverurl, username, password){
	try(oh.logout());
	oh.login(username, password, serverurl);
	campaigns <- oh.campaign.read(output_format="long");
	campaign.names <- names(campaigns$data);
	
	for(i in campaign.names){
		test.campaign(i);		
	}
}

test.campaign <- function(campaign.name){
	cat("testing: ", campaign.name, "\n");
	alldata <- oh.survey_response.read(campaign.name);
	prompt.ids <- grep("prompt.id",names(alldata), value=TRUE);
	prompt.names <- gsub("prompt.id.", "", prompt.ids);
	
	#some campaign plots:
	responseplot(campaign.name);

}
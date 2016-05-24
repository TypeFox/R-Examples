get_PMCIDS = function(abs){
result = NULL;temp2=NULL;for ( i in 1:length(abs@Abstract)){temp = regexpr("PMCID: PMC[1234567890][1234567890][1234567890][1234567890][1234567890][1234567890][1234567890]",abs@Abstract[i]);if (temp != -1){temp1 = list(substr(abs@Abstract[i],temp,temp+17));temp2 = c(temp2,i);result = c(result,temp1);attr(result,"PMID") = abs@PMID[temp2] }}; return(result)}

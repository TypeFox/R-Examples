ecoflora_download_to_local_directory <- function(directory){

    url <- "http://www.ecoflora.co.uk/search_species.php"
    temp_pag <- htmlParse(url)

    query <- paste('//*/a[contains(@href,"search_species2.php?plant_no=")]')
    species <- xpathSApply(temp_pag,query,xmlValue)
    species <- gsub("=.*$","",species,perl=TRUE)
    species <- gsub("\\W+$","",species,perl=TRUE)

    urls<-xpathSApply(temp_pag, query, xmlGetAttr, "href")
    urls<-gsub("search_species2.php?","",urls,fixed=TRUE)

    ECOFLORA_df<-data.frame(species,web_link=urls)
    save(file=file.path(directory,"ECOFLORA_df.Rda"),ECOFLORA_df)
    
}

# a helper method which behaves like dataset_url in Rails
dataset_url <- function(dataset_id, type = c("csv2", "csv", "xml", "xls", "eml", "freeformat"), ...) {
  type = match.arg(type, c("csv2", "csv", "xml", "xls", "eml", "freeformat"))
  seg = switch(type, csv2="/download.csv", csv="/download.csv", xml = ".xml", xls="/download", eml=".eml", freeformat="/freeformats_csv")
  params = Filter(Negate(is.null), list(...))
  if (type %in% c("csv2" ,"eml", "xml")) params$separate_category_columns = TRUE
  query_string = ""
  if (length(params)) query_string = paste("?", paste(names(params), params, sep = "=", collapse = "&"), sep = "")
  url = sprintf("%s/datasets/%d%s%s", bef.options("url"), dataset_id, seg, query_string)
  url = gsub("\\s", "", url)
  return(url)
}

# returns the url of a paper proposal
paperproposal_url <- function(proposal_id, type = c("csv", "xml"),...) {
  type = match.arg(type, c("csv", "xml"))
  seg = switch(type, csv=".csv", xml=".xml")
  params = Filter(Negate(is.null), list(...))
  params$user_credentials = bef.options('user_credentials')
  query_string = ""
  if (length(params)) query_string = paste("?", paste(names(params), params, sep = "=", collapse = "&"), sep = "")
  url = sprintf("%s/paperproposals/%d%s%s", bef.options('url'), proposal_id, seg ,query_string)
  url = gsub("\\s", "", url)
  return(url)
}

# returns the upload url with user credentials
upload_url <- function() {
  base_url = bef.options("url")
  segment = "/datasets/create_with_datafile"
  parameter_sep = "?"
  url = paste0(base_url, segment, parameter_sep, 'user_credentials=', bef.options('user_credentials'))
  return(url)
}

# returns an update url for a dataset
update_url <- function(dataset_id) {
  base_url = bef.options("url")
  segment = paste0("/datasets/", dataset_id)
  parameter_sep = "?"
  url = paste0(base_url, segment, parameter_sep, 'user_credentials=', bef.options('user_credentials'))
  return(url)
}

# handle datagroups related urls
datagroups_url <- function(datagroups_id, type = c("upload", "download")) {
  type = match.arg(type, c("upload", "download"))
  ending = switch(type, download = "/categories.csv", upload = "/update_categories")
  base_url = bef.options("url")
  segment = paste0("/datagroups/", datagroups_id)
  parameter_sep = "?"
  if(type == "upload") {
    url = paste0(base_url, segment, ending, parameter_sep, 'user_credentials=', bef.options('user_credentials'))
  } else {
    url = paste0(base_url, segment, ending, parameter_sep, 'user_credentials=', bef.options('user_credentials'))
  }
  return(url)
}

# a helper that returns the keyword url provided the id
keyword_url <- function(keyword_id) {
  url = sprintf("%s/keywords/%d", bef.options("url"), keyword_id)
  return(url)
}

# function that checks if authentication is given
this_function_requires_api_authentication <- function() {
  user_credentials = bef.options("user_credentials")
  if (is.null(user_credentials) || is.na(user_credentials) || user_credentials == '') {
    stop("This function requires an API key for authentication. Please set your key via bef.options('user_credentials'=yourkey)!")
  }
}

# returns an upload file no matter if it is given as path to a file or data.frame
upload_file <- function(dataset) {
  if(!is.character(dataset) & !is.data.frame(dataset)) {
    stop("Please provide a data frame or a path to a CSV file")
  }
  if (is.data.frame(dataset)) {
    write.csv(dataset, file.path(tempdir(), "upload_dataset.csv"), row.names = FALSE, fileEncoding = "UTF-8", quote = FALSE)
    upload_file = fileUpload(file.path(tempdir(), "upload_dataset.csv"))
  } else {
    upload_file = fileUpload(dataset)
  }
  return(upload_file)
}

# go to the dataset page of the file you just uploaded
bef.goto.dataset_page <- function(id) {
  base_url = bef.options("url")
  segment = "/datasets/"
  id = id
  browseURL(paste0(base_url, segment, id))
}

# go to proposal page
bef.goto.proposal_page <- function(id) {
  base_url = bef.options("url")
  segment = "/paperproposals/"
  id = id
  browseURL(paste0(base_url, segment, id))
}


# Helper that determines internet connection
is_internet_connected <- function() {
  if (.Platform$OS.type == "windows") {
    ipmessage <- system("ipconfig", intern = TRUE)
  } else {
    ipmessage <- system("ifconfig", intern = TRUE)
  }
  validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
  if(!any(grep(validIP, ipmessage[-grep("127.0.0.1", ipmessage)]))) {
    warning("Sorry not internet connection. Please connect first!")
  }
}

# A helper helper that capitalizes the first letter of a string
capitalize <- function(string) {
 first_letter = substr(string, start=0, stop=1)
 rest_of_the_string = substr(string, start=2, stop=nchar(string))
 paste0(toupper(first_letter), rest_of_the_string)
}

# helper to remove nas from lists
remove_na_from_list <- function(object) {
  get_the_nas = lapply(object, function(x) is.na(x))
  remove_the_nas = object[-which(data.frame(do.call(rbind, get_the_nas))$id)]
  return(remove_the_nas)
}

# clean html strings
clean_html_string <- function(string) {
  remove_tags = gsub("<.*?>", "", string)
  remove_newlines = gsub("\\n"," ",remove_tags)
  replace_escaped_quotes = gsub("\"", "'", remove_newlines)
  remove_square_brackets = gsub("\\[.?\\]", "", replace_escaped_quotes)
  remove_empty_brackets = gsub("\\(\\s*\\)", "", remove_square_brackets)
  remove_trailing_ws = gsub("\\s*$", "", remove_empty_brackets)
  remove_leading_ws = gsub("^\\s*", "", remove_trailing_ws)
  remove_unnecessary_ws = gsub("[ ]{2,}", " ", remove_leading_ws)
  cleaned_string = remove_unnecessary_ws
  return(cleaned_string)
}

# get a full list of dataset from portal
all_dataset_of_portal = function() {
  datasets_info_xml = xmlTreeParse(paste0(bef.options("url"), "/datasets.xml"), useInternalNodes = T)
  datasets = xmlToDataFrame(datasets_info_xml, colClasses = c('numeric', 'character'), stringsAsFactors = FALSE)
  return(datasets)
}

# returns all datagroups of the portal
all_datagroups_of_portal = function() {
  this_function_requires_api_authentication()
  datagroup_info_xml = xmlTreeParse(paste0(bef.options("url"), "/datagroups.xml", "?", "user_credentials=", bef.options("user_credentials")), useInternalNodes = T)
  datagroups = xmlToDataFrame(datagroup_info_xml, colClasses = c('numeric', 'character', 'character', 'numeric', 'numeric'), stringsAsFactors = FALSE)
  return(datagroups)
}

# Given a title this returns the dataset id
title_to_dataset_id <- function(dataset_title) {
  all_datasets = all_dataset_of_portal()
  index = match(tolower(dataset_title), tolower(all_datasets$title))
  all_datasets[index, 'id']
}

# Check the portal if a dataset title has been taken or not
the_title_is_taken <- function(dataset_title) {
  all_datasets = all_dataset_of_portal()
  is.element(tolower(dataset_title), tolower(all_datasets$title))
}

# compare two vectors potential candidates for merging
compare_merging_ids <- function(vec1, vec2) {
  matches = (which(!is.na(match(vec1, vec2))))
  percent = (length(matches)/length(vec1)) * 100
  warn = sprintf("Vector one is only an %f percent match to vector 2. Not sufficient for a merge", percent)
  if(length(matches) != length(vec1)) warning(warn)
}

# a helper that returns the metadata attached to a dataset object
bef.extract.metadata <- bef.metadata <- function(dataset) {
  return(attributes(dataset))
}

# a helper function to extract id from a dataset or paperproposal url
url_to_id <- function(url, resource="datasets") {
  pattern = sprintf(".*/%s/(\\d+).*", resource)
  as.numeric(ifelse(grepl(pattern, url), sub(pattern, "\\1", url), NA))
}

# a wrapper function of xpathSApply. which will return NA instead of zero-lenght list when nodes not found
xmlNodesValue <- function(doc, path){
  out = xpathSApply(doc, path, xmlValue, ignoreComments=T)
  out = Filter(function(x) x!="", out)
  if (length(out) == 0) return(NA)
  out
}

# for existing file, append a number to its filename
suggest_filename <- function(filename, dir=getwd()) {
  suggested_filename = filename
  i = 0
  while (file.exists(file.path(dir, suggested_filename))) {
    i = i + 1
    suggested_filename = sub(filename, pattern="(\\.\\w+)?$", replacement=sprintf("(%d)\\1", i))
  }
  return(suggested_filename)
}

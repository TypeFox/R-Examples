#' Upload a dataset to the BEFdata portal
#'
#' This function uploads a dataset from within R to a BEFdata portal set in
#' your options. Before you can do so you need to have an account on the portal
#' and set your credentials via the bef.options() function so the authentication works.
#' The function requires a title for the dataset and either a path to a CSV file for
#' for upload or a dataset object. The dataset must have unique column names.
#'
#' @param dataset The data you like to upload as dataset. This can be a path to a CSV file
#'        or a data frame. Please not that you have to ensure that the column names are unique.
#' @param dataset_title The title of the dataset.
#' @param curl You can pass in a curl handle with additional options. By default a curl handle is
#'        used to improve the memory footprint.
#' @param open_browser If this is set to true the page of the dataset is opened in the
#'        browser after successful upload. This defaults to FALSE.
#' @return Returns a status message with the ID of the dataset.
#' @examples \dontrun{
#'  dataset = data.frame(a = c(1,2,3,4), b = c("a","b","c","d"))
#'  bef.portal.upload.dataset(dataset = dataset, dataset_title = "Testtitle", open_browser = T)
#'       }
#' @import RCurl
#' @import XML
#' @export bef.portal.upload.dataset bef.upload.dataset
#' @aliases bef.upload.dataset

bef.portal.upload.dataset <- bef.upload.dataset <- function(dataset, dataset_title, curl = getCurlHandle(), open_browser = F) {
  this_function_requires_api_authentication()
  if(the_title_is_taken(dataset_title = dataset_title)) {
    stop("The title you have choosen has already been taken. Please choose another one before uploading again!")
  } else {
    postForm(upload_url(), title = dataset_title,  "datafile[file]" = upload_file(dataset), curl = curl)
    if(getCurlInfo(curl)$response.code != 200) {
      stop("Your upload failed. Try again later!")
    } else {
      id = title_to_dataset_id(dataset_title)
      message("Your data has been uploaded successfully!")
      message(paste0("You can find your dataset now under the id:", id))
      if(open_browser) {
	bef.goto.dataset_page(id = id)
      }
    }
  }
}

#' Update a dataset already on BEFdata portal
#'
#' This function updates a dataset from within R. Before you can do so you need to have
#' an account on the portal and set your credentials via the bef.options() function so
#' the authentication works. The function requires an ID of the dataset and either a path
#' to a CSV file or a data frame used to update the old dataset.
#'
#' Please be careful as this function overrides the old dataset (raw data) on the portal
#' completely.
#'
#' @param dataset The data you like to use for upload and to update the old dataset. This can
#' 	  be a either a path to a CSV file or a data frame.
#' @param dataset_id The ID of the dataset you like to update. You can get the ID from the URL
#'        of the dataset you like to update.
#' @param warn Is a warning that should prevent overriding dataset accidentially. Set this to
#'   	  FALSE if you know what you do.
#' @param curl You can pass in a curl handle to reduce memory footprint and to add own options
#' @param open_browser If this is set to true the page of the dataset is opened in the
#'        browser after successful upload. This defaults to FALSE.
#' @return Returns a status message whether the update was successful or not.
#' @examples \dontrun{
#'  dataset = data.frame(a = c(1,2,3,4), b = c("a","b","c","d"))
#'  bef.portal.update.dataset(dataset = dataset, dataset_id = 72)
#'       }
#' @import RCurl
#' @export bef.portal.update.dataset bef.update.dataset
#' @aliases bef.update.dataset

bef.portal.update.dataset <- bef.update.dataset <- function(dataset, dataset_id, warn = TRUE, open_browser = F, curl = getCurlHandle()) {
  if(warn) {
    stop("This function overrides the dataset on the portal! If you are sure set the function option warn = FALSE")
  }
  postForm(paste0("http://befdataproduction.biow.uni-leipzig.de/datasets/", dataset_id, "/update_workbook"),
	   user_credentials = bef.options("user_credentials"),
	   "datafile[file]" = upload_file(dataset),
	   curl = curl)
  if(getCurlInfo(curl)$response.code != 302) {
    stop("Your update failed. Check your access rights!")
  } else {
    if(open_browser) {
      bef.goto.dataset_page(id = dataset_id)
    }
    return(paste("Update of dataset with ID:", dataset_id, "successful!"))
  }
}

#' Attach a freeformat file to a dataset
#'
#' This function attaches a freeformat file to a dataset.
#'
#' @param id The ID of the dataset you like to attach to. You can get the ID from the URL
#'        of the dataset.
#' @param attachment Can be either a path to a file or an R data frame.
#' @param description An optional description to your attachmed freeformat file.
#' @param curl You can pass in a curl handle to reduce memory footprint and to add own options
#' @param open_browser If this is set to true the page of the dataset is opened in the
#'        browser after successful upload. This defaults to FALSE.
#' @return Returns a status message whether the update was successful or not.
#' @examples \dontrun{
#'   bef.portal.attach.to_dataset(id = 72, attachment = dataset, description = "This
#'                is a pure test description for the new attachment file", open_browser = T)
#'       }
#' @import RCurl
#' @export bef.portal.attach.to_dataset bef.attach.to_dataset
#' @aliases bef.attach.to_dataset

bef.portal.attach.to_dataset <- bef.attach.to_dataset <- function(id, attachment, description = "", open_browser = FALSE, curl = getCurlHandle()){
  postForm(paste0(bef.options("url"),"/files/freeformats"),
	   freeformattable_id = id,
	   freeformattable_type = "Dataset",
	   user_credentials = bef.options("user_credentials"),
	   "freeformat[file]" = upload_file(attachment),
	   "freeformat[description]" = description,
	   .opts = curlOptions(
			       referer="http://befdataproduction.biow.uni-leipzig.de",
			       useragent = "rbefdata"),
	   curl = curl)
  if(getCurlInfo(curl)$response.code != 302) {
    stop("Your attachment failed. Check your access rights!")
  } else {
    if(open_browser) {
      bef.goto.dataset_page(id = id)
    }
    return(paste("Attachment to dataset with ID:", id, "successful!"))
  }
}

#' Attach a freeformat file to a paper proposal
#'
#' This function attaches a freeformat file to a proposal
#'
#' @param id The ID of the proposal you like to attach to. You can get the ID from the URL
#'        of the proposal.
#' @param attachment Can be either a path to a file or an R data frame.
#' @param description An optional description to your attachmed freeformat file.
#' @param curl You can pass in a curl handle to reduce memory footprint and to add own options
#' @param is_paper Just to tell the portal if it is a published paper you are attaching as pdf
#'        for example. This parameter can be TRUE or FALSE and defaults to FALSE.
#' @param doi If the attachment has a persistent digital identifyer you can add it here.
#' @param open_browser If this is set to true the page of the proposal is opened in the
#'        browser after successful upload. This defaults to FALSE.
#' @return Returns a status message whether the update was successful or not.
#' @examples \dontrun{
#'   bef.portal.attach.to_proposal(id = 72, attachment = dataset, description = "This
#'                is a pure test description for the new attachment file", open_browser = T)
#'       }
#' @import RCurl
#' @export bef.portal.attach.to_proposal bef.attach.to_proposal
#' @aliases bef.attach.to_proposal

bef.portal.attach.to_proposal <- bef.attach.to_proposal <- function(id, attachment, description = "", is_paper = FALSE, doi = "",  open_browser = FALSE, curl = getCurlHandle()){
  postForm(paste0(bef.options("url"),"/files/freeformats"),
	   freeformattable_id = id,
	   freeformattable_type = "Paperproposal",
	   user_credentials = bef.options("user_credentials"),
	   "freeformat[uri]" = doi,
	   if(is_paper) {"freeformat[is_essential]" = 1} else {"freeformat[is_essential]" = 0},
	   "freeformat[file]" = upload_file(attachment),
	   "freeformat[description]" = description,
	   .opts = curlOptions(referer="http://befdataproduction.biow.uni-leipzig.de",
												 useragent = "rbefdata"),
												 curl = curl)

  if(getCurlInfo(curl)$response.code != 302) {
    stop("Your Attachment failed. Check your access rights!")
  } else {
    if(open_browser) {
      bef.goto.proposal_page(id = id)
    }
    return(paste("Attachment to proposal with ID:", id, "successful!"))
  }
}

#' Upload a category set
#'
#' This function is thought for data admin of a BEFdata portal. It uses the data haromization
#' features provided by the BEFdata portal. In order to do harmonize a datagroup you
#' download dthe category with the command below:
#' bef.portal.get.categories_for(datagroups_id = 22)
#'
#' @param datagroup_id Categories ID you like to harmonize. You can find this in the URL of
#' the category.
#' @param categories Either a data frame with the categories or a csv containing the categories
#' and merge ID and stuff downloaded from the portal.
#' @param curl Pass in a curl handle with own options or to reduce memory footprint.
#' @return Returns a status message whether the upload was successful or not.
#' @examples \dontrun{
#'   dataframe = my_categories_df
#'   bef.portal.upload.categories(categories_id = 22, categories = dataframe)
#'       }
#' @import RCurl
#' @export bef.portal.upload.categories bef.upload.categories
#' @aliases bef.upload.categories

bef.portal.upload.categories <- bef.upload.categories <- function(datagroup_id, categories, curl = getCurlHandle()) {
	# TODO: This is not working and might needs a proper login with session cookie
  postForm(datagroups_url(datagroups_id = datagroup_id , type = "upload"),
					 "csvfile[file]" = upload_file(categories),
					 .opts = curlOptions(referer="http://befdataproduction.biow.uni-leipzig.de",
															 useragent = "rbefdata"),
															 curl = curl)
  if(getCurlInfo(curl)$response.code != 302) {
    stop("Your category upload failed. Check your access rights!")
  } else {
    return(paste("Your Category upload was successful!"))
  }
}

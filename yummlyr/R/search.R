#' Search recipes on Yummly.com
#'
#' Query Yummly.com API to search for recipes with certain parameter. All parameters are optional and can be used in any combination.
#' The criteria you pass via the various parameters are combined with the AND operator (set conjunction). In other words, every recipe has to match the search phrase and satisfy the ingredient, cuisine, course, holiday, time, nutrition, and taste restrictions as described below.
#' If you specify a multi-word phrase to the q parameter, every word has to match something in each matching recipe:
#' @param search_words search phrase, can be supplied in from of vector of words
#' @param require_pictures set to \code{TRUE} if only to return recipes with photos
#' @param allowed_ingredient ingredient that all search results must include
#' @param excluded_ingredient ingredient that all search results should not contain
#' @param allowed_allergy only include recipes whose ingredients are allowed for that allergy
#' @param allowed_diet search results will only include recipes whose ingredients are allowed for that diet
#' @param allowed_cuisine search results will only include recipes with that cuisine
#' @param excluded_cuisine search results will only exclude recipes with that cuisine
#' @param allowed_course search results will only include recipes with that cuisine
#' @param excluded_course search results will only exclude recipes with that cuisine
#' @param allowed_holiday search results will only include recipes with that holiday
#' @param excluded_holiday search results will only exclude recipes with that holiday
#' @param max_total_time search for recipes that do not exceed a specified max total cook + prep time in seconds
#' @param max_results number of results to return
#' @param start start with specific result in search
#' @param nutrition set the range of allowed values for a given nutrition attribute (see below for the list of supported nutrition attributes) by setting a min and/or a max
#' @param flavor set the ranges for taste attributes (this corresponds to the taste sliders on the Yummly.com search page). The values of min and max are between 0 and 1.
#' @param facet_field facet counts for ingredient and diet. When this parameter is called, the response will include a facetCounts object that lists the matching diets or ingredients and how many results match each diet or ingredient.
#' @param app_id application ID
#' @param app_key application key
#' @note This function resembles search query to Yummly API
#' @references \itemize{
#'   \item Yummly Developer Guide \url{https://developer.yummly.com/documentation}
#' }
#' @examples
#' \dontrun{
#' # search for recipes with bacon
#' search_recipes("bacon")
#' 
#' # search for recipes with bacon that have pictures
#' search_recipes("bacon", require_pictures = TRUE)
#' 
#' # search for "Onion Soup" recipes which include garlic and cognac
#' search_recipes("Onion Soup", allowed_ingredient = c("garlic", "cognac"))
#' 
#' # search for "Onion Soup" recipes which do not include "onion soup mix" 
#' search_recipes("Onion Soup", excluded_ingredient = c("onion soup mix"))
#' 
#' # search for "Onion Soup" recipes that are Dairy-Free and Gluten-Free
#' search_recipes("bacon", allowed_allergy =c("Dairy-Free", "Gluten-Free"))
#' 
#' # search for "Onion Soup" recipes that are Pescetarian and Lacto vegetarian
#' search_recipes("bacon", allowed_diet =c("Pescetarian", "Lacto vegetarian")
#' 
#' # search for "Onion Soup" recipes that match American Cuisine
#' search_recipes("bacon", allowed_cuisine =c("American")
#' 
#' # exclude American recipes from a search for "Onion Soup"
#' search_recipes("bacon", excluded_cuisine =c("American")
#' 
#' # search for "Onion Soup" recipes that are Appetizers
#' search_recipes("bacon", allowed_course =c("Appetizers")
#' 
#' # exclude Appetizer recipes from a search for "Onion Soup" 
#' search_recipes("bacon", excluded_course =c("Appetizers")
#' 
#' # search for "Onion Soup" recipes for Thanksgiving 
#' search_recipes("bacon", allowed_holiday =c("Thanksgiving")
#' 
#' # exclude Thanksgiving recipes from a search for "Onion Soup"
#' search_recipes("bacon", excluded_holiday =c("Thanksgiving")
#' 
#' # if you want 20 recipes per page and want to see the second page of results
#' search_recipes("bacon", max_results = 20)
#' 
#' # if you want to start with position 20
#' search_recipes("bacon", start = 20)
#' 
#' # looking for recipes with a lot of Potassium, try setting a min of 3000 mg
#' # and a max of the Daily Suggested Value of 3500 mg
#' search_recipes("bacon", nutrition = list(Calcium=list(min=3, max=3.5)))
#' 
#' #  search for recipes which are very sweet but are not very spicy,
#' search_recipes("bacon", flavor = list(sweet=list(min=0.1, max=1)))
#' }
#' @export
search_recipes <- function(search_words, require_pictures,
                           allowed_ingredient, excluded_ingredient,
                           allowed_diet, allowed_allergy,
                           allowed_cuisine, excluded_cuisine,
                           allowed_course, excluded_course,
                           allowed_holiday, excluded_holiday,
                           max_total_time,
                           max_results, start,
                           nutrition, flavor,
                           facet_field,
                           app_id = auth_cache$APP_ID, app_key = auth_cache$APP_KEY) {
    if (!is.list(search_words) && !is.vector(search_words)) {
        stop("Wrong format of search lists, should be either list or vector")
    }
    if (is.null(app_id) || is.null(app_key)) {
        stop("APP_ID or APP_KEY is not set. Use setup_yummly_credentials or supply appropriate arguments")
    }
    # add search words
    search_words <- paste(search_words, collapse = " ")
    query <- sprintf("%s?_app_id=%s&_app_key=%s&q=%s", URL_SEARCH,
                     app_id, app_key, search_words)
    # add pictures requirement
    if (!missing(require_pictures)) {
        if (is.logical(require_pictures)) {
            query <- sprintf("%s&requirePictures=%s", query, tolower(require_pictures[1]))  
        } else {
            warning("require_pictures argument is not logical, it will be discarded")
        }
    }
    # add different parameters
    query <- add_argument(allowed_ingredient, "allowedIngredient", "ingredient", query)
    query <- add_argument(excluded_ingredient, "excludedIngredient", "ingredient", query)
    query <- add_argument(allowed_allergy, "allowedAllergy", "allergy", query)
    query <- add_argument(allowed_diet, "allowedDiet", "diet", query)
    query <- add_argument(allowed_cuisine, "allowedCuisine", "cuisine", query)
    query <- add_argument(excluded_cuisine, "excludedCuisine", "cuisine", query)
    query <- add_argument(allowed_course, "allowedCourse", "course", query)
    query <- add_argument(excluded_course, "excludedCourse", "course", query)
    query <- add_argument(allowed_holiday, "allowedHoliday", "holiday", query)
    query <- add_argument(excluded_holiday, "excludedHoliday", "holiday", query)
    # add maxTotalTime, maxResult, start
    if (!missing(max_total_time) && is.numeric(max_total_time)) {
        query <- paste(query, "&maxTotalTimeInSeconds=", max_total_time[1], sep="")
    } 
    if (!missing(max_results) && is.numeric(max_results)) {
        query <- paste(query, "&maxResult=", max_results[1], sep="")
    }
    if (!missing(start) && is.numeric(start)) {
        query <- paste(query, "&start=", start[1], sep="")
    }
    # add NUTRITION attribute
    if (!missing(nutrition)) {
        nutrition_search_values <- check_arguments(names(nutrition), "nutrition")
        incorrect_value <- which(!sapply(nutrition, function(x) is.numeric(x[[1]])))
        if (length(incorrect_value)) {
            stop(sprintf("For %s nutrition arguments, value parameter is not correct",
                         paste(names(nutrition)[incorrect_value]), collapse = ", "))
        } 
        incorrect_type <- which(!sapply(nutrition, function(x) names(x) %in% c("max", "min")))
        if (length(incorrect_type)) {
            stop(sprintf("For %s nutrition arguments, type parameter is not correct",
                         paste(names(nutrition)[incorrect_type]), collapse = ", "))
        }
        nutrition_argument <- sapply(names(nutrition), 
                                     function(x) {
                                         min <- sprintf("nutrition.%s.%s=%s",
                                                        nutrition_search_values[x],
                                                        "min",
                                                        nutrition[[x]]$min)
                                         max <- sprintf("nutrition.%s.%s=%s",
                                                        nutrition_search_values[x],
                                                        "max",
                                                        nutrition[[x]]$max)
                                         c(min, max)
                                     })
        query <- add_argument(unlist(nutrition_argument), argument_name = "", check = FALSE, query = query)
    }
    # add flavor attribute
    if (!missing(flavor)) {
        flavor_search_values <- check_arguments(names(flavor), "flavor")
        incorrect_value <- which(!sapply(flavor, function(x) is.numeric(x[[1]])))
        if (length(incorrect_value)) {
            stop(sprintf("For %s flavor arguments, value parameter is not correct",
                         paste(names(flavor)[incorrect_value]), collapse = ", "))
        } 
        incorrect_type <- which(!sapply(flavor, function(x) names(x) %in% c("max", "min")))
        if (length(incorrect_type)) {
            stop(sprintf("For %s flavor arguments, type parameter is not correct",
                         paste(names(flavor)[incorrect_type]), collapse = ", "))
        }
        incorrect_value <- which(!sapply(flavor, function(x) unlist(x) >= 0 && unlist(x) <= 1.0))
        if (length(incorrect_value)) {
            stop(sprintf("For %s flavor arguments, value parameter is not correct. It should be between 0 and 1",
                         paste(names(flavor)[incorrect_value]), collapse = ", "))
        }
        flavor_argument <- sapply(names(flavor), 
                                     function(x) {
                                         min <- sprintf("flavor.%s.%s=%s",
                                                        flavor_search_values[x],
                                                        "min",
                                                        flavor[[x]]$min)
                                         max <- sprintf("flavor.%s.%s=%s",
                                                        flavor_search_values[x],
                                                        "max",
                                                        flavor[[x]]$max)
                                         c(min, max)
                                     })
        query <- add_argument(unlist(flavor_argument), argument_name = "", check = FALSE, query = query)
    }
    if (!missing(facet_field)) {
        if (!all(facet_field %in% c("ingredient", "diet"))) {
            stop("Wrong facetField argument. Only diet and ingredient are supported")
        }
        query <- add_argument(facet_field, "facetField", "", query, check = FALSE)
    }
    content <- perform_query(utils::URLencode(query))
    jsonlite::fromJSON(content)
}

#' Add argument to a query
#' @param argument_values value of the argument
#' @param argument_name name of the argument
#' @param type argument type from metadata
#' @param query existing query to append to
#' @param check if to check again metadata values
add_argument <- function(argument_values, argument_name, type, query, check = TRUE) {
    if (missing(argument_values)) {
        return(query)
    }
    if (check) {
        argument_values <- check_arguments(argument_values, type)
    }
    arg <- prepare_array_parameter(argument_values, argument_name)
    paste(query, arg, sep = "&")       
}

#' Prepare search parameter
#'
#' Prepare search parameter from direction
#' @param param vector parameter to use
#' @param name name for parameter to use
prepare_array_parameter <- function(param, name) {
    if (name != "") {
        name <- paste(name, "[]=", sep="")
    } else {
        name <- paste(name, sep="")
    }
    paste(name, param, sep= "", collapse ="&")
}

#' Check ingredients
#'
#' Check ingredients list against predifined ingredients by Yummly
#' @param arguments arguments from metadata to check to check
#' @param type type of arguments from metadata
#' @note Predifined list is downloaded from Metadata Dictionaries
check_arguments <- function(arguments, type) {
    metadata <- metadata[[type]]
    field <- "description"
    available_arguments <- metadata[[field]]
    if (is.null(available_arguments)) {
        field <- "longDescription"
        available_arguments <- metadata[[field]]
    }
    result <- sapply(arguments, function(argument) {
        exact_match <- which(argument == available_arguments)
        possible_matches <- which(grepl(argument, available_arguments))
        if (length(exact_match) || length(possible_matches)) {
            if (length(exact_match)) {
                metadata[exact_match, ]$searchValue
            } else {
                warning(sprintf("Multiple arguments match %s (no exact match found), choosing %s",
                                argument, metadata[possible_matches[1], ][[field]]))
                metadata[possible_matches[1], ]$searchValue
            }
        } else {
            stop(sprintf("%s argument is not found (directly or loosely)", argument))
        }
    })
    result
}

# cfDatatype Class --------------------------------------------------------

#' @importFrom methods setClass
setClass(Class = "cfDatatype",
         slots = c(
           dt_name = "character",
           dt_type = "character",
           dt_sel_option_names = "list",
           dt_sel_combo_name = "character",
           dt_param = "character",
           dt_sel_option_params = "list",
           dt_selected_options = "list",
           dt_option_length = "numeric"
         ))

#' @importFrom methods setMethod
setMethod("initialize", "cfDatatype", function(.Object, dt_name, dt_type,
                                               dt_sel_option_names,
                                               dt_sel_combo_name,
                                               dt_param,
                                               dt_sel_option_params,
                                               dt_selected_options,
                                               dt_option_length){
  
  if (anyDuplicated(dt_param)){
    dt_name = dt_name[!duplicated(dt_param)]
    dt_type = dt_type[!duplicated(dt_param)]
    dt_sel_option_names = dt_sel_option_names[!duplicated(dt_param)]
    dt_sel_combo_name = dt_sel_combo_name[!duplicated(dt_param)]
    dt_param = dt_param[!duplicated(dt_param)]
    dt_sel_option_params = dt_sel_option_params[!duplicated(dt_param)]
    dt_selected_options = dt_selected_options[!duplicated(dt_param)]
    dt_option_length = dt_option_length[!duplicated(dt_param)]
  }
  
  match_dt = order(match(dt_name, unique(dt_name)))
  dt_name = dt_name[match_dt]
  dt_type = dt_type[match_dt]
  dt_sel_option_names = dt_sel_option_names[match_dt]
  dt_sel_combo_name = dt_sel_combo_name[match_dt]
  dt_param = dt_param[match_dt]
  dt_sel_option_params = dt_sel_option_params[match_dt]
  dt_selected_options = dt_selected_options[match_dt]
  dt_option_length = dt_option_length[match_dt]
  
  lengths = dt_option_length
  n.dt = length(lengths)
  
  selections = dt_selected_options
  
  lengths_list = lapply(lengths, seq, from = 1)
  prm_list = mapply("+", lengths_list, c(0, cumsum(head(lengths, -1))), 
                    SIMPLIFY = FALSE)
  
  dt_param = paste(dt_param,
                   sapply(prm_list, paste0, collapse = ","), 
                   sep = ",")
  names(dt_param) = paste0("dt", seq_along(lengths))
    
  prm_sel_list = mapply("[", prm_list, selections, SIMPLIFY = FALSE)
  prm_names = lapply(prm_sel_list, function(x) paste0("prm", x))
  
  dt_sel_option_params = mapply("names<-", dt_sel_option_params, prm_names,
                                SIMPLIFY = FALSE)
  
  .Object@dt_name = dt_name
  .Object@dt_type = dt_type
  .Object@dt_sel_option_names = dt_sel_option_names
  .Object@dt_sel_combo_name = dt_sel_combo_name
  .Object@dt_param = dt_param
  .Object@dt_sel_option_params = dt_sel_option_params
  .Object@dt_selected_options = dt_selected_options
  .Object@dt_option_length = dt_option_length
  
  return(.Object)  
})

# Internals ---------------------------------------------------------------

# Return datatype information for a given stage of the selection process. 
#
# This function returns the name of the chosen datatype and the href to link it 
# to the next branch of the datatype selection tree.
#
# doc: the XML to extract the information for the datatypes
# ...: passed to menu
#' @importFrom XML xmlValue xmlGetAttr
#' @importFrom utils menu
dt_href = function(doc, ...){
  # "cloud \n          cover" --> "cloud cover"
  choices = gsub("\\n", "", gsub(" {2,}", "", sapply(doc, xmlValue)))
  hrefs = sapply(doc, xmlGetAttr, "href")
  dt = menu(choices, ...)
  if (dt)
    c(hrefs[dt], choices[dt])
}

# Return the href and datatype name from the first stage datatype selection
# 
# This function is only intended for use within the cf_datatype function.
# The arguments are passed from the cf_datatype arguments to this function.
#
# selection: passed from the select_1 argument
# g        : logical passed to the graphics argument of the menu function
#' @importFrom selectr querySelectorAll
#' @importFrom XML htmlParse xmlGetAttr xmlValue
first_stage_selection = function(selection, g, iter){
  domain = "http://cliflo.niwa.co.nz/pls/niwp/"
  full_path = paste0(domain, "wgenf.choose_datatype?cat=cat1")
  
  datatypes_xml = querySelectorAll(htmlParse(full_path),
                                   "table.header td.popup a.top")
  
  if (!is.na(selection) && selection > 9){
    stop(paste("the first selection can only be between 1 and 9 for datatype", 
               iter), call. = FALSE)
  }
  
  if (is.na(selection)){
    dt_href(datatypes_xml, graphics = g, 
            title = "Daily and Hourly Observations")
  } else {
    dt_name = sapply(datatypes_xml, xmlValue)[selection]
    c(sapply(datatypes_xml, xmlGetAttr, "href")[selection], dt_name)
  }
}

# Return the href and datatype type from the second stage datatype selection
# 
# This function is only intended for use within the cf_datatype function. The 
# arguments are passed from the cf_datatype arguments to this function.
#
# href_1    : href from the first stage selection
# selection : second stage selection passed from the select_2 argument
# dt_name   : the name of the datatype from the first stage selection. This is
#             used for menu titles and warnings
# g         : logical passed to the graphics argument of the menu function
#' @importFrom selectr querySelectorAll
#' @importFrom XML htmlParse xmlGetAttr xmlValue
second_stage_selection = function(href_1, selection, dt_name, g){
  domain = "http://cliflo.niwa.co.nz/pls/niwp/"
  full_path = paste0(domain, href_1)
  datatypes_xml = querySelectorAll(htmlParse(full_path), "a.dt")
  
  if (is.na(selection)){
    gsub("\\n", "", gsub(" {2,}", "", dt_href(datatypes_xml, g, dt_name)))
  } else {
    if (selection > length(datatypes_xml))
      stop(paste("second selection (select_2) is out of range for", 
                 dt_name), call. = FALSE)
    dt_type = gsub("\\n", "", gsub(" {2,}", "", 
                                   sapply(datatypes_xml, xmlValue)[selection]))
    c(sapply(datatypes_xml, xmlGetAttr, "href")[selection], dt_type)
  }
}

# Choose datatype options interactively
#
# This function is only used when the user does not supply check box options for
# a given datatype and is only intended for use within the option_selections
# internal function.
#
# datatype_name    : the name of the datatype used for the title
# datatype_options : the possible options for the given datatype to be displayed 
#                    in the menu
# g                : logical passed to the graphics argument of the menu function
#' @importFrom utils menu
choose_dt_options = function(datatype_name, datatype_options, g){
  selected_options = menu(datatype_options, g, paste(datatype_name, "options"))
  
  if (length(datatype_options) > 1){
    again = menu(c("yes", "no"), g, "Choose another option?")
    finished = FALSE
    
    while(again == 1 && !finished){
      selected_options = c(selected_options,
                           menu(datatype_options, g, 
                                paste(datatype_name, "options")))
      finished = length(unique(selected_options)) == length(datatype_options)
      if (!finished)
        again = menu(c("yes", "no"), g, "Choose another option?")
    }
  }
  
  selected_options
}

# Return all other options to produce an instance of the cfDatatype class
# 
# This function is only intended for use within the cf_datatype function. The 
# arguments are passed from the cf_datatype arguments to this function.
#
# href_2          : href from the second stage selection
# selection_check : the users check box selections passed from check_box
# selection_combo : the users combo box selections passed from combo_box
# dt_type         : the datatype from the second stage selection. This is
#                   used for menu titles and warnings
# g               : logical passed to the graphics argument of the menu function
#' @importFrom selectr querySelectorAll
#' @importFrom XML htmlParse xmlGetAttr xmlValue xmlSApply
#' @importFrom utils menu
option_selections = function(href_2, selection_check, selection_combo, 
                             dt_type, g){
  selection_check = unique(selection_check)
  domain = "http://cliflo.niwa.co.nz/pls/niwp/"
  full_path = paste0(domain, href_2)
  dt_options_xml = querySelectorAll(htmlParse(full_path), 
                                    "td.selected table tr td.selected")
  dt_params_xml = querySelectorAll(htmlParse(full_path), 
                                   "td.selected table tr td input")
  dt_param_values = xmlSApply(dt_params_xml, xmlGetAttr, "value")
  dt_combo_xml = querySelectorAll(htmlParse(full_path), 
                                  "td.selected table tr td select option")
  
  if (length(dt_combo_xml) == 0)
    dt_options = xmlSApply(dt_options_xml, xmlValue)
  else{
    dt_options_inc_combo = xmlSApply(dt_options_xml, xmlValue)
    dt_options_inc_combo = dt_options_inc_combo[dt_options_inc_combo != ""]
    combo_name = tail(dt_options_inc_combo, 1)
    dt_options = head(dt_options_inc_combo, -1)
    dt_combo_param_names = xmlSApply(dt_combo_xml, xmlValue)
    dt_combo_param_values = xmlSApply(dt_combo_xml, xmlGetAttr, "value")
  }
  
  if (any(is.na(selection_check)))
    selected_options = choose_dt_options(dt_type, dt_options, g)
  else{
    if (length(selection_check) > length(dt_options))
      stop(paste("the number of check box options is too many for datatype", 
                 dt_type), call. = FALSE)
    
    if (any(selection_check > length(dt_options)))
      stop(paste("the check box options for datatype", dt_type, 
                 "must be between 1 and", length(dt_options)), call. = FALSE)
    selected_options = selection_check
  }
  
  selected_params = dt_param_values[selected_options]
  selected_param_names = dt_options[selected_options]
  
  combo_names = NA
  if (length(dt_combo_xml) != 0){
    if (is.na(selection_combo))
      combo_selected = menu(dt_combo_param_names, g, combo_name)
    else{
      if (length(selection_combo) != 1)
        stop(paste("you can only choose one combo box option for datatype", 
                   dt_type), call. = FALSE)
      
      if (selection_combo > length(dt_combo_param_names))
        stop(paste("the combo box option for datatype", dt_type, 
                   "must be between 1 and", length(dt_combo_param_names)), 
             call. = FALSE)
      combo_selected = selection_combo
    }
    
    selected_params = c(selected_params, dt_combo_param_values[combo_selected])
    selected_options = c(selected_options, length(dt_options_xml))
    combo_names = dt_combo_param_names[combo_selected]
  } else
    if (!is.na(selection_combo))
      message(paste("combo options are not required for", dt_type))
  
  dt_param = gsub("wgenf.genform1\\?cdt=|\\&.*", "", href_2)
  
  list(dt_param, 
       selected_params, selected_param_names, combo_names,
       selected_options, length(dt_options_xml))
}

# Add selected datatypes to the curl session. 
#
# Equivalent to updating the cliflo page. Adds each datatype 
# to the session and saves the cookies in the temporary directory for future 
# use.
#
# object : a cfDatatype object
# user   : a cfUser object
#' @importFrom RCurl getCurlHandle getForm
cf_update_dt = function(object, user = cf_user()){
  cookies = file.path(tempdir(), user@username)
  curl = getCurlHandle(followlocation = TRUE,
                       timeout = 100, 
                       useragent = 
                         paste("clifro", R.Version()$version.string),
                       cookiefile = cookies, 
                       cookiejar = cookies)
  
  all_dt_params = c(object@dt_param, unlist(object@dt_sel_option_params))
  
  postForm("http://cliflo.niwa.co.nz/pls/niwp/wgenf.genform1_proc",
           cselect = "wgenf.genform1?fset=defdtype",
           auswahl = "wgenf.genform1?fset=defagent",
           agents = "",
           dateauswahl = "wgenf.genform1?fset=defdate",
           date1_1="2014",
           date1_2="05",
           date1_3="25",
           date1_4="00",
           date2_1="2014",
           date2_2="05",
           date2_3="28",
           date2_4="00",
           formatselection = "wgenf.genform1?fset=deffmt",
           TSselection = "NZST",
           dateformat = "0",
           Splitdate = "N",
           mimeselection = "htmltable",
           cstn_id = "A",
           cdata_order = "DS",
           .params = all_dt_params,
           curl = curl)
}


# cfDatatype constructor --------------------------------------------------

#' The Clifro Datatype Object
#' 
#' Create a \code{cfDatatype} object by selecting one or more CliFlo datatypes 
#' to build the \pkg{clifro} query.
#' 
#' An object inheriting from the \code{\link{cfDatatype}} class is created by 
#' the constructor function \code{\link{cf_datatype}}. The function allows the 
#' user to choose datatype(s) interactively (if no arguments are given), or to 
#' create datatypes programatically if the tree menu nodes are known a priori 
#' (see examples). This function uses the same nodes, check box and combo box 
#' options as CliFlo and can be viewed at the 
#' \href{http://cliflo.niwa.co.nz/pls/niwp/wgenf.choose_datatype?cat=cat1}{datatype selection page}.
#' 
#' @param select_1 a numeric vector of first node selections
#' @param select_2 a numeric vector of second node selections
#' @param check_box a list containing the check box selections
#' @param combo_box a numeric vector containing the combo box selection 
#' (if applicable)
#' @param graphics a logical indicating whether a graphics menu should be used,
#' if available
#' 
#' @note For the 'public' user (see examples) only the Reefton Ews station data 
#' is available. 
#' 
#' @note Currently clifro does not support datatypes from the special datasets
#' (Ten minute, Tier2, Virutal Climate, Lysimeter) or upper air measurements
#' from radiosondes and wind radar.
#' 
#' @importFrom methods new
#' @name cfDatatype-class
#' @rdname cfDatatype-class
#' @aliases cfDatatype
#' @export
#' @return \code{cfDatatype} object
#' @seealso \code{\link{cf_user}} to create a \pkg{clifro} user, 
#'   \code{\link{cf_station}} to choose the CliFlo stations and 
#'   \code{vignette("choose-datatype")} for help choosing \code{cfDatatype}s.
#' @examples
#' \dontrun{
#' # Select the surface wind datatype manually (unknown tree nodes)
#' hourly.wind.dt = cf_datatype()
#' #  2  --> Datatype:        Wind
#' #  1  --> Datatype 2:      Surface Wind
#' #  2  --> Options:         Hourly Wind
#' # (2) --> Another option:  No
#' #  3  --> Units:           Knots
#' hourly.wind.dt
#'
#' # Or select the datatype programatically (using the selections seen above)
#' hourly.wind.dt = cf_datatype(2, 1, 2, 3)
#' hourly.wind.dt
#' }
cf_datatype = function(select_1 = NA, 
                       select_2 = NA, 
                       check_box = NA,
                       combo_box = NA, 
                       graphics = FALSE){
  
  select_1 = unlist(select_1)
  select_2 = unlist(select_2)
  combo_box = unlist(combo_box)
  
  ## Based on modified code from the R Core Team - with thanks
  is_wholenumber = function(x){
    if (is.na(x))
      TRUE
    else{
      if (!is.numeric(x))
        FALSE
      else
        abs(x - round(x)) < .Machine$double.eps^0.5
    }
  }
  
  if (!all(sapply(select_1, is_wholenumber),
           sapply(select_2, is_wholenumber),
           unlist(lapply(check_box, function(x) sapply(x, is_wholenumber))),
           sapply(combo_box, is_wholenumber))
  )
    stop("arguments must be NA's or integers")
  
  if (!is.list(check_box))
    check_box = list(check_box)
  
  ## Check argument lengths coincide with number of datatypes
  arg.lengths = c(length(select_1), length(select_2), length(combo_box))
  unexp.arg.length = arg.lengths != max(arg.lengths)
  
  if (any(unexp.arg.length)){
    which.unexp = which(unexp.arg.length)
    stop(paste(names(match.call())[which.unexp[1] + 1], 
               "argument has unexpected length", 
               arg.lengths[which.unexp[1]]))
  }
  
  n.dt = length(select_1)
  dt_name = dt_type = dt_param = dt_sel_combo_name = character(n.dt)
  dt_option_length = numeric(n.dt)
  dt_selected_options = dt_sel_option_names = 
    dt_sel_option_params = vector("list", n.dt)
  
  for (i in seq_along(select_1)){
    
    if (is.na(select_1[i]) && any(!is.na(select_2[i]),
                                  !is.na(check_box[[i]]),
                                  !is.na(combo_box[i])))
      stop(paste("first selection must be known before other options are",
                 "chosen for datatype", i))
    
    if (any(is.na(select_1[i]), is.na(select_2[i])) && 
          any(!is.na(check_box[[i]]), !is.na(combo_box[i])))
      stop(paste("first and second selections must be known before other",
                 "options are chosen for datatype", i))
    
    if (any(is.na(select_1[i]), is.na(select_2[i]), 
            is.na(check_box[[i]])) && !is.na(combo_box[i]))
      stop(
        paste("all selections must be known before combo box options are",
              "chosen for datatype", i))
    
    href_1 = first_stage_selection(select_1[i], graphics, i)
    dt_name[i] = href_1[2]
    
    href_2 = second_stage_selection(href_1[1], select_2[i], 
                                    dt_name[i], graphics)
    dt_type[i] = href_2[2]
    
    options_list = option_selections(href_2[1], check_box[[i]], combo_box[i], 
                                     dt_type[i], graphics)
    
    dt_sel_option_names[[i]] = options_list[[3]]
    dt_sel_combo_name[i] = options_list[[4]]
    dt_param[i] = options_list [[1]]
    dt_sel_option_params[[i]] = options_list[[2]]
    dt_selected_options[[i]] = options_list[[5]]
    dt_option_length[i] = options_list[[6]]
  }
  
  new("cfDatatype", 
      dt_name = dt_name,
      dt_type = dt_type,
      dt_sel_option_names = dt_sel_option_names,
      dt_sel_combo_name = dt_sel_combo_name,
      dt_param = dt_param,
      dt_sel_option_params = dt_sel_option_params,
      dt_selected_options = dt_selected_options,
      dt_option_length = dt_option_length)
}


# Methods -----------------------------------------------------------------

#' @importFrom methods new setMethod
#' @rdname Extract
#' @aliases [,cfDatatype,ANY,missing,missing
setMethod("[", 
          signature(x = "cfDatatype",
                    i = "ANY",
                    j = "missing",
                    drop = "missing"),
          function(x, i, j, drop){
            
            dt_param = sapply(x@dt_param, function(x) strsplit(x, ",")[[1]][1])
            
            new("cfDatatype", 
                dt_name = x@dt_name[i], 
                dt_type = x@dt_type[i], 
                dt_sel_option_names = x@dt_sel_option_names[i], 
                dt_sel_combo_name = x@dt_sel_combo_name[i], 
                dt_param = dt_param[i],
                dt_sel_option_params = x@dt_sel_option_params[i], 
                dt_selected_options = x@dt_selected_options[i],
                dt_option_length = x@dt_option_length[i])
          })

#' Arithmetic Operators for Clifro Objects
#' 
#' This operator allows you to add more datatypes or stations to 
#' \code{cfDatatype} and \code{cfStation} objects respectively.
#' 
#' @param e1 a \code{cfDatatype} or \code{cfStation} object
#' @param e2 an object matching the class of e1
#' 
#' @rdname clifroAdd
#' @aliases +,cfDatatype,cfDatatype-method
#' @importFrom methods new setMethod
setMethod("+", signature(e1 = "cfDatatype",
                         e2 = "cfDatatype"),
          function(e1, e2){
            
            dt_param = sapply(c(e1@dt_param, e2@dt_param), 
                              function(x) strsplit(x, ",")[[1]][1])
            
            new("cfDatatype",
                dt_name = c(e1@dt_name, e2@dt_name), 
                dt_type = c(e1@dt_type, e2@dt_type), 
                dt_sel_option_names = 
                  c(e1@dt_sel_option_names, e2@dt_sel_option_names), 
                dt_sel_combo_name = 
                  c(e1@dt_sel_combo_name, e2@dt_sel_combo_name), 
                dt_param = dt_param,
                dt_sel_option_params = 
                  c(e1@dt_sel_option_params, e2@dt_sel_option_params), 
                dt_selected_options = 
                  c(e1@dt_selected_options, e2@dt_selected_options),
                dt_option_length = c(e1@dt_option_length, e2@dt_option_length))
          })

# Show method
#' @importFrom methods new setMethod show
setMethod("show", signature(object = "cfDatatype"), function(object){
  first_sel_names = object@dt_name
  second_sel_names = object@dt_type
  
  check_select = sapply(object@dt_sel_option_names, 
                        function(x) paste0("[", paste(x, collapse = ","), "]"))
    
  combo_select = object@dt_sel_combo_name
  if (any(is.na(combo_select)))
    combo_select = replace(combo_select, which(is.na(combo_select)), "")
  show(data.frame(dt.name = first_sel_names,
                  dt.type = second_sel_names, 
                  dt.options = check_select, 
                  dt.combo = combo_select, 
                  row.names = names(object@dt_param)))
})
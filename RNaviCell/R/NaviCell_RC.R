#    RNaviCell package
#
#    Copyright (C) {2015} {Institut Curie, 26 rue d'Ulm, 75005 Paris}
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
#    USA

#' NaviCell reference class
#'
#' NaviCell (https://navicell.curie.fr) is a web-based environment for
#' browsing, commenting and analyzing very large biological molecular networks
#' using Google Maps and for visualizing 'omics' data on top of the network
#' maps.
#'
#' NaviCell can also act as a server allowing to be remotely controlled through
#' a REST API. A python and a R language bindings have been developped on top
#' of the REST API to hide technical details and to provide users and
#' programmers a friendly interface. A Java binding has been initiated. 
#' 
#' This is the R binding implementation. For more information about the
#' NaviCell Web Service and a tutorial on how to use it, see
#' https://navicell.curie.fr/pages/nav_web_service.html and
#' https://github.com/eb00/RNaviCell.
#'
#' @import RCurl
#' @import RJSONIO
#' @importFrom methods new
#' @export NaviCell
#' @examples \dontrun{
#'	### Opens a communication with web service, build does not finish if example is tested
#'	file<-system.file("extdata", "script.R", package = "RNaviCell")
#'	source("file")
#'	}
NaviCell <- setRefClass(
    # class name
    "NaviCell",

    # Define the fields 
    fields = list( 
        proxy_url = "character",
        map_url = "character",
        msg_id = "numeric",
        session_id = "character",
        hugo_list = "vector",
        packsize = "numeric"
    ),

    # Set default values
    methods = list(
        initialize = function(...) {
            proxy_url <<- "https://navicell.curie.fr/cgi-bin/nv_proxy.php"
            map_url <<- "https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php"
            msg_id <<- 1000
            session_id <<- ""
            packsize <<- 500000
        }
    )
)


#------------------------------------------------------------------------------
#
#  Session and utility functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    incMessageId = function(...) {
    "Increase message ID counter." 
        msg_id <<- msg_id + 1
    }
)

NaviCell$methods(
    generateSessionId =  function(...) {
    "Generate a session ID."    
        .self$incMessageId()
        response = postForm(.self$proxy_url, style = 'POST', id = "1", perform = "genid", msg_id = .self$msg_id, mode = "session", .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        if (response != "") {
            response <- .self$formatResponse(response)
            .self$session_id <- response
        }
    }
)

NaviCell$methods(
    serverIsReady =  function(...) {
    "Test if NaviCell server is ready (internal utility)."    
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_is_ready')
        str_data <- .self$makeData(.self$formatJson(list_param))
        .self$incMessageId()
        response <- postForm(.self$proxy_url, style='POST', id = .self$session_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
        return(fromJSON(response)$data)
    }
)

NaviCell$methods(
    waitForReady =  function(...) {
    "Wait until NaviCell server is ready (internal utility)."  
        for (i in 1:50) {
            if (.self$serverIsReady() == TRUE) {
                break
            }
            else {
                message("wait for NaviCell server to be ready..")
                Sys.sleep(1)
            }
        }
    }
)

NaviCell$methods(
    isImported =  function(...) {
    "Test if data table is imported (internal utility)."    
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_is_imported')
        str_data <- .self$makeData(.self$formatJson(list_param))
        .self$incMessageId()
        response <- postForm(.self$proxy_url, style='POST', id = .self$session_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
        return(fromJSON(response)$data)
    }
)

NaviCell$methods(
    waitForImported =  function(...) {
    "Wait until data is imported (internal utility)."  
        for (i in 1:50) {
            if (.self$isImported() == TRUE) {
                message("data imported.")
                break
            }
            else {
                message("waiting for data to be imported...")
                Sys.sleep(0.5)
            }
        }
    }
)


NaviCell$methods(
    launchBrowser = function(...) {
    "Launch client browser and points to the default NaviCell map."
        .self$incMessageId()
        if (.self$session_id == "") {
            .self$generateSessionId()
        }
        url <- paste(.self$map_url, '?id=', .self$session_id, sep = '')
        browseURL(url)
        .self$waitForReady()
    }
)

NaviCell$methods(
    listSessions = function(...) {
    "List all NaviCell server sessions."
        response <- postForm(.self$proxy_url, style='POST', id="1", perform="list", msg_id = .self$msg_id, mode="session", .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        message(response)
    }
)

NaviCell$methods(
    attachSession = function(session_id) {
    "Attach a NaviCell server session ID."
        if (.self$session_id != "") {
            warning("Session id already set.")
            return()
        }
        .self$incMessageId()
        # check session id on NaviCell server
        response <- postForm(.self$proxy_url, style='POST', id = session_id, msg_id = .self$msg_id, mode = 'session', perform = 'check', .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        if (response == "ok") {
            .self$session_id <- session_id
        }
        else {
            warning("Wrong session id.")
        }
    }
)

NaviCell$methods(
    attachLastSession = function(...) {
    "Attach NaviCell handle to the last existing NaviCell Web Service session."
        response <- postForm(.self$proxy_url, style='POST', id = '1', msg_id = .self$msg_id, mode = 'session', perform = 'get', which = '@@', .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        message(response)
        .self$attachSession(response)
    }
)

NaviCell$methods(
    formatResponse = function(response) {
    "Format response obtained from the RCurl 'postForm' command (internal utility)."
        ret = ''
        if (class(response) == 'raw') {
            ret <- rawToChar(response)
        }
        else if (class(response) == "character") {
            ret <- response[1]
        }
        return(ret)
    }
)

NaviCell$methods(
    makeData = function(json_string) {
    "Create NaviCell server command string from a list of parameters (internal utility)."
         ret <- paste("@COMMAND ", json_string, sep = "")
         return(ret)
    }
)

NaviCell$methods(
    formatJson = function(list_param) {
    "Format list of parameters to NaviCell server compatible JSON format (internal utility)."
        data <- toJSON(list_param)
        # remove unnecessary characters 
        data <- gsub("\n", '', data)
        data <- gsub(" ", "", data)
        return(data)
    }
)

NaviCell$methods(
    file2dataString = function(fileName) {
    "Load the content of a text file as tab-delimited string. Convert to NaviCell compatible format."
        data_string <- NULL
        data_string <- paste(readLines(fileName, warn=F),collapse='\n')
        if (!is.null(data_string)) {
            if (substr(data_string, nchar(data_string), nchar(data_string)) != "\n") {
                data_string <- paste(data_string, '\n', sep="")
            }
            data_string <- paste("@DATA\n", data_string, sep="")
        }
        return(data_string)
    }
)


#------------------------------------------------------------------------------
#
#  Navigation and Zooming functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    setZoom = function(zoom_level) {
    "Set a given zoom level on associated NaviCell map in browser. zoom_level = integer value."
        .self$incMessageId()
        list_param <- list(module='', args = array(zoom_level), msg_id = .self$msg_id, action = 'nv_set_zoom')
        str_data <- .self$makeData(.self$formatJson(list_param))
        .self$incMessageId()
        response <- postForm(.self$proxy_url, style='POST', id = .self$session_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    setMapCenter = function(location) {
    "Set the relative position of the map center. location = 'MAP_CENTER' or 'MAP_EAST' or 'MAP_SOUTH' or MAP_NORTH' or 'MAP_SOUTH_WEST' or 'MAP_SOUTH_EAST' or 'MAP_NORTH_EAST'."
        .self$incMessageId()
        list_param <- list(module='', args = array(location), msg_id = .self$msg_id, action = 'nv_set_center')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    setMapCenterAbsolute = function(pos_x, pos_y) {
    "Set the absolute position of the map center. x = x coordinate (integer), y = y coordinate (integer)." 
        .self$incMessageId()
        list_param <- list(module='', args = list('ABSOLUTE', pos_x, pos_y), msg_id = .self$msg_id, action = 'nv_set_center')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    moveMapCenter = function(x, y) {
    "Move the map center (relative). x = x coordinate (integer), y = y coordinate (integer)." 
        .self$incMessageId()
        list_param <- list(module='', args = list('RELATIVE', x, y), msg_id = .self$msg_id, action = 'nv_set_center')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)


#------------------------------------------------------------------------------
#
# Entity selection functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    selectEntity = function(entity) {
    "Select an entity on the map. entity = entity's name (string)" 
        .self$incMessageId()
        list_param <- list(module='', args = array(entity), msg_id = .self$msg_id, action = 'nv_find_entities')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    findEntities = function(entity, bubble) {
    "Find one or more entities on the map. entity = entity's name pattern (string), bubble = TRUE or FALSE." 
        .self$incMessageId()
        list_param <- list(module='', args = array(c(entity, bubble)), msg_id = .self$msg_id, action = 'nv_find_entities')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    uncheckAllEntities = function(...) {
    "Uncheck all entities on the map."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_uncheck_all_entities')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)

NaviCell$methods(
    unhighlightAllEntities = function(...) {
    "Uncheck all entities on the map."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_unhighlight_all_entities')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #print(response)
    }
)


#------------------------------------------------------------------------------
#
#  Get info from NaviCell server functions 
#
#------------------------------------------------------------------------------


NaviCell$methods(
    getHugoList = function(...) {
    "Get the list of the HUGO gene symbols for the current map (the list is stored in the object field hugo_list."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_hugo_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        #message(response)
        response <- fromJSON(response)
        .self$hugo_list <- response$data
        return(response$data)
    }
)

NaviCell$methods(
    getBiotypeList = function(...) {
    "Return the list of biotypes understood by NaviCell Web Service."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_biotype_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        response <- fromJSON(response)
        return(response$data)
    }
)

NaviCell$methods(
    getModuleList = function(...) {
    "Return the module list of the current NaviCell map."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_module_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        response <- fromJSON(response)
        return(response$data)
    }
)

NaviCell$methods(
    getImportedDatatables = function(...) {
    "Return the list of datatables imported in the current NaviCell session."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_datatable_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        response <- fromJSON(response)
        return(response$data)
    }
)

NaviCell$methods(
    getImportedSamples = function(...) {
    "Return the list of samples from all the datatables imported in the current NaviCell session."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_datatable_sample_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        response <- fromJSON(response)
        return(response$data)
    }
)

NaviCell$methods(
    getImportedGenes = function(...) {
    "Return the list of genes from all the datatables imported in the current NaviCell session."
        .self$incMessageId()
        list_param <- list(module='', args = array(), msg_id = .self$msg_id, action = 'nv_get_datatable_gene_list')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        response <- .self$formatResponse(response)
        response <- fromJSON(response)
        return(response$data)
    }
)

#------------------------------------------------------------------------------
#
#  Datatable import functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    importDatatable = function(datatable_biotype, datatable_name, mat) {
    "Import a datatable (matrix) in the current map session."

        # check if the field hugo_list is set
        if (length(.self$hugo_list) == 0) {
            hl <- .self$getHugoList()
        }

        # filter matrix on hugo_list
        # abort if there is no overlap
        rownames(mat) %in% .self$hugo_list -> idx 
        if (sum(idx) < 1) {
            warning("Error: no overlap between map and matrix HUGO gene symbols.")
            return()
        }

        data_string <- NULL
        # case 1: simple gene list
        if (dim(mat)[2] == 0) {
            # convert to list
            gene_list <- rownames(mat)
            sel <- gene_list[idx]
            data_string <- .self$geneList2string(sel) 
        }

        # case 2: one or more columns in the matrix
        if (dim(mat)[2] > 0) {
            # select rows on hugo_list
            # watch out: if matrix has 1 col, return type is vector, not matrix, so cast the return to matrix.
            # watch out: when matrix as only one col. colname is lost with subselect, so set it back
            mat_select <- as.matrix(mat[idx,])
            if (ncol(mat) == 1) {
                colnames(mat_select) <- colnames(mat)
            }
            data_string <- .self$matrix2string(mat_select)
        }

        if (!is.null(data_string)) {
            if (nchar(data_string) < .self$packsize) {
                #print(data_string)
                .self$incMessageId()
                list_param <- list(module='', args = list(datatable_biotype, datatable_name, "", data_string, emptyNamedList), msg_id = .self$msg_id, action = 'nv_import_datatables')
                str_data <- .self$makeData(.self$formatJson(list_param))

                .self$incMessageId()
                response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
                .self$waitForImported()
                #print(.self$formatResponse(response))
            }
            else {
                #print("size > packsize")
                list_param <- list(module='', args = list(datatable_biotype, datatable_name, "", data_string, emptyNamedList), action = 'nv_import_datatables')
                fill_cmd <- .self$makeData(.self$formatJson(list_param))
                .self$sendBigData(fill_cmd)
                .self$waitForImported()
            }
        }
    }
)


NaviCell$methods(
    sendBigData = function(fill_cmd) {
        "slice data in packets big data string and send it to server (internal utility)"
        cmd_len <- nchar(fill_cmd)
        cmd_packcount <- as.integer(cmd_len / .self$packsize) + 1

        # slice data in packets and send them to server
        for (i in 0:(cmd_packcount-1)) {
            cmd_packnum <- i+1
            
            stop <- (i+1) * .self$packsize
            start <- 0
            if (i > 0) {
                start <- (i * packsize) + 1
            }
            if (stop > cmd_len) {
                stop <- cmd_len - 1
            }
            start <- start + 1
            stop <- stop + 1

            response <- postForm(.self$proxy_url, style = 'POST', perform = "filling", data = substr(fill_cmd, start, stop), id = .self$session_id, packcount = cmd_packcount, packnum = cmd_packnum, mode='cli2srv', .opts=curlOptions(ssl.verifypeer=F)) 
            #print(.self$formatResponse(response))
        }

        # end message to trigger data re-composition on server side
        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', data="@@", id = .self$session_id, perform = "send_and_rcv", packcount = cmd_packcount, msg_id = .self$msg_id, mode='cli2srv', .opts=curlOptions(ssl.verifypeer=F)) 

    }
)


NaviCell$methods(
    readDatatable = function(fileName) {
    "Read a data file and create an R matrix. Returns a matrix object."
        mat <- as.matrix(read.table(fileName, header=T, row.names=1))
        return(mat)
    }
)


NaviCell$methods(
    matrix2string = function(mat) {
    "Convert an R matrix object to a formatted string (internal utility)."
        header = ""
        if (nrow(mat) == 1) {
            header <- paste(colnames(mat), sep="") 
        }
        else {
            header <- paste(colnames(mat), collapse='\t', sep="") 
        }

        string <- paste('@DATA\ngenes\t', header, sep="")
        string <- paste(string, '\n', sep="")
        nb_row = nrow(mat) 
        for (row in 1:nb_row) {
            gene_name <- paste(rownames(mat)[row], '\t', sep="")  
            row_string = ""
            if (nb_row == 1) {
                row_string <- paste(mat[row], sep="")
            }
            else {
                row_string <- paste(mat[row,], collapse='\t', sep="")
            }
            row_string <- paste(row_string, '\n', sep="")
            string <- paste(string, gene_name, row_string, sep="")
        }
        return(string)
    }
)

NaviCell$methods(
    geneList2string = function(gene_list) {
    "Convert a gene list R object to a formatted string (internal utility)."
        string <- paste('@DATA\ngenes\n')
        gene_string <- paste(gene_list, sep="", collapse="\n")
        string <- paste(string, gene_string, sep="")
        string <- paste(string, "\n", sep="")
        return(string)
    }
)
 

#------------------------------------------------------------------------------
#
#  Drawing Configuration Dialog functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    drawingConfigOpen = function(...) {
    "Open drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('open'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigClose = function(...) {
    "Close drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('close'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigApply = function(...) {
    "Apply changes to drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigApplyAndClose = function(...) {
    "Apply changes to drawing configuration and close dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply_and_close'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigCancel = function(...) {
    "Cancel changes to drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('cancel'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectHeatmap = function(checked) {
    "Select heatmap display in drawing configuration dialog. checked = TRUE or FALSE."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_heatmap', checked)), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectBarplot = function(checked) {
    "Select barplot display in drawing configuration dialog. checked = TRUE or FALSE."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_barplot', checked)), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectGlyph = function(glyph_num, checked) {
    "Select glyph display in drawing configuration dialog. glyph_num = glyph number, checked = TRUE or FALSE."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_glyph', glyph_num, checked)), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectMapStaining = function(checked) {
    "Select map staining display in drawing configuration dialog. checked = TRUE or FALSE."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_map_staining', checked)), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectDisplayAllGenes= function(...) {
    "Select 'Display all genes' option in drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('display_all_genes'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    drawingConfigSelectDisplaySelectedGenes= function(...) {
    "Select 'Display selected genes' option in drawing configuration dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('display_selected_genes'), msg_id = .self$msg_id, action = 'nv_drawing_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

#------------------------------------------------------------------------------
#
#  MyData Dialog functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    mydataDialogOpen = function(...) {
    "Open MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('open'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogClose = function(...) {
    "Close MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('close'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogSetDatatables = function(...) {
    "Set Datatables tab active for MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('select_datatables'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogSetSamples = function(...) {
    "Set Samples tab active for MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('select_samples'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogSetGenes = function(...) {
    "Set Genes tab active for MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('select_genes'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogSetGroups = function(...) {
    "Set Groups tab active for MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('select_groups'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mydataDialogSetModules = function(...) {
    "Set Modules tab active for MyData Dialog."
        .self$incMessageId()
        list_param <- list(module='', args = array('select_modules'), msg_id = .self$msg_id, action = 'nv_mydata_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


#------------------------------------------------------------------------------
#
# Glyph Editor functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    glyphEditorOpen = function(glyph_num) {
    "Open the glyph editor. glyph_num = glyph number"
        .self$incMessageId()
        list_param <- list(module='', args = array(c('open', glyph_num)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorClose = function(glyph_num) {
    "Close the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('close', glyph_num)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorApply = function(glyph_num) {
    "Apply changes in the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply', glyph_num)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorApplyAndClose = function(glyph_num) {
    "Apply changes and close the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply_and_close', glyph_num)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorCancel = function(glyph_num) {
    "Cancel changes and close the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('cancel', glyph_num)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorSelectSample = function(glyph_num, sample_name) {
    "Select sample or group in the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_sample', glyph_num, sample_name)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorSelectShapeDatatable = function(glyph_num, datatable_name) {
    "Select datatable for glyph shape in the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable_shape', glyph_num, datatable_name)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorSelectColorDatatable = function(glyph_num, datatable_name) {
    "Select datatable for glyph color in the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable_color', glyph_num, datatable_name)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorSelectSizeDatatable = function(glyph_num, datatable_name) {
    "Select datatable for glyph size in the glyph editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable_size', glyph_num, datatable_name)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    glyphEditorSetTransparency = function(glyph_num, value) {
    "Set transparency parameter in the glyph editor. value = integer between 1 and 100."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_transparency', glyph_num, value)), msg_id = .self$msg_id, action = 'nv_glyph_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)



#------------------------------------------------------------------------------
#
# Barplot Editor functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    barplotEditorOpen = function(...) {
    "Open the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('open'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorClose = function(...) {
    "Close the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('close'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorApply = function(...) {
    "Apply changes for the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorApplyAndClose = function(...) {
    "Apply changes and close the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply_and_close'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorCancel = function(...) {
    "Cancel changes and close the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('cancel'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorSelectSample = function(col_num, sample_name) {
    "Select a sample or a group in the barplot editor. col_num = column index number, sample_name = sample or group name"
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_sample', col_num, sample_name)), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorSelectDatatable = function(datatable_name) {
    "Select a datatable in the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable', datatable_name)), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorClearSamples = function(...) {
    "Clear all samples in the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('clear_samples'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorSelectAllSamples = function(...) {
    "Select all samples in the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('all_samples'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorSelectAllGroups = function(...) {
    "Select all groups in the barplot editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('all_groups'), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    barplotEditorSetTransparency = function(value) {
    "Select transparency parameter in the barplot editor. value = integer between 1 and 100"
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_transparency', value)), msg_id = .self$msg_id, action = 'nv_barplot_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


#------------------------------------------------------------------------------
#
# Heatmap Editor functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    heatmapEditorOpen = function(...) {
    "Open the heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('open'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorClose = function(...) {
    "Close the heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('close'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorCancel = function(...) {
    "Cancel changes and close the heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('cancel'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorApply = function(...) {
    "Apply changes for the heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorApplyAndClose = function(...) {
    "Apply changes and close the heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('apply_and_close'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorSelectSample = function(col_num, sample_name) {
    "Select sample or group in heatmap editor. col_num = editor column number, sample_name = sample or group name."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_sample', col_num, sample_name)), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorSelectDatatable = function(row_num, datatable_name) {
    "Select datatable in heatmap editor. row_num = editor row number."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable', row_num, datatable_name)), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorClearSamples = function(...) {
    "Clear all samples in heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('clear_samples'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorSelectAllSamples = function(...) {
    "Select all samples in heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('all_samples'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorSelectAllGroups = function(...) {
    "Select all groups in heatmap editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('all_groups'), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    heatmapEditorSetTransparency = function(value) {
    "Set transparency parameter in heatmap editor. value = integer between 1 and 100."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_transparency', value)), msg_id = .self$msg_id, action = 'nv_heatmap_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

#------------------------------------------------------------------------------
#
# Unordered Discrete Configuration Editor functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    unorderedConfigOpen = function(datatable_name, datatable_parameter) {
    "Open unordered discrete configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('open', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigClose = function(datatable_name, datatable_parameter) {
    "Open unordered discrete configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('close', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigCancel= function(datatable_name, datatable_parameter) {
    "Cancel changes for unordered discrete configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('cancel', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


NaviCell$methods(
    unorderedConfigApply = function(datatable_name, datatable_parameter) {
    "Apply changes to unordered discrete configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigApplyAndClose = function(datatable_name, datatable_parameter) {
    "Apply changes to unordered discrete configuration editor for a given type of parameter, and close the window. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply_and_close', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetAdvancedConfig = function(datatable_name, datatable_parameter, checked) {
    "Open/close advanced configuration for unordered discrete configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_advanced_configuration', datatable_name, datatable_parameter, checked)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetDiscreteValue= function(datatable_name, datatable_parameter, sample_or_group, index, value) {
    "Set discrete value for unordered discrete configuration editor for a given
    type of parameter. datatable_parameter = 'shape' or
    'color' or 'size', sample_or_group = 'sample' or 'group', index = integer,
    value = double."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_discrete_value', datatable_name, datatable_parameter, sample_or_group, index, value)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetDiscreteColor= function(datatable_name, sample_or_group, index, color) {
    "Set color value for unordered discrete configuration editor. sample_or_group = 'sample' or 'group', index = integer, color = string hex code color value, e.g. 'FF0000'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_discrete_color', datatable_name, "color", sample_or_group, index, color)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetDiscreteSize= function(datatable_name, sample_or_group, index, size) {
    "Set size value for unordered discrete configuration editor. sample_or_group = 'sample' or 'group', index = integer, size = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_discrete_size', datatable_name, "size", sample_or_group, index, size)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetDiscreteShape = function(datatable_name, sample_or_group, index, shape) {
    "Set shape value for unordered discrete configuration editor. sample_or_group = 'sample' or 'group', index = integer, shape = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_discrete_shape', datatable_name, "shape", sample_or_group, index, shape)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSetDiscreteCondition = function(datatable_name, datatable_parameter, sample_or_group, index, condition) {
    "Set condition value for unordered discrete configuration editor. datatable_parameter = 'size' or 'shape' or 'color'. sample_or_group = 'sample' or 'group', index = integer, condition = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_discrete_cond', datatable_name, datatable_parameter, sample_or_group, index, condition)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSwitchSampleTab = function(datatable_name, datatable_parameter) {
    "Switch to sample tab for unordered discrete configuration editor. datatable_parameter = 'size' or 'shape' or 'color'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('switch_sample_tab', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    unorderedConfigSwitchGroupTab = function(datatable_name, datatable_parameter) {
    "Switch to group tab for unordered discrete configuration editor. datatable_parameter = 'size' or 'shape' or 'color'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('switch_group_tab', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_unordered_discrete_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


#------------------------------------------------------------------------------
#
# Continuous Configuration Editor functions 
#
#------------------------------------------------------------------------------

NaviCell$methods(
    continuousConfigOpen = function(datatable_name, datatable_parameter) {
    "Open continuous configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('open', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigClose = function(datatable_name, datatable_parameter) {
    "Close continuous configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('close', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigApplyAndClose = function(datatable_name, datatable_parameter) {
    "Apply changes and close continuous configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply_and_close', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigCancelAndClose = function(datatable_name, datatable_parameter) {
    "Cancel changes and close continuous configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('cancel', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetAbsVal = function(datatable_parameter, datatable_name, checked) {
    "Set absolute value mode for continuous configuration editor for a given type of parameter. datatable_parameter = 'shape' or 'color' or 'size', checked = TRUE or FALSE."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_sample_absval', datatable_parameter, datatable_name, checked)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetSampleMethod = function(datatable_parameter, datatable_name, method_index) {
    "Set the method used when multiple symbols map to the same entity. datatable_parameter = 'shape' or 'color' or 'size', method_index = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_sample_method', datatable_parameter, datatable_name, method_index)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetGroupMethod = function(datatable_parameter, datatable_name, method_index) {
    "Set the method used when multiple symbols map to the same entity. datatable_parameter = 'shape' or 'color' or 'size', method_index = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_group_method', datatable_parameter, datatable_name, method_index)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetSelectionSize = function(datatable_name, sample_or_group, index, size) {
    "Set the size selection to a given value for the 'size' parameter. sample_or_group = 'sample' or 'group', index = integer, size = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_select_size', datatable_name, "size", sample_or_group, index, size)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


NaviCell$methods(
    continuousConfigSetSelectionShape = function(datatable_name, sample_or_group, index, shape) {
    "Set the shape selection to a given value for the 'shape' parameter. sample_or_group = 'sample' or 'group', index = integer, shape = integer."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_select_shape', datatable_name, "shape", sample_or_group, index, shape)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


NaviCell$methods(
    continuousConfigSwitchSampleTab = function(datatable_name, datatable_parameter) {
    "Switch continuous configuration editor window to 'sample' tab. Parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('switch_sample_tab', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSwitchGroupTab = function(datatable_name, datatable_parameter) {
    "Switch continuous configuration editor window to 'group' tab. Parameter = 'shape' or 'color' or 'size'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('switch_group_tab', datatable_name, datatable_parameter)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


NaviCell$methods(
    continuousConfigSetStepCount = function(sample_or_group, datatable_parameter, datatable_name,  step_count) {
    "Set continuous configuration step count parameter to a given value. sample_or_group = 'sample' or 'group'. parameter = 'shape' or 'color' or 'size' step_count = integer value."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('step_count_change', sample_or_group, datatable_parameter, datatable_name, step_count)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetColorAt = function(datatable_name, sample_or_group, index, color_hex_value) {
    "Set continuous configuration color value. sample_or_group = 'sample' or 'group'."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_input_color', datatable_name, 'color', sample_or_group, index, color_hex_value)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigSetValueAt = function(datatable_name, parameter_type, sample_or_group, index, continuous_value) {
    "Set continuous configuration continuous value at a given index. sample_or_group = 'sample' or 'group'. parameter_type = 'size' or 'shape' or 'color'. "
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_input_value', datatable_name, parameter_type, sample_or_group, index, continuous_value)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    continuousConfigApply = function(datatable_name, parameter_type) {
    "Apply changes to the continuous configuration editor. parameter_type = 'size' or 'shape' or 'color'. "
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply', datatable_name, parameter_type)), msg_id = .self$msg_id, action = 'nv_display_continuous_config_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


#------------------------------------------------------------------------------
#
# Map Staining  functions 
#
#------------------------------------------------------------------------------


NaviCell$methods(
    mapStainingEditorSelectDatatable = function(datatable_name) {
    "Select a datatable for the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_datatable', datatable_name)), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorSelectSample = function(sample_name) {
    "Select a sample for the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('select_sample', sample_name)), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorApply = function(...) {
    "Apply modifications for the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply')), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorCancel = function(...) {
    "Cancel modifications and close the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('cancel')), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorOpen = function(...) {
    "Open the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array('open'), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorClose = function(...) {
    "Close the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('close')), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorApplyAndClose = function(...) {
    "Apply changes and close the map staining editor."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('apply_and_close')), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    mapStainingEditorSetTransparency = function(transparency_value) {
    "Set the transparency value parameter for the map staining editor (integer value between 1 and 100)."
        .self$incMessageId()
        list_param <- list(module='', args = array(c('set_transparency'), transparency_value), msg_id = .self$msg_id, action = 'nv_map_staining_editor_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F))
        #print(.self$formatResponse(response))
    }
)


#------------------------------------------------------------------------------
#
# Sample Annotation functions 
#
#------------------------------------------------------------------------------


NaviCell$methods(
    importSampleAnnotationFromFile = function(fileName) {
        data_string <- .self$file2dataString(fileName)
        if (!is.null(data_string)) {
            if (nchar(data_string) < .self$packsize) {
                .self$incMessageId()
                list_param <- list(module='', args = list("import", data_string), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
                str_data <- .self$makeData(.self$formatJson(list_param))
                .self$incMessageId()
                response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
                #print(.self$formatResponse(response))
            }
            else {
                list_param <- list(module='', args = list("import", data_string), action = 'nv_sample_annotation_perform')
                fill_cmd <- .self$makeData(.self$formatJson(list_param))
                .self$sendBigData(fill_cmd)
                .self$waitForImported()
            }
        }
    }
)

NaviCell$methods(
    sampleAnnotationSelectAnnotation = function(annotation_name, true_or_false) {
    "Select or un-select an annotation for creating groups from a sample annotation table. true_or_false = TRUE, select, true_or_false = FALSE, un-select."
        .self$incMessageId()
        list_param <- list(module='', args = list("select_annotation", annotation_name, true_or_false), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))

        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    sampleAnnotationApply = function(...) {
    "Apply the modifications done on a sample annotation table."
        .self$incMessageId()
        list_param <- list(module='', args = list("apply"), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))

        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
        #print(.self$formatResponse(response))

    }
)

NaviCell$methods(
    sampleAnnotationOpen = function(...) {
    "Open sample annotation dialog."
        .self$incMessageId()
        list_param <- list(module='', args = list("open"), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))

        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
        #print(.self$formatResponse(response))

    }
)

NaviCell$methods(
    sampleAnnotationClose = function(...) {
    "Close sample annotation dialog."
        .self$incMessageId()
        list_param <- list(module='', args = list("close"), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))

        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
        #print(.self$formatResponse(response))
    }
)

NaviCell$methods(
    sampleAnnotationCancel = function(...) {
    "Cancel changes and close sample annotation dialog."
        .self$incMessageId()
        list_param <- list(module='', args = list("cancel"), msg_id = .self$msg_id, action = 'nv_sample_annotation_perform')
        str_data <- .self$makeData(.self$formatJson(list_param))

        .self$incMessageId()
        response <- postForm(.self$proxy_url, style = 'POST', id = .self$session_id, msg_id = .self$msg_id, mode='cli2srv', perform='send_and_rcv', data=str_data, .opts=curlOptions(ssl.verifypeer=F)) 
        #print(.self$formatResponse(response))
    }
)

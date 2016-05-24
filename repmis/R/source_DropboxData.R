#' No longer supported
#'
#' @export

source_DropboxData <-function()
{
    stop('Unfortunately, source_DropboxData is no longer supported due to changes in the Dropbox API.', call. = FALSE)
#    url <- paste0('https://dl.dropboxusercontent.com/s/',
#                    key, '/', file)

#    temp_file <- tempfile()
#    on.exit(unlink(temp_file))

#  key <- NULL
#    key <- as.list(url)
#    if (isTRUE(clearCache)){
#        Found <- findCache(key = key)
#        if (is.null(Found)){
#            message('Data not in cache. Nothing to remove.')
#        }
#        else if (!is.null(Found)){
#            message('Clearing data from cache.')
#            file.remove(Found)
#        }
#    }

#    if (isTRUE(cache)){
#        data <- loadCache(key)
#        if (!is.null(data)){
#            message('Loading cached data.\n')
#            return(data);
#        }
#        data <- download_data_intern(url = url, sha1 = sha1,
#                                    temp_file = temp_file)
#        data <- fread(data, sep = sep, header = header, data.table = F, ...)
#        saveCache(data, key = key)
#        data;
#    }
#    else if (!isTRUE(cache)){
#        data <- download_data_intern(url = url, sha1 = sha1,
#                                    temp_file = temp_file)
#        data <- fread(data, sep = sep, header = header, data.table = F, ...)
#        return(data)
#    }
}

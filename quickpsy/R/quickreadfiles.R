#' Reads several files
#'
#' \code{quickreadfiles} builts a data frame from several txt files. It
#' assumes that in each file, the first row has the names of the variables.
#' @param path Path of the file (default is the working directory).
#' @param ... arguments of the form name_var = c('value1', 'value2',..).
#' A new column with variable name name_var is addes to the data frame.
#' @examples
#' # download the 3 files in
#' # https://github.com/danilinares/quickpsy/tree/master/inst/extdata/example1
#' # and add them to your working directory
#' # dat <- quickreadfiles(subject = c('aa', 'bb', 'cc'), session = c('1', '2'))
#' # fit <- quickpsy(dat, phase, resp, grouping=.(subject), lapses = T, guess = T)
#' # plotcurves(fit)
#' @import dplyr
#' @export

quickreadfiles <- function(path = getwd(), ...) {

  arguments <- c(as.list(environment()), list(...))
  arguments[1] <- NULL
  exist<-NULL # Joan added 1-4-2015

  namesfun <- function(d) {
    namefile <- paste0(path,'/', paste(unlist(d), collapse = ''), '.txt')
    data.frame(namefile, exist = file.exists(namefile), stringsAsFactors=F)
  }


  namefiles <- expand.grid(arguments) %>%
    group_by_(.dots = names(arguments)) %>%
    do(namesfun(.))

  namefiles %>% filter(exist) %>%
    group_by_(.dots = names(arguments)) %>%
    do(read.table(.$namefile, header = T))
}

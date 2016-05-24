##' This function checks whether there are overlapping between the
##' transitional countries.
##'
##'
##' @param old The FAOST_CODE of the old countries
##' @param new The FAOST_CODE of the new countries
##' @param var The variable to be checked
##' @param year The column which index the time.
##' @param data The data frame
##' @param take The type of check/replacement to be done.
##'
##' @export
##'


overlap = function(old, new, var, year = "Year", data, take){

    n.var = length(var)
    for(i in 1:n.var){
        old_years = unique(data[data$FAOST_CODE %in% old &
                           !is.na(data[, var[i]]), year])
        new_years = unique(data[data$FAOST_CODE %in% new &
                           !is.na(data[, var[i]]), year])
        int_years = intersect(old_years, new_years)
        if(length(int_years) > 0)
            cat(paste0(var[i], " has overlap in:\n\t Year = ",
                       paste0(int_years, collapse = ", "),
                       "\n\tfor FAOST_CODE = ",
                       paste0(c(old, new), collapse = ", "), "\n"))
        switch(take,
               simpleCheck = {},
               takeNew = {
                 data[data$FAOST_CODE %in% old & 
                        data[, year] %in% int_years, var[i]] = NA
               },
               takeOld = {
                 data[data$FAOST_CODE %in% new & 
                        data[, year] %in% int_years, var[i]] = NA
               },
               complete = {
                 for(j in sort(unique(int_years))){
                   if(NROW(data[data[, "FAOST_CODE"] %in% new & 
                                  !is.na(data[, var[i]]) & data[, year] == j, "FAOST_CODE"]) == length(new)) {
                     ## All the new entities are there. Remove the old.
                     data[data[, "FAOST_CODE"] %in% old & data[, year] == j, var[i]] = NA
                   } else {
                     ## Not all the new entities are there. Remove the new.
                     data[data[, "FAOST_CODE"] %in% new & data[, year] == j, var[i]] = NA
                   }
                 }
               }
        )
    }
    data
}

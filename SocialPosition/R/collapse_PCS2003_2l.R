collapse_PCS2003_2l <-
function(PCS2003_2l, data){
#collapsing in 1 levels
data$PCS2003_1l[PCS2003_2l %in% c(10)] <- 1
data$PCS2003_1l[PCS2003_2l %in% c(21,22,23)] <- 2
data$PCS2003_1l[PCS2003_2l %in% c(31,32,36)] <- 3
data$PCS2003_1l[PCS2003_2l %in% c(41,46,47,48)] <- 4
data$PCS2003_1l[PCS2003_2l %in% c(51,54,55,56)] <- 5
data$PCS2003_1l[PCS2003_2l %in% c(61,66,69)] <- 6
return(data)
}

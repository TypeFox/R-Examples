collapse_PCS2003_3l <-
function(PCS2003_3l, data){
#collapsing in 2 levels
data$PCS2003_2l[PCS2003_3l %in% c(11,12,13)] <- 10
data$PCS2003_2l[PCS2003_3l %in% c(21)] <- 21
data$PCS2003_2l[PCS2003_3l %in% c(22)] <- 22
data$PCS2003_2l[PCS2003_3l %in% c(23)] <- 23
data$PCS2003_2l[PCS2003_3l %in% c(31)] <- 31
data$PCS2003_2l[PCS2003_3l %in% c(33,34,35)] <- 32
data$PCS2003_2l[PCS2003_3l %in% c(37,38)] <- 36
data$PCS2003_2l[PCS2003_3l %in% c(42,43,44,45)] <- 41
data$PCS2003_2l[PCS2003_3l %in% c(46)] <- 46
data$PCS2003_2l[PCS2003_3l %in% c(47)] <- 47
data$PCS2003_2l[PCS2003_3l %in% c(48)] <- 48
data$PCS2003_2l[PCS2003_3l %in% c(52,53)] <- 51
data$PCS2003_2l[PCS2003_3l %in% c(54)] <- 54
data$PCS2003_2l[PCS2003_3l %in% c(55)] <- 55
data$PCS2003_2l[PCS2003_3l %in% c(56)] <- 56
data$PCS2003_2l[PCS2003_3l %in% c(62,63,64,65)] <- 61
data$PCS2003_2l[PCS2003_3l %in% c(67,68)] <- 66
data$PCS2003_2l[PCS2003_3l %in% c(69)] <- 69
#collapsing in 1 levels
data$PCS2003_1l[data$PCS2003_2l %in% c(10)] <- 1
data$PCS2003_1l[data$PCS2003_2l %in% c(21,22,23)] <- 2
data$PCS2003_1l[data$PCS2003_2l %in% c(31,32,36)] <- 3
data$PCS2003_1l[data$PCS2003_2l %in% c(41,46,47,48)] <- 4
data$PCS2003_1l[data$PCS2003_2l %in% c(51,54,55,56)] <- 5
data$PCS2003_1l[data$PCS2003_2l %in% c(61,66,69)] <- 6
return(data)
}

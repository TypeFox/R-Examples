check.table <- function (ecopath) 
{
if (is.data.frame(ecopath) == FALSE) 
ecopath <- as.data.frame(ecopath)
names <- names(ecopath)
wanted <- c("group_name", "TL", "biomass", "prod", "accessibility", "OI")
verif <- is.na(match(wanted, names))
pb <- paste(wanted[which(verif == TRUE)], collapse = " ")
if (pb != "") 
cat(paste("The column(s) ", pb, " is(are) not present (if the column OI is not present, it's not an issue but the use of the OI smooth (smooth_type=3) won't be possible).\n"))
if (length(grep("catch", names)) == 0) 
print("No fleet/catches are detected. Even if no catches are made (MPA or no data), a column with 0 value must be entered. The name of the column has to be written as 'catch.something'.")

cat("\nCHECK OF THE NA VALUES, IF A MESSAGE APPEARS, PLEASE READ IT CAREFULLY AND FOLLOW THE INSTRUCTIONS. NO NA IS ACCEPTED AS INPUT DATA.\n\n")
for (i in wanted[which(verif == FALSE)]) {
if (match(NA, ecopath[, i], FALSE) > 0) {
if (i == "prod") {
pb <- paste(ecopath[is.na(ecopath$prod), ]$group_name, collapse = " ")
cat(paste("The column prod contains a NA for the group", pb, ", that's an issue if the concerned groups are not detritus. If the concerned group is a detritus group, change the NA value by 0; otherwhise check your table and correct by the right value (use fix() for example or reload it after correction).\n"))
}
else {
toto <- paste(i)
cat(paste("The column", i, "contains a NA, that's a problem! Check and correct the table (you can use fix()).\n"))
}}
}

for (pecheries in colnames(ecopath)[grep("catch", colnames(ecopath))]) {
if (match(NA, ecopath[, pecheries], FALSE) > 0) {
cat(paste("The column", pecheries, "contains NA. That's a problem! Even if no catches are made (MPA or no data), 0 value must be entered. Use fix() to change the NA by the proper value or reload the dataset after correction.\n"))
}}

cat("\nCHECK OF THE ACCESSIBILITY VALUES, IF A MESSAGE APPEARS, PLEASE READ IT CAREFULLY AND FOLLOW THE INSTRUCTIONS.\n\n")
for (pecheries in colnames(ecopath)[grep("catch", colnames(ecopath))]) {
for (i in 1:length(ecopath$group_name)) {
if ((ecopath[i, pecheries] > 0) && (ecopath[i,"accessibility"] == 0)){
cat(paste("The group", ecopath[i,"group_name"], "has a null accessibility, that's a problem! He is fished by the fishery", pecheries, "and consequently needs a positive value.\n"))
}}}

for (i in 1:length(ecopath$group_name)) {
if ((ecopath[i,"accessibility"] > 1)){
cat(paste("The group", ecopath[i,"group_name"], "has an accessibility value higher to 1, that's a problem! Choose a value between 0 and 1.\n"))
}}
}
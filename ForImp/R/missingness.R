missingness <-
function(mat)
{
mu<-table(rowSums(is.na(mat)))
#row.names(mu)<-"#units"
mv<-colSums(is.na(mat))
list(number_of_missing_values=sum(mv),missing_values_per_unit=mu,missing_values_per_variable=mv)
}


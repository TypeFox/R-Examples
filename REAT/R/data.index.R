data.index <-
function (dataset, col_index, col_ref, value_ref) {

refcolumn <- dataset[[col_ref]]
 
refvaluerows <- which (dataset[[col_ref]] == value_ref)

index_values <- dataset[[col_index]]/dataset[[col_index]][refvaluerows]*100
return(index_values)
}

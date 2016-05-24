LAI_function <-
function(Number_of_leaves_by_plant,Leaf_Area,row_distance,in_row_distance){

 LAI <- Number_of_leaves_by_plant * Leaf_Area / (row_distance * in_row_distance)

return(LAI)

}


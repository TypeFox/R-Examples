`array.id` <-
function(x){
## x = slidedescription table from read.slidedescription

          ## paste the arraydescribing columns to one vector
            arrays <- paste (x[,"pad"],x[,"slide"],x[,"incubation_run"],x[,"spotting_run"],sep="-")
            print (paste("identified",length(arrays),"arrays"))
          ## return character vector with array identifier
            return(arrays)
}


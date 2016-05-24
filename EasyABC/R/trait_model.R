## Toy model
trait_model <- function(input = c(1, 1, 1, 1, 1, 1)) {
    .C("trait_model", input = c(input[1], 500, input[2:3], 1, input[4:5]), stat_to_return = array(0, 
        4))$stat_to_return
} 



setMethod ('print', 'EotMode', 
           function(x, ...) {
             if (inherits(x, 'EotMode')) {
               show(x)
             }
           }
)


setMethod ('print', 'EotStack', 
           function(x, ...) {
             if (inherits(x, 'EotStack')) {
               show(x)
             }
           }
)
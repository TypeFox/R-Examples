setClass("MyContControl",
    # definition of additional properties
    contains = "VirtualContControl")

setMethod("contaminate",
    signature(x = "data.frame", control = "MyContControl"),
    function(x, control, i) {
        # method definition
    })
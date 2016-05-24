setClass("MyNAControl",
    # definition of additional properties
    contains = "VirtualNAControl")

setMethod("setNA",
    signature(x = "data.frame", control = "MyNAControl"),
    function(x, control, i) {
        # method definition
    })
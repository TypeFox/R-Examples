setClass("MyDataControl",
    # class definition
    contains = "VirtualDataControl")

setMethod("generate",
    signature(control = "MyDataControl"),
    function(control) {
        # method definition
    })
getFonts <- function(font){

    ## helper function used by table{Continuous, Nominal}
    fontify <- if(font == ""){
        function(label, prefix){label}} else {
        function(label, prefix){paste("\\", prefix, font, "{", label, "}", sep = "")}
        }

    return(list(math =
                ## function to format math in the header line 
                function(label)math(fontify(label, prefix="math")),
                text=
                ## function to format text in the header line 
                function(label) fontify(label, prefix="text")))    
}

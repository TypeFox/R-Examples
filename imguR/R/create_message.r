create_message <-
function(recipient,
         body,
         ...){
    out <- imgurPOST('conversations/', 
                     body = list(recipient = recipient,
                                 body = body),
                     ...)
    structure(out, class = 'imgur_basic')
}

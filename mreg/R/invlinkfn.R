"invlinkfn" <-
function(x, family, ...){
args <- list(...)
switch( family,
"identity"=x,
"log"=log(x),
"logit"=log(x/(1-x))       ,
"hyper"=1/x-1 
)       

}


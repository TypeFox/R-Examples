"linkfn" <-
function(x, family, ...){
args <- list(...)
switch( family,
"identity"=x,
"log"=exp(x),
"logit"=exp(x)/(1+exp(x))       ,
"hyper"=1/(1+x)   
)       

}


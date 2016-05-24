myfn<-function(prm, y = NULL, t = NULL){;
for (i in 1:length(prm) ){
joe<-paste(names(prm)[[i]],"<-",prm[[i]]);
 eval(parse(text=joe));
};
 eval(crossprod(b1/(1 + b2 * exp(-1 * b3 * t)) - y))
 }

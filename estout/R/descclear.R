`descclear` <-
function(store="default"){
    prev.list <- paste(store,".dcl",sep="")
assign(prev.list,NULL,estout:::estoutstorage)
}

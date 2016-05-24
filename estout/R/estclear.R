`estclear` <-
function(store="default"){
    prev.list <- paste(store,".ccl",sep="")
assign(prev.list,NULL,estout:::estoutstorage)
}

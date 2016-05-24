add.col <-
function(df, new.col, MaxTime) {n.row<-dim(df)[1]
	if(is.null(n.row)){n.row<-MaxTime}
           length(new.col)<-n.row
           cbind(df, new.col)
    }


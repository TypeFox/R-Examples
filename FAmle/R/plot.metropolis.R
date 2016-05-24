plot.metropolis <-
function(x,plot.type='carlin',pos=1:x$iter,...)
	switch(plot.type,ts=plot.ts(x$sims,...),pairs=pairs(x$sims[pos,],...),hist=hist(x,...),
	post.pred=Plot.post.pred(x,...),carlin=Carlin(x))


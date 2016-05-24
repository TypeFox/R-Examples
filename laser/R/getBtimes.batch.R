`getBtimes.batch` <-
function(fname, basal = NULL)
{
	
	btfx <- function(x){
		return(rev(sort(as.numeric(branching.times(x)))));
	}
	
	scaleFx <- function(x, z){
		sf <- z/max(x);
		return(x*sf);
	}
	v <- read.tree(fname);
	if (length(v) ==1)
		stop('Error\nThis function cannot be used on files containing a single tree\n');
	bt <- lapply(v, branching.times);
	if(!is.null(basal))
		bt <- lapply(bt, scaleFx, z=basal);
	bt <- as.matrix(t(as.data.frame(bt)));
	return(bt);
}


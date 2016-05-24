#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
follow_points <- function(series, m, d, ref, k, s, nearest) {
	res <- numeric(s)
	nearest[is.na(nearest)] <- -1
	.C("follow_points",series=as.double(series), m=as.integer(m), d=as.integer(d), length=as.integer(length(series)), nref=as.integer(length(ref)), nrow=as.integer(nrow(nearest)), k=as.integer(k), s=as.integer(s), nearest=as.integer(nearest), ref=as.integer(ref), res=as.double(res),PACKAGE="tseriesChaos")[["res"]]
}

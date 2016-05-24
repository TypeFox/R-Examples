PlusMinus <- function(x, y, head, xhead, yhead, digits=2, ...) {
    head <- if (missing(head)) quote(Heading())
	else substitute(Heading(head), list(head=substitute(head)))
    xhead <- if (missing(xhead)) quote(Heading())
	else substitute(Heading(head), list(head=substitute(xhead)))
    yhead <- if (missing(yhead)) quote(Heading())
	else substitute(Heading(head), list(head=substitute(yhead)))
    fmt <- function(x){
		 s <- format(x, digits=digits, ...)
		 is_stderr <- (1:length(s)) > length(s) %/% 2
		 s[is_stderr] <- sprintf("$%s$", s[is_stderr])
		 s[!is_stderr] <- latexNumeric(s[!is_stderr])
		 s
	       }
    substitute( head*1*Format(fmt())*(xhead*Justify("r@{}")*x 
                                    + yhead*Justify("@{ $\\pm$ }l")*y),
               list(x=substitute(x), y=substitute(y), 
               head=head, xhead=xhead, yhead=yhead, fmt=fmt))
}

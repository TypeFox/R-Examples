# Evaluates the significance of rank correlations for a pair of rankings
rankor<-function(p, q, index = "spearman", approx = "exact", tiex = "woodbury", CC = FALSE, type = "two-sided", print = TRUE, sizer = 100000, repgin = 1000){
	 if (!is.numeric(p)) {stop("Non-numeric argument to mathematical function")}
	 if (!is.matrix(p) & !is.numeric(q)) {stop("Non-numeric argument to mathematical function")}
	 nx<-length(p);ny<-length(q);n<-max(nx,ny)
	 if (!(nx==ny)) {stop("not all arguments have the same length")}
	 ain<-c("spearman","kendall","gini","r4","fy","filliben")
	 apx<-c("gaussian","student","vggfr","exact")
	 alter<- c("two-sided", "greater","less")
	 tos<-c("woodbury","gh","wgh","midrank","dubois","untied rankings")
	 index<-tolower(index)
	 index<-match.arg(index, ain, several.ok = TRUE) 
	 approx<-tolower(approx);approx<-match.arg(approx, apx, several.ok = TRUE)
	 tiex<-tolower(tiex);tiex<-match.arg(tiex, tos, several.ok = TRUE)
	 type=tolower(type);type<-match.arg(type, alter, several.ok = TRUE)
   	 ra<-comprank(p,q,index,tiex,sizer,repgin,print)
   	 if(print) {cat("Method for breaking ties:",tos[ra$ities],"\n")}
   	 if(ra$ities<6 & print) {print("p-values of rank correlation tests are only indicative in presence of ties.")}
   	 if(ra$ities<6 & print ) {print("The continuity correction should not be considered with ties present.")}
   	 a<-ranktes(ra$r, n, index, approx, CC, type, print)
	return(a)}
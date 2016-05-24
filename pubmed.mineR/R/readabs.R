readabs <-
function(x)
{
xa <- readLines(con = x); 
xa <- c("PMID: MOCK","",xa);
xb <- regexpr("PMID", xa); 

xc = unlist(lapply(1:length(xb), function(x){ if (xb[x] == 1) return(x)}))

jn <- NULL; abs <- NULL; pmid <- NULL; 

xd = lapply(1:length(xc), function(x){k = x + 1; if(k <= length(xc)) {tempjn <- xc[x] + 3; jn = xa[tempjn]; tempabs1 <- tempjn + 1; 

tempabs2 <- xc[k]-1; abs = paste(xa[tempabs1:tempabs2], collapse = " "); temppmid <- xc[k];temppmid2 <- unlist(strsplit(xa[temppmid], " ", 

fixed = TRUE)); pmid = temppmid2[2]; return(c(jn,abs,pmid))}})

jn = unlist(lapply(xd,function(x){return(x[1])})); abs = unlist(lapply(xd,function(x){return(x[2])})); pmid = 

unlist(lapply(xd,function(x){return(x[3])}));


write.table(cbind(jn,abs,pmid), file="newabs.txt", quote=FALSE,sep = "\t", row.names = FALSE); resultabs <- new("Abstracts",Journal = jn, 

Abstract = abs, PMID = as.numeric(pmid)); return(resultabs)};


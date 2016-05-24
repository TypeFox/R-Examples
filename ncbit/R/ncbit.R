print.taxdump <-
function(x, ...){
	if(nrow(x)>6){
		cat(paste("\nNCBI GenBank taxonomy assembled ", attributes(x)$date, "\n\n ...showing the first several entries...\n", sep=""))
	} else {
		cat(paste("\nNCBI GenBank taxonomy assembled ", attributes(x)$date, "\n\n", sep=""))
	}
	print(attributes(x)$header)	
}

ncbit <-
function(update=FALSE, ...){
	
	gb_path=.path_to_gb.dmp()
	rda=as.list(gb_path)$ncbi
	
	build=update
	if(!file.exists(rda) | build){
		build=TRUE
	} else {
		build=FALSE
		if(exists("ncbi")) return(ncbi)
	}
	
	if(build){
		if(file.exists("taxdump.tar.gz")) unlink("taxdump.tar.gz")
		cat("Please be patient as 'taxdump' is built from NCBI:\n\t*** download may take several minutes ***\n")
		if(!system("which curl", ignore.stdout=TRUE)==0) stop("Install 'curl' before proceeding.")
		if(!system("which gunzip", ignore.stdout=TRUE)==0) stop("Install 'gunzip' before proceeding.")
		if(!system("which tar", ignore.stdout=TRUE)==0) stop("Install 'tar' before proceeding.")
		if(!system("which perl", ignore.stdout=TRUE)==0) stop("Install 'perl' before proceeding.")
        
		system("curl -OL ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
		system("gunzip taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
		system("tar -xvf taxdump.tar", ignore.stderr=TRUE, wait=TRUE)
		system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' names.dmp")
		system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' nodes.dmp")
		
        FUN=function(rm=TRUE){
            mv=TRUE
            if(mv){
                if(!system("which mv", ignore.stdout=TRUE)==0) stop("Install 'mv' before proceeding.")
                system("rm -rf /tmp/idx")
                if(all(sapply(ff<-c("nodes.dmp", "names.dmp"), file.exists))){
                    system(paste("mv nodes.dmp", gb_path[["nodes"]], sep=" "))
                    system(paste("mv names.dmp", gb_path[["names"]], sep=" "))
                    cleanup=c("taxdump.tar", "readme.txt", "gc.prt", "merged.dmp", "division.dmp", "delnodes.dmp", "citations.dmp", "gencode.dmp")
                    cleanup=cleanup[cleanup%in%dir()]
                    if(rm) {
                        if(!system("which rm", ignore.stdout=TRUE)==0) stop("Install 'rm' before proceeding.")
                        system(paste("rm -f", paste(cleanup, collapse=" "), sep=" "))
                    }
                } else {
                    stop("Error encountered from NCBI: 'nodes.dmp' and (or) 'names.dmp' cannot be located.")
                }
                
                ## read in NCBI data
                names.dmp<-read.table(gb_path[["names"]],header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE)
                names.dmp<-names.dmp[,1:4]
                names(names.dmp)<-c("id", "node", "unique", "type")
                nodes.dmp<-read.table(gb_path[["nodes"]],header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE)
                nodes.dmp<-nodes.dmp[,c(1:5,11:13)]
                names(nodes.dmp)<-c("id","parent_id","rank","embl_code","division_id","GenBank_hidden_flag","hidden_subtree_root_flag","comments")
                ncbi=cbind(names.dmp, nodes.dmp[match(names.dmp$id, nodes.dmp$id), c("parent_id", "rank")])
                rownames(ncbi)=NULL
                ncbi[which(ncbi[,"node"]=="root"),"rank"]="root"
                attr(ncbi, "date")=Sys.Date()
                class(ncbi)=c("taxdump", class(ncbi))
                attr(ncbi, "header")=as.data.frame(head(ncbi))
                save(ncbi, file=rda)
                unlink(gb_path[["nodes"]])
                unlink(gb_path[["names"]])
            } else {
                warning("Attempting to revert to previously installed version", immediate.=TRUE)
            }
            FUN(...)
        }
	}
    
	ncbi=get(load(rda))
	return(ncbi)
}

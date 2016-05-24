setGeneTable <-function(geneInfo,table='gene',db='snplistdb') {
    if ( !all(c('gene','chr','start','end') %in% names(geneInfo)) ) {
        return(print("The column names should have 'gene', 'chr', 'start', 'end'."))
    }
    
    drv   <- dbDriver("SQLite")
    dbFile<- sprintf("%s.sqlite",db)
    conn  <- dbConnect(drv, dbname = dbFile)
    
    tables<-dbGetQuery(conn,"SELECT name FROM sqlite_master WHERE type='table';")
    if(nrow(tables)!=0 && any(table==tables)) {
        cat("remove the existing table :",table, "\n\n")
        dbGetQuery(conn,sprintf("DROP TABLE %s;",table))
    }  

    dbGetQuery(conn,sprintf("CREATE TABLE %s(gene,chr,start,end, primary key(gene))",table))
    dbWriteTable(conn,table,geneInfo,row.names=FALSE,append=TRUE)

    cat("Create Table :",table,"\n")
    print( dbGetQuery(conn,sprintf("SELECT * FROM %s LIMIT 10;",table)) )
    cat(".....\n\n")

    r<-dbGetQuery(conn,sprintf("SELECT COUNT(*) FROM %s;",table))
    dbDisconnect(conn)
    
    return(r)
}


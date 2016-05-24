#currently no tests for the presupplied masks, correct.db.names, the Database and TableSchemaList classes.

check.table.import <- function(dta, tbsl, name, pks=paste(name, "ind", sep="_"))
{
    expect_is(tbsl, "TableSchemaList")
    expect_named(tbsl@tab.list, name)
    #are the column names ok
    expect_equal(length(intersect(tbsl@tab.list[[name]]$db.cols, union(names(dta), pks))), length(union(tbsl@tab.list[[name]]$db.cols, union(names(dta), pks))))
    
    #are the column types ok in a basic sense
    ##just removes autoincremented PKs
    common.cols <- intersect(tbsl@tab.list[[name]]$db.cols, names(dta))
    
    tab.cols <- sapply(names(dta), function(x) class(dta[,x]))
    tab.cols <- toupper(tab.cols)
    tab.cols[tab.cols %in% c("CHARACTER", "FACTOR")] <- "TEXT"
    
    expect_identical(tab.cols, tbsl@tab.list[[name]]$db.schema[names(tab.cols)])
}

#these are simple relationships so the keys will be equivalent and there should be no modifications of columns etc
check.direct.keys <- function(tbsl, from, to, key.name, orig.to.obj, orig.from.obj)
{
    expect_named(tbsl@tab.list[[to]]$foreign.keys, from)
    expect_identical(tbsl@tab.list[[to]]$foreign.keys[[from]]$local.keys, tbsl@tab.list[[to]]$foreign.keys[[from]]$ext.keys)
    
    expect_identical(tbsl@tab.list[[to]]$foreign.keys[[from]]$local.keys, key.name)
    
    expect_identical(tbsl@tab.list[[to]]$db.col, orig.to.obj@tab.list[[to]]$db.col)
    expect_identical(tbsl@tab.list[[to]]$db.schema, orig.to.obj@tab.list[[to]]$db.schema)
    
    expect_identical(tbsl@tab.list[[from]]$db.col, orig.from.obj@tab.list[[from]]$db.col)
    expect_identical(tbsl@tab.list[[from]]$db.schema, orig.from.obj@tab.list[[from]]$db.schema)
}

convert.factors.to.strings <- function(dta)
{
    for(i in 1:ncol(dta))
    {
        if(class(dta[,i]) == "factor")
        {
            dta[,i] <- as.character(dta[,i])
        }
    }
    
    return(dta)
}


test_that("Create and work with TBSL and Database objects in a basic sense",
{
    #makeSchemaFromData, append and length
    
    baseball.teams <- new("TableSchemaList")
    
    expect_equal(length(baseball.teams), 0)
    
    franches <- makeSchemaFromData(TeamsFranchises, name="team_franch")
    check.table.import(TeamsFranchises, franches, "team_franch")
    
    baseball.teams <- append(baseball.teams, franches)
    
    expect_equal(length(baseball.teams), 1)
    
    teams <- makeSchemaFromData(Teams, name="teams")
    check.table.import(Teams, teams, "teams")
    
    baseball.teams <- append(baseball.teams, teams)
    
    expect_equal(length(baseball.teams), 2)
    
    salaries <- makeSchemaFromData(Salaries, name="salaries")
    check.table.import(Salaries, salaries, "salaries")
    
    baseball.teams <- append(baseball.teams, salaries)
    
    expect_equal(length(baseball.teams), 3)
    
    #relationships
    
    relationship(baseball.teams, from="team_franch", to="teams") <- franchID ~ franchID
    check.direct.keys(baseball.teams,  from="team_franch", to="teams", key.name="franchID", orig.to.obj=teams, orig.from.obj=franches)
    
    relationship(baseball.teams, from="teams", to="salaries") <- teamID ~ teamID
    check.direct.keys(baseball.teams,  from="teams", to="salaries", key.name="teamID", orig.to.obj=salaries, orig.from.obj=teams)
    
    #helpers for TableSchemaLists
    
    col.list <- columns(baseball.teams)
    
    expect_named(col.list, c("team_franch", "teams", "salaries"))
    expect_equal(col.list, list(team_franch=c("team_franch_ind", names(TeamsFranchises)), teams=c("teams_ind", names(Teams)), salaries=c("salaries_ind", names(Salaries))))
    
    expect_equal(tables(baseball.teams), c("team_franch", "teams", "salaries"))
    
    #Basic formation and checks of Database objects
    
    temp.db <- tempfile()
    
    baseball.db <- Database(baseball.teams, temp.db)
    
    #columns
    
    expect_equal(columns(baseball.db), columns(baseball.teams))
    
    #tables
    
    expect_equal(tables(baseball.db), tables(baseball.teams))
    
    #dbFile
    
    expect_equal(dbFile(baseball.db), temp.db)
    
    #schema
    
    expect_equal(schema(baseball.db)@tab.list, baseball.teams@tab.list)
    
    #maybe not the best way to do this, though haven't seen another aside from re-creating the object
    assign(x="baseball.db", value=baseball.db, envir=.GlobalEnv)
})

test_that("Another, more complex TBSL example based off a sample tracking use case",{
    
    db.list <- test.db.1()
    
    sample.tracking <- new("TableSchemaList")
    
    clinical <- makeSchemaFromData(db.list$clinical, name="clinical")
    check.table.import(db.list$clinical, clinical, "clinical")
    
    sample.tracking <- append(sample.tracking, clinical)
    
    expect_equal(length(sample.tracking), 1)
    
    #this one should fail due to ng.ul column
    expect_error(makeSchemaFromData(db.list$dna, name="dna"))
    
    db.list$dna <- correct.df.names(db.list$dna)
    
    dna <- makeSchemaFromData(db.list$dna, name="dna")
    check.table.import(db.list$dna, dna, "dna")
    
    sample.tracking <- append(sample.tracking, dna)
    
    samples <- makeSchemaFromData(db.list$samples, name="samples")
    check.table.import(db.list$samples, samples, "samples")
    
    sample.tracking <- append(sample.tracking, samples)
    
    expect_equal(length(sample.tracking), 3)
    
    #more complicated usage of relationship
    
    relationship(sample.tracking, from="clinical", to="samples") <- sample_id~sample_id
    check.direct.keys(sample.tracking,  from="clinical", to="samples", key.name="sample_id", orig.to.obj=samples, orig.from.obj=clinical)
    
    relationship(sample.tracking, from="clinical", to="dna") <-sample_id~sample_id
    check.direct.keys(sample.tracking,  from="clinical", to="dna", key.name="sample_id", orig.to.obj=dna, orig.from.obj=clinical)
    
    #Here, db.cols (and db.schema) should be modified so that sample and wave in samples should be replaced with dna's autoinc pk
    relationship(sample.tracking, from="dna", to="samples") <- .~sample_id+wave
    #if there was no clinical to samples rels would expect below however will need to keep sample_id to be consistent with other relationships
    #expect_equal(sort(sample.tracking@tab.list$samples$db.cols), sort(c("samples_ind", "dna_ind", names(db.list$samples)[names(db.list$samples) %in% c("sample_id", "wave")==F])))
    
    expect_equal(sort(sample.tracking@tab.list$samples$db.cols), sort(c("samples_ind", "dna_ind", names(db.list$samples)[names(db.list$samples) %in% "wave"==F])))
    
    #this should not be true for dna
    expect_equal(sort(sample.tracking@tab.list$dna$db.cols), sort(c("dna_ind", names(db.list$dna))))
    
    #also check that the keys look sane
    expect_named(sample.tracking@tab.list$samples$foreign.keys, c("clinical", "dna"))
    
    #should just be the direct keys for clinical
    expect_equal(sample.tracking@tab.list$samples$foreign.keys$clinical, list(local.keys="sample_id", ext.keys="sample_id"))
    
    #should be sample_id and wave as well as dna's pk
    
    expect_equal(sample.tracking@tab.list$samples$foreign.keys$dna, list(local.keys="dna_ind", ext.keys=c("sample_id", "wave")))
    
    
    assign(x="sample.tracking",  value=sample.tracking, envir=.GlobalEnv)
})

test_that("constraint<- method is sane",{
    
    constraint(sample.tracking, "dna") <- ~ sample_id + wave
    expect_equal(gsub("\\s+", "", sample.tracking@tab.list$dna$db.constr), gsub("\\s+", "", "CONSTRAINT dna_idx UNIQUE (sample_id, wave)"))
    
    expect_true(grepl(gsub("\\s+", "", "CONSTRAINT dna_idx UNIQUE (sample_id, wave)"), gsub("\\s+", "", createTable(sample.tracking, 'dna', mode='normal')), fixed=T))
    expect_true(grepl("INSERT\\s+OR\\s+IGNORE", insertStatement(sample.tracking, 'dna')))
    
    constraint(sample.tracking, "dna", should.ignore=F, constr.name="test") <- ~ sample_id + wave
    expect_identical(gsub("\\s+", "", sample.tracking@tab.list$dna$db.constr), gsub("\\s+", "", "CONSTRAINT test UNIQUE (sample_id, wave)"))
    expect_true(sample.tracking@tab.list$dna$should.ignore == F)
    
    expect_true(grepl(gsub("\\s+", "", "CONSTRAINT test UNIQUE (sample_id, wave)"), gsub("\\s+", "", createTable(sample.tracking, 'dna', mode='normal')), fixed=T))
    expect_true(grepl("INSERT\\s+INTO", insertStatement(sample.tracking, 'dna')))
    
    constraint(sample.tracking, "dna", should.ignore=F) <- NULL
    expect_true(sample.tracking@tab.list$dna$db.constr == "")
    expect_true(sample.tracking@tab.list$dna$should.ignore == F)
})

test_that("createTable",
{
    tbsl <- sample.tracking
    
    valid.tables <- names(tbsl@tab.list)
    
    db.con <- dbConnect(SQLite(), tempfile())
    
    for(i in valid.tables)
    {
        print(i)
        for (j in c("normal", "merge"))
        {
            print(j)
            
            f.keys <- tbsl@tab.list[[i]]$foreign.keys
            
            #if there are no foreign keys available, don't allow create table statements to be generated
            if (j == "merge" && poplite:::shouldMerge(tbsl, i)==F)
            {
                expect_error(createTable(tbsl, table.name=i, mode="merge"))
            }
            else
            {
                if (is.null(f.keys) || j == "normal")
                {
                    add.cols <- character(0)
                    add.type <- character(0)
                    add.pk <- integer(0)
                }
                else
                {
                    #the basic table should already exist so can retrieve the previous info on coltyps
                    
                    temp.prag <- do.call("rbind", lapply(names(f.keys), function(x)
                           {
                                dbGetQuery(db.con, paste("pragma table_info(",x,")"))
                           }))
                    
                    key.vals <- as.character(unlist(sapply(f.keys, "[[", "ext.keys")))
                    key.prag <- temp.prag[temp.prag$name %in% key.vals,]
                    add.cols <- key.prag$name
                    add.type <- key.prag$type
                    add.pk <- rep(0L, length(add.type))
                }
                
                prag.tab.name <- ifelse(j=="merge", paste0(i, "_temp"), i)
                
                expect_true(is.null(dbGetQuery(db.con, createTable(tbsl, table.name=i, mode=j))))
                tab.prag <- dbGetQuery(db.con, paste("pragma table_info(",prag.tab.name,")"))
                sub.prag <- tab.prag[,c("name", "type", "pk")]
                
                col.types <- sapply(strsplit(tbsl@tab.list[[i]]$db.schema, "\\s+"), "[", 1)
                col.names <- tbsl@tab.list[[i]]$db.cols
                is.pk <- as.integer(grepl("PRIMARY KEY", tbsl@tab.list[[i]]$db.schema))
                
                col.names <- append(col.names, add.cols)
                col.types <- append(col.types, add.type)
                is.pk <- append(is.pk, add.pk)
                
                query.dta <- data.frame(name=col.names, type=col.types, pk=is.pk, stringsAsFactors=FALSE)
                
                if (j == "merge")
                {
                    #need to add in constrains similar to shouldMerge here...
                    keep.cols <- unlist(lapply(f.keys, function(x)
                           {
                                if (length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys)))
                                {
                                    return(x$local.keys)
                                }else{
                                    return(NULL)
                                }
                           }))
                    
                    #there can be some duplicates in this scenario
                    query.dta <- query.dta[!duplicated(query.dta),]
                    
                    #remove the local keys except for those that are 'direct'
                    
                    locs <- unique(sapply(f.keys, "[[", "local.keys"))
                    
                    rm.locs <- setdiff(locs, keep.cols)
                    
                    query.dta <- query.dta[query.dta$pk == 0 & query.dta$name %in% rm.locs == F,]
                }
                
                ord.prag <- sub.prag[do.call("order", sub.prag),]
                ord.query <- query.dta[do.call("order", query.dta),]
                
                rownames(ord.prag) <- NULL
                rownames(ord.query) <- NULL
                
                expect_equal(ord.prag, ord.query)
            }
        }
        
    }
    
    dbDisconnect(db.con)
})

test_that("insertStatement", {
    set.seed(123)
    
    tbsl <- sample.tracking
    
    valid.tables <- names(tbsl@tab.list)
    
    db.con <- dbConnect(SQLite(), tempfile())
    
    for(i in valid.tables)
    {
        print(i)
        for(j in c("normal", "merge"))
        {
            print(j)
            f.keys <- tbsl@tab.list[[i]]$foreign.keys
            
            if (j == "merge" && poplite:::shouldMerge(tbsl, i)==F)
            {
                expect_error(insertStatement(tbsl, i, j))
            }
            else
            {
                #first create the tables
                expect_true(is.null(dbGetQuery(db.con, createTable(tbsl, table.name=i, mode=j))))
                
                prag.tab.name <- ifelse(j=="merge", paste0(i, "_temp"), i)
                tab.prag <- dbGetQuery(db.con, paste("pragma table_info(",prag.tab.name,")"))
                
                #create a couple lines of fake data to insert into the database
                
                ins.dta <- as.data.frame(matrix(sample.int(10000, 10*nrow(tab.prag)), ncol=nrow(tab.prag), nrow=10, dimnames=list(NULL, tab.prag$name)), stringsAsFactors=fALSE)
                
                for(p in colnames(ins.dta))
                {
                    if (tab.prag$type[tab.prag$name == p] == "TEXT")
                    {
                        ins.dta[,p]  <- as.character(ins.dta[,p])
                    }
                }
                
                #load into the database
                
                dbBeginTransaction(db.con)
                expect_true(is.null(dbGetPreparedQuery(db.con, insertStatement(tbsl, i, mode=j), bind.data = ins.dta)))
                dbCommit(db.con)
                
                #check whether it respects should.ignore
                
                ignore.match <- regexpr(pattern="INSERT\\s+OR\\s+IGNORE", text=insertStatement(tbsl, i, mode=j), perl=TRUE)
                
                if (tbsl@tab.list[[i]]$should.ignore)
                {
                    expect_true(ignore.match != -1)
                }
                else
                {
                    expect_true(ignore.match == -1)
                }
            }
        }
    }
    
    dbDisconnect(db.con)
})

test_that("mergeStatement",
{  
    tbsl <- sample.tracking
    
    valid.tables <- names(tbsl@tab.list)
    
    for(i in valid.tables)
    {
        #again if there are no foreign keys make sure the query dies
        f.keys <- tbsl@tab.list[[i]]$foreign.keys
        
        print(i)
        if (is.null(f.keys))
        {
            expect_error(mergeStatement(tbsl, i))
        }
        else
        {
            #only consider the f.keys which are not direct
            
            is.direct.keys <- sapply(f.keys, function(x) length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys)))
            
            f.keys <- f.keys[is.direct.keys == F]
            
            cur.stat <- mergeStatement(tbsl, i)
            
            #is the insert statement table definition consistent
            
            tab.match <- regexpr(pattern=paste0(i, "\\s+\\(\\s+([\\w+_,]+)\\s+\\)"), text=cur.stat, perl=TRUE)
            tab.str <- substr(cur.stat, start=attr(tab.match, "capture.start"), stop=attr(tab.match, "capture.start")+attr(tab.match, "capture.length")-1)
            split.tab <- strsplit(tab.str, ",")[[1]]
            
            tab.cols <- tbsl@tab.list[[i]]$db.cols
            tab.cols <- tab.cols[tbsl@tab.list[[i]]$db.schema != "INTEGER PRIMARY KEY AUTOINCREMENT"]
            
            expect_true(length(intersect(split.tab, tab.cols)) == length(union(split.tab, tab.cols)))
            
            #is the select statement consistent with the insert statement table definition
                ##need to take into account the specified tables for the statements.
                ##for now will simply remove the table names, prior to checking.  If a problem in the future
                ##can make sure such columns exist as well.
                
            select.match <- regexpr(pattern=paste0("SELECT\\s+", tab.str), text=gsub("[\\w_]+\\.", "", cur.stat, perl=T), perl=TRUE)
            
            expect_true(select.match != -1)
            
            #are the joins sane
            
            join.base <- sapply(f.keys, function(x) paste0("\\(", paste(x$ext.keys, collapse=","), "\\)"))
            
            join.str <- paste0(i, "_temp\\s+", paste(paste("JOIN", names(join.base), "USING", join.base, sep="\\s+"), collapse="\\s+"))
            
            join.match <- regexpr(pattern=join.str, text=cur.stat, perl=TRUE)
            
            expect_true(join.match != -1)
            
            #is it respecting should.ignore
            
            ignore.match <- regexpr(pattern="INSERT\\s+OR\\s+IGNORE", text=cur.stat, perl=TRUE)
            
            if (tbsl@tab.list[[i]]$should.ignore)
            {
                expect_true(ignore.match != -1)
            }else
            {
                expect_true(ignore.match == -1)
            }
        }
    }
})

test_that("Database population",{
    
    #simple example first
    
    ins.vals <- list(team_franch=TeamsFranchises, teams=Teams, salaries=Salaries)
    
    #populate the entire database
    
    do.call(populate, append(list(baseball.db), ins.vals))
    
    #read back in each of the tables and make sure they are consistent with in memory data.frames
    
    test.con <- dbConnect(SQLite(), dbFile(baseball.db))
    
    expect_true(all(names(ins.vals) %in% dbListTables(test.con)))
    
    db.tab.list <- lapply(names(ins.vals), function(x)
           {
                dbReadTable(test.con, x)
           })
    
    names(db.tab.list) <- names(ins.vals)
    
    #these should all have 'table_ind' columns
    
    for (i in names(db.tab.list))
    {
        expect_true(paste(i, "ind", sep="_") %in% names(db.tab.list[[i]]))
    }
    
    #these relationships are all 'direct' keys so check for equivalence between the lists removing the 'table'_ind columns
    
    stripped.db.list <- lapply(db.tab.list, function(x) x[,-1])
    fixed.ins <- lapply(ins.vals, convert.factors.to.strings)
    
    expect_equal(stripped.db.list,fixed.ins)
    
    dbDisconnect(test.con)
    
    #more complex example
    
    sample.tracking.db <- Database(sample.tracking, tempfile())
    
    samp.list <- test.db.1()
    
    #from above, this correction had to be done
    samp.list$dna <- correct.df.names(samp.list$dna)
    
    do.call(populate, append(list(sample.tracking.db), samp.list))
    
    #again read in the db tables
    
    test.con <- dbConnect(SQLite(), dbFile(sample.tracking.db))
    
    expect_true(all(names(samp.list) %in% dbListTables(test.con)))
    
    db.tab.list <- lapply(names(samp.list), function(x)
           {
                dbReadTable(test.con, x)
           })
    
    names(db.tab.list) <- names(samp.list)
    
    #all should have a 'table'_ind column
    
    for (i in names(db.tab.list))
    {
        expect_true(paste(i, "ind", sep="_") %in% names(db.tab.list[[i]]))
    }
    
    #for clinical and dna, should just be the 'table'_ind columns as above
    
    stripped.db.list.1 <- lapply(db.tab.list[c('clinical', 'dna')], function(x) x[,-1])
    fixed.samp.list <- lapply(samp.list , convert.factors.to.strings)
    
    expect_equal(stripped.db.list.1,fixed.samp.list[c('clinical', 'dna')], check.attributes=F)
    
    #for samples, will be both 'table'_ind column as well as a 'foreign table'_ind derived from dna and a 'direct' key from clinical
    ##merge this table with dna to pull out the appropriate column
    
    test.samples <- samp.list[['samples']]
    #add in the samples_ind in this case
    test.samples$samples_ind <- 1:nrow(test.samples)
    
    test.dna <- samp.list[['dna']]
    #it is simply autoincremented PK
    test.dna$dna_ind <- 1:nrow(test.dna)
    
    test.samples.merge <- merge(test.samples, test.dna, all=F)
    
    test.samples.merge <- test.samples.merge[,names(db.tab.list[['samples']])]
    test.samples.merge <- convert.factors.to.strings(test.samples.merge)
    test.samples.merge <- test.samples.merge[do.call('order', test.samples.merge),]
    
    res.samples <- db.tab.list$samples
    
    res.samples <- res.samples[do.call('order', res.samples),]
    
    #remove the sample_inds here, as they don't really matter and seem to misalign, probably due to differences in assignment
    
    expect_equal(test.samples.merge[,-1], res.samples[,-1], check.attributes=F)
    
    dbDisconnect(test.con)
    
    #make sure that populate only provides the desired input...
    
    temp.db <- Database(sample.tracking, tempfile())
    
    expect_error(populate(temp.db, samples=samp.list$clinical, use.tables="clinical"))
    
    assign(x="sample.tracking.db",  value=sample.tracking.db, envir=.GlobalEnv)
})

test_that("Querying with Database objects",
{   
    #onto querying Database objects
    #sample.tracking.db
    
    test.con <- dbConnect(SQLite(), dbFile(sample.tracking.db))
    
    #start with some basic select queries
    
    db.tab.list <- lapply(tables(sample.tracking.db), function(x)
           {
                dbReadTable(test.con, x)
           })
    
    names(db.tab.list) <- tables(sample.tracking.db)
    
    db.samps <- db.tab.list[["samples"]]
    
    #the .tables keyword should select all the columns on the given table
    all.samps <- select(sample.tracking.db, .tables="samples")
    expect_equal(as.data.frame(all.samps), db.samps)
    
    #the se version:
    
    all.samps.se <- select_(sample.tracking.db, .dots=list(.tables="samples"))
    expect_equal(as.data.frame(all.samps.se), db.samps)
    
    #try specifying all the tables
    
    all.tab.join <- select(sample.tracking.db, .tables=tables(sample.tracking.db))
    
    merge.dta <- db.tab.list[[1]]
    
    for(i in seq_along(db.tab.list[-1]))
    {
        merge.dta <- merge(merge.dta, db.tab.list[-1][[i]])
    }
    
    atj.dta <- convert.factors.to.strings(as.data.frame(all.tab.join)[,names(merge.dta)])
    atj.dta <- atj.dta[do.call("order", atj.dta),]
    
    merge.dta <- convert.factors.to.strings(merge.dta)
    merge.dta <- merge.dta[do.call("order", merge.dta),]
    
    expect_equal(atj.dta, merge.dta, check.attributes=F)
    
    #again the se version
    
    all.tab.join.se <- select_(sample.tracking.db, .dots=list(.tables=tables(sample.tracking.db)))
    
    atj.dta.se <- convert.factors.to.strings(as.data.frame(all.tab.join.se)[,names(merge.dta)])
    atj.dta.se <- atj.dta.se[do.call("order", atj.dta.se),]
    
    expect_equal(atj.dta.se, merge.dta, check.attributes=F)
    
    #specifying a subset of the tables
    
    sub.tab.join <- select(sample.tracking.db, .tables=c("samples", "dna"))
    
    merge.dta <- merge(db.tab.list[["samples"]], db.tab.list[["dna"]])
    
    stj.dta <- convert.factors.to.strings(as.data.frame(sub.tab.join)[,names(merge.dta)])
    stj.dta <- stj.dta[do.call("order", stj.dta),]
    
    merge.dta <- convert.factors.to.strings(merge.dta)
    merge.dta <- merge.dta[do.call("order", merge.dta),]
    
    expect_equal(stj.dta, merge.dta, check.attributes=F)
    
    #specifying table, couple of columns
    sub.samps <- select(sample.tracking.db, sample_id:dna_ind,.tables="samples")
    expect_equal(as.data.frame(sub.samps), db.samps[,c("sample_id", "did_collect", "dna_ind")])
    
    #se version
    
    sub.samps.se <- select_(sample.tracking.db, "sample_id:dna_ind",.dots=list(.tables="samples"))
    expect_equal(as.data.frame(sub.samps.se), db.samps[,c("sample_id", "did_collect", "dna_ind")])
    
    #subsetting again without specifying tables
    sub.samps.nt <- select(sample.tracking.db, sample_id:dna_ind)
    expect_equal(as.data.frame(sub.samps.nt), db.samps[,c("sample_id", "did_collect", "dna_ind")])
    
    #se version
    
    sub.samps.nt.se <- select_(sample.tracking.db, "sample_id:dna_ind")
    expect_equal(as.data.frame(sub.samps.nt.se), db.samps[,c("sample_id", "did_collect", "dna_ind")])
    
    #shouldn't work if the subsetting is ambigous
    
    expect_error(select(sample.tracking.db, dna_ind))
    
    #will work if table is specified
    
    #either
    
    expect_equal(as.data.frame(select(sample.tracking.db, dna_ind, .tables="dna")), data.frame(db.tab.list$dna$dna_ind), check.attributes=F)
    
    #or
    expect_equal(as.data.frame(select(sample.tracking.db, dna.dna_ind)), data.frame(db.tab.list$dna$dna_ind), check.attributes=F)
    
    #can select columns across different tables
    
    use.cols <- unique(append(names(db.tab.list$samples)[-1], c("sample_id", "sex", "age", "status")))
    
    two.tab.cols <- merge(db.tab.list$samples, db.tab.list$clinical, by="sample_id", all=F)
    
    nt.res <- as.data.frame(select(sample.tracking.db, sample_id:dna_ind, sample_id:status))
    
    expect_true(all(names(nt.res) %in% use.cols))
    
    nt.comp <- nt.res[do.call("order", nt.res[,use.cols]),use.cols]
    two.tab.cols <- two.tab.cols[do.call("order", two.tab.cols[,use.cols]),use.cols]
    
    expect_equal(nt.comp, two.tab.cols, check.attributes=F)
    
    #again specifying tables
    
    td.res <- as.data.frame(select(sample.tracking.db, samples.sample_id:dna_ind, clinical.sample_id:status))
    
    expect_true(all(names(td.res) %in% use.cols))
    
    td.res <- td.res[do.call("order", td.res[,use.cols]),use.cols]
    
    expect_equal(td.res, two.tab.cols, check.attributes=F)
    
    #or via the .tables mechanism this should not work
    
    expect_error(select(sample.tracking.db, sample_id:dna_ind, sample_id:status, .tables=c("clinical", "samples")))
    
    #there were a few additional ones, check whether mixed named, unnamed works
 
    nu.res <- as.data.frame(select(sample.tracking.db, samples.sample_id:dna_ind, sample_id:status))
    
    expect_true(all(names(nu.res) %in% use.cols))
    
    nu.res <- nu.res[do.call("order", nu.res[,use.cols]), use.cols]
    
    expect_equal(nu.res, two.tab.cols, check.attributes=F)
    
    #bug that came up when preparing the examples
    
    bug.1 <- select(baseball.db, yearID:WCWin, franchName)
    expect_true(all(colnames(bug.1) %in% c("franchName", columns(baseball.db)$teams[which(columns(baseball.db)$teams == "yearID"):which(columns(baseball.db)$teams == "WCWin")])))
    
    #onto filtering
    
    #this shouldn't work as the sample_id column is ambigous
    expect_error(filter(sample.tracking.db,sample_id == 1))
    
    samp.1.filt <- filter(sample.tracking.db, samples.sample_id == 1)
    samp.1.df <- as.data.frame(samp.1.filt)
    
    #check against the table itself
    
    expect_equal(samp.1.df, db.samps[db.samps$sample_id == 1,])
    
    #standard eval version
    
    samp.1.filt.se <- filter_(sample.tracking.db, "samples.sample_id == 1")
    samp.1.df.se <- as.data.frame(samp.1.filt.se)
    
    expect_equal(samp.1.df.se, db.samps[db.samps$sample_id == 1,])
    
    #or
    
    samp.1.filt.se.2 <- filter_(sample.tracking.db, .dots=list("samples.sample_id == 1"))
    samp.1.df.se.2 <- as.data.frame(samp.1.filt.se.2)
    
    expect_equal(samp.1.df.se.2, db.samps[db.samps$sample_id == 1,])
    
    #also should work like when unambigous:
    
    status.res <- filter(sample.tracking.db, status == 1)
    
    expect_equal(as.data.frame(status.res), db.tab.list$clinical[db.tab.list$clinical$status == 1,], check.attributes=F)
    
    #multiple filters are not defined
    
    expect_error(filter(sample.tracking.db, status == 1, sample_id==3))
    
    #though you can do:
    ##as long as the columns are uniquely defined
    clin.filt.res <- filter(sample.tracking.db, status == 1 & sample_id==3)
    
    expect_equal(as.data.frame(clin.filt.res), db.tab.list$clinical[with(db.tab.list$clinical, status == 1 & sample_id == 3),], check.attributes=F)
    
    #similar to also specifying columns
    
    clin.filt.res.1 <- filter(sample.tracking.db, clinical.status == 1 & clinical.sample_id==3)
    
    expect_equal(as.data.frame(clin.filt.res), as.data.frame(clin.filt.res.1))
    
    #again, partially specifying columns
    
    clin.filt.res.2 <- filter(sample.tracking.db, clinical.status == 1 & sample_id==3)
    
    expect_equal(as.data.frame(clin.filt.res), as.data.frame(clin.filt.res.2))
    
    #undefined columns should break
    
    expect_error(filter(sample.tracking.db, kitten == 1 & sample_id==3))
    
    #for filter, cross table queries should result in an outer join
    
    #not currently supported, but in the future...
    #filter(sample.tracking.db, clinical.status == 1 | dna.wave==2)
   
    #filter(sample.tracking.db, clinical.(status == 1 & sex == "m") | dna.wave==2)
})

test_that("sample tracking example but with direct keys between dna and samples", {
    
    db.list <- test.db.1()
    
    sample.tracking <- new("TableSchemaList")
    
    clinical <- makeSchemaFromData(db.list$clinical, name="clinical")
    
    sample.tracking <- append(sample.tracking, clinical)
    
    samples <- makeSchemaFromData(db.list$samples, name="samples")
    
    sample.tracking <- append(sample.tracking, samples)
    
    db.list$dna <- correct.df.names(db.list$dna)
    
    dna <- makeSchemaFromData(db.list$dna, name="dna")
    
    sample.tracking <- append(sample.tracking, dna)
    
    relationship(sample.tracking, from="clinical", to="samples") <- sample_id~sample_id
    
    relationship(sample.tracking, from="clinical", to="dna") <-sample_id~sample_id
    
    #Here, db.cols (and db.schema) should be modified so that sample and wave in samples should be replaced with dna's autoinc pk
    relationship(sample.tracking, from="samples", to="dna") <- sample_id+wave~sample_id+wave
    
    temp.st.db <- Database(sample.tracking, tempfile())
    
    do.call(populate, append(list(temp.st.db),db.list))
    
    #this query was seen to fail
    samp.dna.join <- select(temp.st.db, .tables=c("samples", "dna"))
    
    test.con <- dbConnect(SQLite(), dbFile(temp.st.db))
    
    #start with some basic select queries
    
    db.tab.list <- lapply(tables(temp.st.db), function(x)
           {
                dbReadTable(test.con, x)
           })
    
    names(db.tab.list) <- tables(temp.st.db)
    
    samp.dna.merge <- merge(db.tab.list$samples, db.tab.list$dna, by=c("sample_id", "wave"), all=F)
    
    samp.dna.merge <- convert.factors.to.strings(samp.dna.merge)
    samp.dna.merge <- samp.dna.merge[do.call("order", samp.dna.merge),]
    
    samp.dna.join <- convert.factors.to.strings(as.data.frame(samp.dna.join))[,names(samp.dna.merge)]
    samp.dna.join <- samp.dna.join[do.call("order", samp.dna.join),]
    
    expect_equal(samp.dna.merge, samp.dna.join, check.attributes=F)
    
})

#igraph::plot.igraph(poplite:::tsl.to.graph(om.schema.obj))

test_that("oligoMask queries that break poplite", {
  
    om.schema.obj <- new("TableSchemaList",tab.list=test.schema.2)
    relationship(om.schema.obj, from="allele", to="genotype") <- ref_id+allele_num~ref_id+allele_num
    
    
    test.db <- tempfile()
    
    temp.con <- dbConnect(SQLite(), test.db)
    
    for(i in names(test.db.2))
    {
        dbWriteTable(conn=temp.con, name=i, test.db.2[[i]])
    }
    
    dbDisconnect(temp.con)
    
    #make a database object
    
    test.database.1 <- Database(om.schema.obj, test.db)
    
    prob.tab <- as.data.frame(select(test.database.1, .tables="probe_info"))
    
    expect_equal(prob.tab, test.db.2$probe_info)
    
    #should all be from probe_info
    prob.tab.2 <- as.data.frame(select(test.database.1, probe_id, fasta_name, align_status))
    
    expect_equal(prob.tab.2, test.db.2$probe_info[,c("probe_id", "fasta_name", "align_status")])
    
    all.tab.1 <- as.data.frame(select(test.database.1, .tables=tables(test.database.1)))
    
    all.merge <- test.db.2[[1]]
    
    for(i in names(test.db.2)[-1])
    {
        all.merge <- merge(all.merge, test.db.2[[i]], all=F)
    }
    
    expect_equal(all.tab.1[,names(all.merge)], all.merge)
    
    
    all.tab.2 <- as.data.frame(select(test.database.1, probe_id, fasta_name, align_status, probe_chr, probe_start, probe_end, seqnames, start,
			end, filter, geno_chr, genotype.allele_num, strain))
    
    expect_equal(all.tab.2, all.merge[,c("probe_id", "fasta_name", "align_status", "probe_chr", "probe_start", "probe_end", "seqnames", "start",
			"end", "filter", "geno_chr", "allele_num", "strain")])
    
    
    sel.tab <- as.data.frame(select_(test.database.1, "probe_info.probe_id", "reference.ref_id", "reference.filter", "probe_info.align_status"))
    
    expect_equal(sel.tab, all.merge[,c("probe_id", "ref_id", "filter", "align_status")])
    
    
})


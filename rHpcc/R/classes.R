# The base class that generates the ECL code and executes it on the HPCC cluster
# EXAMPLE R Code:
# ecl1 <- ECL$new(hostName="127.0.0.1")
# recPerson <- ECLRecord$new(name="Person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)

ECL <- setRefClass("ECL", 
  fields = list(hostName = "character", port="character", eclCode = "character", clusterName = "character"),
  methods= list(                   
    add = function(obj) {
      eclCode <<- paste(eclCode, obj$print()) 
    },
                     
    addImport = function(value) {
      eclCode <<- paste(value, eclCode) 
    },
                     
    print = function() { eclCode },
                     
    execute = function() {
      if(length(clusterName) == 0){
        clusterName <<- "thor"
      }
      if(length(port) == 0){
        port <<- "8008"
      }

      result <- eclDirectCall(hostName, port, clusterName, eclCode)    
      result
    }              
  )
)

# Creates an ECL "Record Set" definition. 
# EXAMPLE R Code:
# ecl1 <- ECL$new(hostName="127.0.0.1")
# recPerson <- ECLRecord$new(name="Person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
ECLRecord <- setRefClass("ECLRecord", 
  fields = list(name = "character", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                           
    addField = function(fieldName, fieldValue, seperator="") {
      def <<- paste(def, fieldName, seperator, fieldValue, ";", sep=" ")
    },
                           
    print = function() {
      result <- paste (name, " := RECORD ", def, " END;")
    }                    
  )
)

# Creates a DATASET definition. The DATASET declaration defines a file of records, on disk or in memory.
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1")
# recPerson <- ECLRecord$new(name="Person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# dsPerson <- ECLDataset$new(name="ds_person", datasetType = recPerson, logicalFileName ="~ds::person", fileType="CSV")
# ecl1$add(dsPerson)
ECLDataset <- setRefClass("ECLDataset", 
  fields = list(name = "character", datasetType="ECLRecord", logicalFileName="character", fileType="character", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
    
    getDatasetType = function() {
      result <- datasetType
    },
                            
    addExpression = function(fieldName) {
      def <<- paste(def, fieldName, sep=" ")
    },
                            
    print = function() {
      if(length(def) > 0){
        result <- paste(name, ":=", def, ";",sep=" ")
      } else {
          result <- paste (name, " := DATASET('",logicalFileName,"',", datasetType$getName(), ",", fileType,");")
      }       
    }     
   )
)

# Creates a ECL action "Output". The OUTPUT  action produces a recordset result from the supercomputer, based on which form and options you choose.
# If no file to write to is specified, the result is stored in the workunit and returned to the calling program as a data stream.
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1")
# recPerson <- ECLRecord$new(name="Person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# dsPerson <- ECLDataset$new(name="ds_person", datasetType = recPerson, logicalFileName ="~ds::person", fileType="CSV")
# ecl1$add(dsPerson)
# outputPerson <- ECLOutput$new(name="outputPerson", def = dsPerson$getName())
# ecl1$add(outputPerson)
# xmlContent <- ecl1$execute()
# parseResults(xmlContent)
ECLOutput <- setRefClass("ECLOutput", 
  fields = list(name = "character", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                           
    print = function(){
      result <- paste("OUTPUT(", def, ");", sep="")
    }
  )
)

# Creates an ECL "PROJECT" definition. The PROJECT function processes through all records in the recordset performing the 
# transform function on each record in turn.
# EXAMPLE R CODE:
# ecl <- ECL$new(hostName="127.0.0.1")
# ecl$addImport("IMPORT STD;")
# person <- ECLRecord$new(name="Person")
# person$addField("STRING", "code")
# person$addField("STRING", "firstName")
# person$addField("STRING", "lastName")
# ecl$add(person)
# 
# personOut <- ECLRecord$new(name="PersonOut")
# personOut$addField("STRING", "code")
# personOut$addField("STRING", "firstName")
# personOut$addField("STRING", "lastName")
# ecl$add(personOut)
# 
# personDS <- ECLDataset$new(name="personDS", datasetType = person, logicalFileName ="~ds::person", fileType="CSV")
# ecl$add(personDS)
# 
# personProject <- ECLProject$new(name="PersonProject", inDataset=personDS, outECLRecord=personOut);
# personProject$addField("SELF.firstName", "Std.Str.ToUpperCase(LEFT.firstName)");
# personProject$addField("SELF", "LEFT");
# ecl$add(personProject)
# outputProject <- ECLOutput$new(name="outputProject", def = personProject$getName())
# ecl$add(outputProject)
# ecl$print()
# xmlContent <- ecl$execute()
# data <- parseResults(xmlContent)
# data
ECLProject <- setRefClass("ECLProject", 
  fields = list(name = "character", inDataset="ECLDataset", outECLRecord="ECLRecord", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                            
    addField = function(id, value) {
      def <<- paste(def, id, ":=", value, ";", sep=" ")
    },
                            
    print = function() {
      result <- paste (name, " := PROJECT( ", inDataset$getName(), ", TRANSFORM(", outECLRecord$getName() , ",",def, " ));")
    }     
  )
)

# Creates an ECL "TRANSFORM" definition. A TRANSFORM defines the specific operations that must occur on a record-by-record basis.
# EXAMPLE R CODE:
# transfrm <- ECLTransform$new(name="transfrm", outECLRecord=rec_revenueDef);
# transfrm$addField("SELF.orderNumber", "RIGHT.orderNumber");
# transfrm$addField("SELF.prodCode", "LEFT.productCode");
# transfrm$addField("SELF.prodName", "LEFT.productName");
# transfrm$addField("SELF.revenue", "RIGHT.priceEach * RIGHT.quantityOrdered");
# 
# joinCondition <- "LEFT.productCode=RIGHT.productCode"
# ds_revenue <- ECLJoin$new(name="ds_revenue", leftRecordSet= ds_products, rightRecordSet=ds_orderDetails, joinCondition = joinCondition, joinType = "INNER", def=transfrm$print());
# ecl1$add(ds_revenue)
# output <- ECLOutput$new(name="output", def = ds_revenue$getName())
# ecl$add(output)
# ecl$print()
# xmlContent <- ecl$execute()
# data <- parseResults(xmlContent)
# data

ECLTransform <- setRefClass("ECLTransform", 
  fields = list(name = "character", outECLRecord="ECLRecord", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                              
    addField = function(id, value) {
      def <<- paste(def, id, ":=", value, ";", sep=" ")
    },
                              
    print = function() {
      result <- paste ("TRANSFORM(", outECLRecord$getName() , ",",def, " )")
    }
  )
)

# Creates an ECL "JOIN" definition.
# A inner join if omitted, else one of the listed types in the JOIN Types
# JOIN Types: INNER,LEFT OUTER,RIGHT OUTER,FULL OUTER,LEFT ONLY,RIGHT ONLY,FULL ONLY
# EXAMPLE R CODE:
# transfrm <- ECLTransform$new(name="transfrm", outECLRecord=rec_revenueDef);
# transfrm$addField("SELF.orderNumber", "RIGHT.orderNumber");
# transfrm$addField("SELF.prodCode", "LEFT.productCode");
# transfrm$addField("SELF.prodName", "LEFT.productName");
# transfrm$addField("SELF.revenue", "RIGHT.priceEach * RIGHT.quantityOrdered");
# 
# joinCondition <- "LEFT.productCode=RIGHT.productCode"
# ds_revenue <- ECLJoin$new(name="ds_revenue", leftRecordSet= ds_products, rightRecordSet=ds_orderDetails, joinCondition = joinCondition, joinType = "INNER", def=transfrm$print());
# ecl1$add(ds_revenue)
# output <- ECLOutput$new(name="output", def = ds_revenue$getName())
# ecl$add(output)
# ecl$print()
# xmlContent <- ecl$execute()
# data <- parseResults(xmlContent)
# data
ECLJoin <- setRefClass("ECLJoin", 
  contains="ECLDataset", 
  fields=list(name = "character", leftRecordSet = "ECLDataset", 
              rightRecordSet = "ECLDataset", joinCondition = "character", joinType = "character", def = "character"),
  methods = list(                      
    getName = function() {
      result <- name
    },          
    
    setName = function(value) {
      name <<- value
    },
                         
    print = function() {
      if(length(def) > 0) {                           
        result <- paste (name, " := JOIN( ", leftRecordSet$getName(),",",rightRecordSet$getName(),",", joinCondition, ",", def ,",", joinType ," );")
      } else {
        result <- paste (name, " := JOIN( ", leftRecordSet$getName(),",",rightRecordSet$getName(),",", joinCondition,",", joinType , " );")
      }                           
    }
  )
)

# Creates an ECL "SORT" definition.
# The SORT function sorts the recordset according to the values specified.
# EXAMPLE R CODE: 
# sort <- ECLSort$new(name="sortedTable", inDataset = tblCatalog)
# sort$addField("ProdLine")
# sort$addField("ProdName")
# ecl1$add(sort)
ECLSort <- setRefClass("ECLSort", 
  contains="ECLDataset", 
  fields = list(name = "character", inDataset="ECLDataset", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                         
    addField = function(value) {
      if(length(def) > 0){
        def <<- paste(def, ",")
      }
      def <<- paste(def, value, sep=" ")
    },
                         
    print = function() {
      result <- paste (name, " := SORT( ", inDataset$getName() ,",",def," );")
    }
  )
)

# Creates an ECL "TABLE" definition.
# The TABLE function is similar to OUTPUT, but instead of writing records to a file, it outputs those records in a new
# table (a new dataset in the supercomputer), in memory. The new table is temporary and exists only while the specific
# query that invoked it is running.
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1")
# recPerson <- ECLRecord$new(name="rec_person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# 
# dsPerson <- ECLDataset$new(name="ds_person", datasetType = recPerson, logicalFileName ="~ds::person", fileType="CSV")
# ecl1$add(dsPerson)
# 
# recPersonTable <- ECLRecord$new(name="personNewTableFormat")
# recPersonTable$addField(dsPerson$getName(), "code", seperator=".")
# recPersonTable$addField(dsPerson$getName(), "firstName", seperator=".")
# recPersonTable$addField(dsPerson$getName(), "lastName", seperator=".")
# 
# ecl1$add(recPersonTable)
# 
# tblPerson <- ECLTable$new(name="PersonNewTable", inDataset = dsPerson, format= recPersonTable)
# ecl1$add(tblPerson)

ECLTable <- setRefClass("ECLTable", 
  contains="ECLDataset", fields = list(name = "character", inDataset="ECLDataset", format = "ECLRecord", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                          
    print = function() {
      if(length(def) > 0) {
        result <- paste (name, " := TABLE( ", inDataset$getName() ,",",format$getName(),",", def,");")
      } else {
        result <- paste (name, " := TABLE( ", inDataset$getName() ,",",format$getName()," );")
      }                       
    }
  )
)

# Creates an ECL "DEDUP" definition.
# The DEDUP function evaluates the recordset for duplicate records, as defined by the condition parameter, and returns
# a unique return set. This is similar to the DISTINCT statement in SQL. The recordset should be sorted, unless ALL is specified
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1", port="8008")
# recPerson <- ECLRecord$new(name="rec_person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# 
# dsPerson <- ECLDataset$new(name="ds_person", datasetType = recPerson, logicalFileName ="~ds::person", fileType="CSV")
# ecl1$add(dsPerson)
# 
# recPersonTable <- ECLRecord$new(name="personNewTableFormat")
# recPersonTable$addField(dsPerson$getName(), "code", seperator=".")
# recPersonTable$addField(dsPerson$getName(), "firstName", seperator=".")
# recPersonTable$addField(dsPerson$getName(), "lastName", seperator=".")
# 
# ecl1$add(recPersonTable)
# 
# tblPerson <- ECLTable$new(name="PersonNewTable", inDataset = dsPerson, format= recPersonTable)
# ecl1$add(tblPerson)
# 
# PersonNewTableSorted <- ECLSort$new(name="PersonNewTableSorted", inDataset = tblPerson)
# PersonNewTableSorted$addField("lastName")
# ecl1$add(PersonNewTableSorted)
# 
# mySets <- ECLDedUp$new(name="mySets", inDataset = PersonNewTableSorted)
# mySets$addField("lastName")
# ecl1$add(mySets)
# ecl1$print()
ECLDedUp <- setRefClass("ECLDedUp", 
  fields = list(name = "character", inDataset="ECLDataset", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                          
    addField = function(value) {
      if(length(def) > 0){
        def <<- paste(def, ",")
      }
      def <<- paste(def, value, sep=" ")
    },
                          
    print = function() {
      result <- paste (name, " := DEDUP( ", inDataset$getName() ,",",def," );")
    }
  )
)

# Creates an ECL "ROLLUP" definition.
# The ROLLUP function is similar to the DEDUP function with the addition of the call to the transform function to
# process each duplicate record pair. This allows you to retrieve valuable information from the "duplicate" record before it's thrown away
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1", port="8008")
# recPerson <- ECLRecord$new(name="rec_person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# 
# dsPerson <- ECLDataset$new(name="ds_person", datasetType = recPerson, logicalFileName ="~ds::person", fileType="CSV")
# ecl1$add(dsPerson)
# 
# recPersonContact <- ECLRecord$new(name="rec_myRec")
# recPersonContact$addField(dsPerson$getName(), "firstName", seperator=".")
# recPersonContact$addField(dsPerson$getName(), "lastName", seperator=".")
# 
# ecl1$add(recPersonContact)
# 
# tblPerson <- ECLTable$new(name="LnameTable ", inDataset = dsPerson, format= recPersonContact)
# ecl1$add(tblPerson)
# 
# sort <- ECLSort$new(name="sortedTable", inDataset = tblPerson)
# sort$addField("firstName")
# ecl1$add(sort)
# 
# condition <- "LEFT.firstName = RIGHT.firstName"
# rollUp <- ECLRollUp$new(name="TransformPersons ", inDataset=sort, outECLRecord=recPersonContact, condition = condition);
# rollUp$addField("SELF", "LEFT");
# ecl1$add(rollUp)
# ecl1$print()

ECLRollUp <- setRefClass("ECLRollUp", 
  fields = list(name = "character", inDataset="ECLDataset", outECLRecord="ECLRecord", condition = "character", def = "character"),
  methods = list(                        
    getName = function() {
      result <- name
    },                     
    
    addField = function(id, value) {
      def <<- paste(def, id, ":=", value, ";", sep=" ")
    },
                           
    print = function() {
      result <- paste (name, " := ROLLUP( ", inDataset$getName(),",", condition ,", TRANSFORM(", outECLRecord$getName() , ",",def, " ));")
    }
  )
)

# Creates an ECL "DISTRIBUTE" definition.
# The DISTRIBUTE function re-distributes records from the recordset across all the nodes of the cluster
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="127.0.0.1", port="8008")
# recPerson <- ECLRecord$new(name="rec_person")
# recPerson$addField("STRING", "code")
# recPerson$addField("STRING", "firstName")
# recPerson$addField("STRING", "lastName")
# recPerson$addField("STRING", "address")
# recPerson$addField("STRING", "stateCode")
# recPerson$addField("STRING", "city")
# recPerson$addField("STRING", "zip")
# ecl1$add(recPerson)
# 
# condition <- "SKEW(0.1)"
# distribute <- ECLDistribute$new(inECLRecord=recPerson, condition=condition)
# ecl1$add(distribute)
# ecl1$print()
ECLDistribute <- setRefClass("ECLDistribute", 
  fields = list(name="character", inECLRecord="ECLRecord", condition = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                               
    print = function() {
      if(length(condition) > 0) {
        result <- paste ("DISTRIBUTE( ", inECLRecord$getName(),",", condition ," );")
      } else {
        result <- paste("DISTRIBUTE( ", inECLRecord$getName()," );")
      }
    }
  )
)

# Creates an ECL "ITERATE" definition.
# The ITERATE function processes through all records in the recordset one pair of records at a time, performing the
# transform function on each pair in turn.
# EXAMPLE R CODE:
# ecl1 <- ECL$new(hostName="192.168.217.128", port="8008")
# resType <- ECLRecord$new(name="rec_resType")
# resType$addField("INTEGER1", "Val")
# resType$addField("INTEGER1", "Rtot")
# ecl1$add(resType)
# 
# dsRecords <- ECLDataset$new(name="ds_records", datasetType = resType, logicalFileName ="~ds::iterate", fileType="CSV")
# ecl1$add(dsRecords)
# 
# iterate <- ECLIterate$new(name="ECLIterate", inDataset=dsRecords, outECLRecord=resType);
# iterate$addField("SELF.Rtot", "LEFT.Rtot+RIGHT.Val");
# iterate$addField("SELF", "RIGHT");
# ecl1$add(iterate)
# 
# outputIterate <- ECLOutput$new(name="outputIterate", def = iterate$getName())
# ecl1$add(outputIterate)
# ecl1$print()
# 
# xmlContent <- ecl1$execute()
# parseResults(xmlContent)
ECLIterate <- setRefClass("ECLIterate", 
  fields = list(name = "character", inDataset="ECLDataset", outECLRecord="ECLRecord", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                            
    addField = function(id, value) {
      def <<- paste(def, id, ":=", value, ";", sep=" ")
    },
                            
    print = function() {
      result <- paste (name, " := ITERATE( ", inDataset$getName(),", TRANSFORM(", outECLRecord$getName() , ",",def, " ));")
    }
   )
)

# The TOPN function returns the first count number of records in the sorts order from the recordset.
# topn <- ECLTOPN$new(name="T1", inDataset = dsRecords, count="5")
# topn$addField("-Rtot")
# ecl1$add(iterate)
ECLTOPN <- setRefClass("ECLTOPN", 
  contains="ECLDataset", 
  fields = list(name = "character", inDataset="ECLDataset", count="character", def = "character"),
  methods = list(
    getName = function() {
      result <- name
    },
                         
    addField = function(value) {
      if(length(def) > 0){
        def <<- paste(def, ",")
      }
      def <<- paste(def, " ", value)
    },
                         
    print = function() {
      result <- paste (name, " := TOPN( ", inDataset$getName() ,",",count,",",def," );")
    }
  )
)
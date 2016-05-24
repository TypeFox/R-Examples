library(testthat)
library(XML)
library(rgeos)

# Some functions have different names between GEOS and JTS

process_testxml = function(xmlfile, n, total, descSkip) {
    
    funcTranslate=list( "getboundary"       = list(func=gBoundary,res=readWKT,arg1=readWKT),
                        "getCentroid"       = list(func=gCentroid,res=readWKT,arg1=readWKT),
                        "convexhull"        = list(func=gConvexHull,res=readWKT,arg1=readWKT),
                        "getInteriorPoint"  = list(func=gPointOnSurface,res=readWKT,arg1=readWKT),

                        "isSimple"          = list(func=gIsSimple,res=as.logical,arg1=readWKT),
                        "isValid"           = list(func=gIsValid,res=as.logical,arg1=readWKT),

                        "isWithinDistance"  = list(func=gWithinDistance,res=as.logical,arg1=readWKT,arg2=readWKT,arg3=as.numeric),
                        "intersects"        = list(func=gIntersects,res=as.logical,arg1=readWKT,arg2=readWKT),
                        "contains"          = list(func=gContains,res=as.logical,arg1=readWKT,arg2=readWKT),
                        "within"            = list(func=gWithin,res=as.logical,arg1=readWKT,arg2=readWKT),


                        "intersection"      = list(func=gIntersection,res=readWKT,arg1=readWKT,arg2=readWKT),
                        "union"             = list(func=gUnion,res=readWKT,arg1=readWKT,arg2=readWKT),
                        "difference"        = list(func=gDifference,res=readWKT,arg1=readWKT,arg2=readWKT),
                        "symdifference"     = list(func=gSymdifference,res=readWKT,arg1=readWKT,arg2=readWKT),

                        "relate"            = list(func=gRelate,res=as.logical,arg1=readWKT,arg2=readWKT,arg3=as.character),

                        "covers"            = list(func=gCovers,res=as.logical,arg1=readWKT,arg2=readWKT),
                        "coveredBy"         = list(func=gCoveredBy,res=as.logical,arg1=readWKT,arg2=readWKT))
    
    
    context(paste('(',n,'/',total,')',basename(xmlfile)))
    #x = xmlRoot(xmlTreeParse(I(readLines(xmlfile)),ignoreBlanks=TRUE))
    x = xmlRoot(xmlTreeParse(readLines(xmlfile),ignoreBlanks=TRUE))
    
    nodes = xmlSApply(x,xmlName)

    test_that("valid node types",{
        validNodeTypes = c("precisionModel","case","comment")
        expect_that( all(nodes %in% validNodeTypes), is_true() )
    })


    #Handle precisionModel nodes - only use the first model
    pmAttrs =  xmlAttrs( x[[ which(nodes == "precisionModel")[1] ]] )

    test_that("precisionModel attribute tests", {
        expect_that( length(pmAttrs) == 1 | length(pmAttrs) == 3, is_true() )

        if (length(pmAttrs) == 1) {
            type = pmAttrs[["type"]]
        } else if (length(pmAttrs) == 3) {
            setScale(as.numeric( pmAttrs[["scale"]] ))

            expect_that( pmAttrs[["offsetx"]], equals("0.0") )
            expect_that( pmAttrs[["offsety"]], equals("0.0") )
        } 
    })

    #Handle case nodes
    for ( i in which(nodes == "case") ) {
        caseNodes = xmlSApply(x[[i]],xmlName)

        whichDesc = which(caseNodes == "desc")
        whichTests = which(caseNodes == "test")

        desc = xmlValue( x[[i]][[ whichDesc[1] ]] )
    
        if (desc %in% descSkip)
            next

        whichArgs = which(caseNodes != "desc" & caseNodes != "test")
    
        args = rep( NA,length(whichArgs) )
        # argument nodes can either contain the value or have a file attribute
        for ( j in whichArgs) {
            if (is.null( xmlAttrs(x[[i]][[j]]) )) {
                args[[ xmlName(x[[i]][[j]]) ]] = xmlValue(x[[i]][[j]])
            } else {
                file = xmlAttrs(x[[i]][[j]])[["file"]]
                args[[ xmlName(x[[i]][[j]]) ]] = paste( readLines(file), collapse="" )
            }
        }
    
        #make sure the arg names are lowercase for the sake of consistency
        names(args) = tolower(names(args))
    
        for ( j in whichTests ) {
        
            test_that(paste(desc,'- test nodes in proper format') , {
                expect_that( xmlSize( x[[i]][[j]] ), equals(1) )
                expect_that( xmlName( x[[i]][[j]][[1]] ), equals("op") )
            })
        
            if ( xmlSize( x[[i]][[j]] ) == 1 & xmlName( x[[i]][[j]][[1]] ) == "op" ) {

                opAttrs = xmlAttrs( x[[i]][[j]][[1]] )
                opReturn = xmlValue( x[[i]][[j]][[1]] )
                opNArgs = length(opAttrs)-1

                # some ops seem to have a pattern argument that is not used
                if ( 'pattern' %in% names(opAttrs) )
                    opNArgs = opNArgs-1

                opName = opAttrs[['name']]
            
               
            
                test_that(paste(desc,'-',opName), {
                
                    funcdetails = funcTranslate[[opName]]
                    expect_that( is.null(funcdetails), is_false() )
            
                    if ( !is.null(funcdetails) ) {
                        funcNArgs = length( funcdetails )-2
                        expect_that(funcNArgs==opNArgs, is_true())
                    
                        funcArgs = list()
                        for (k in 1:funcNArgs) {
                            argName = paste("arg",k,sep='')

                            argVal = tolower(opAttrs[[argName]])
                            if (argVal %in% names(args))
                                argVal = args[[ argVal ]]
                            funcArgs[k] =  funcdetails[[argName]](argVal)    
                        }
                    
                        funcReturn = do.call(funcdetails[["func"]], funcArgs)
                        expectedReturn = funcdetails[["res"]](opReturn)
                    
                        if (is.logical(funcReturn)) {
                            expect_that(funcReturn == expectedReturn, is_true())
                        } else if (is.null(funcReturn)) {
                            expect_that(is.null(funcReturn) & is.null(expectedReturn), is_true())
                        } else if (gIsEmpty(expectedReturn)) {
                            expect_that(identical(funcReturn,expectedReturn), is_true())
                        } else { # if it isn't logical or NULL it should be a geometry
                            expect_that(gEquals(funcReturn,expectedReturn),is_true())
                        }
                    }
                }) 
            }                
        } 
    }
}
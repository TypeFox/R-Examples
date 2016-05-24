source("process_testxml.R")

# xml files to skip
fileSkip = c("TestValid2.xml","ExternalRobustness.xml","TestRobustOverlayFloat.xml","TestRobustOverlayFixed.xml")

# tests to skip
descSkip = c(# TestFunctionLAPrec.xml
             "LA - line and sliver intersecting, dimensional collapse",
             # TestFunctionAAPrec.xml
             "AA - intersecting slivers, dimensional collapse",
             "AA - sliver triangle with multiple intersecting boxes",
             "AA - sliver triangles, at angle to each other",
             "AA - B sliver crossing A triangle in line segment with length < 1",
             "AA - A hole close to shell, B coincident with A",
             "AA - hole close to shell, B coincident with A",
             "AA - sliver triangle, cut by polygon",
             "AA - polygon with outward sliver, cut by polygon",
             "AA - hole close to shell",
             "AA - A sliver triangle cutting all the way across B",
             "AA - A polygon with sliver cutting all the way across B",
             "AA - B hole close to shell, A coincident with B",
             # TestInteriorPoint.xml
             "L - linestring with single segment")

xmldir = 'tests/testxml'
testdirs = list.files(system.file(xmldir,package="rgeos"))

for (d in testdirs) {
    testfiles = list.files(system.file(file.path(xmldir,d),package="rgeos")) 
    
    for (f in testfiles)  {

        # files to skip
        if (f %in% fileSkip)
            next

        xmlfile =  system.file(file.path(xmldir,d,f),package="rgeos")

        n = which(f==testfiles)
        total = length(testfiles)
        
        process_testxml(xmlfile, n, total, descSkip)
    }
}


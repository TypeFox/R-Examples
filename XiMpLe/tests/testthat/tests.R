context("XML generation")

test_that("generate empty XML node", {
	sampleXMLStandard <- dget("sample_XML_node_empty_dput.txt")
	expect_that(
		XMLNode("empty"),
		equals(sampleXMLStandard)
	)
})

test_that("generate closed XML node", {
	sampleXMLStandard <- dget("sample_XML_node_closed_dput.txt")
	expect_that(
		XMLNode("empty", ""),
		equals(sampleXMLStandard)
	)
})

test_that("generate closed XML node with attributes", {
	# re-create object sampleXMLnode.attrs
	load("sample_XML_node_attrs.RData")
	expect_that(
		XMLNode("empty", "test", attrs=list(foo="bar")),
		equals(sampleXMLnode.attrs)
	)
})

test_that("generate nested XML tag tree", {
	# re-create object sampleXMLTree
	load("sample_XML_tree.RData")

	sampleXMLnode.empty <- XMLNode("empty")
	sampleXMLnode.closed <- XMLNode("empty", "")
	sampleXMLnode.attrs <- XMLNode("empty", "test", attrs=list(foo="bar"))
	sampleXMLTree.test <- XMLTree(
		XMLNode("tree",
			sampleXMLnode.empty,
			sampleXMLnode.closed,
			sampleXMLnode.attrs
		)
	)

	expect_that(
		sampleXMLTree.test,
		equals(sampleXMLTree)
	)
})


context("XML parsing")

test_that("parse XML file", {
	# re-create object sampleXMLparsed
	load("sample_RSS_parsed.RData")

	sampleXMLFile <- normalizePath("koRpus_RSS_sample.xml")
	XMLtoParse <- file(sampleXMLFile, encoding="UTF-8")
	sampleXMLparsed.test <- parseXMLTree(XMLtoParse)
	close(XMLtoParse)

	expect_that(
		sampleXMLparsed.test,
		equals(sampleXMLparsed))
})


context("extracting nodes")

test_that("extract node from parsed XML tree", {
	# re-create object sampleXMLparsed
	load("sample_RSS_parsed.RData")
	# re-create object sampleXMLnode.extracted
	load("sample_XML_node_extracted.RData")

	sampleXMLnode.test <- node(sampleXMLparsed, node=list("rss","channel","atom:link"))

	expect_that(
		sampleXMLnode.test,
		equals(sampleXMLnode.extracted))
})


context("changing node values")

test_that("change attribute values in XML node", {
	# re-create object sampleXMLparsed
	load("sample_RSS_parsed.RData")
	# re-create object sampleXMLnode.extracted
	load("sample_XML_tree_changed.RData")

	# replace URL
	node(sampleXMLparsed,
		node=list("rss","channel","atom:link"),
		what="attributes", element="href") <- "http://example.com"

	# remove "rel" attribute
	node(sampleXMLparsed,
		node=list("rss","channel","atom:link"),
		what="attributes", element="rel") <- NULL

	expect_that(
		sampleXMLparsed,
		equals(sampleXMLparsed.changed))
})

test_that("change nested text value in XML node", {
	# re-create object sampleXMLparsed
	load("sample_RSS_parsed.RData")
	# re-create object sampleXMLnode.extracted
	load("sample_XML_tree_changed_value.RData")

	# change text
	node(sampleXMLparsed,
		node=list("rss","channel","item","title"),
		what="value",
		cond.value="Changes in koRpus version 0.04-30") <- "this value was changed!"

	expect_that(
		sampleXMLparsed,
		equals(sampleXMLparsed.changed.value))
})

context("getter/setter methods")

test_that("set and get XML node name", {
	sampleXMLnode <- XMLNode("name", attrs=list(atr="test"))
	# set node name
	XMLName(sampleXMLnode) <- "changed"
	sampleXMLnode.name <- XMLName(sampleXMLnode)

	expect_that(
		sampleXMLnode.name,
		equals("changed"))
})

test_that("set and get XML node attributes", {
	sampleXMLnode <- XMLNode("name", attrs=list(atr="test"))
	# set node attributes
	XMLAttrs(sampleXMLnode) <- list()
	sampleXMLnode.attrs <- XMLAttrs(sampleXMLnode)

	expect_that(
		sampleXMLnode.attrs,
		equals(list()))
})

test_that("set and get XML node text value", {
	sampleXMLnode <- XMLNode("name", attrs=list(atr="test"))
	# set node name
	XMLValue(sampleXMLnode) <- "added value"
	sampleXMLnode.value <- XMLValue(sampleXMLnode)

	expect_that(
		sampleXMLnode.value,
		equals("added value"))
})

test_that("set and get XML node children", {
	sampleXMLnode <- XMLNode("name", attrs=list(atr="test"))
	# children will be returned as a list
	sampleXMLChild <- list(XMLNode("child", attrs=list(atr="test")))
	# set node children
	XMLChildren(sampleXMLnode) <- sampleXMLChild
	sampleXMLnode.children <- XMLChildren(sampleXMLnode)

	expect_that(
		sampleXMLnode.children,
		equals(sampleXMLChild))
})

test_that("set and get XML tree children", {
	load("sample_XML_tree.RData")
	# children will be returned as a list
	sampleXMLChild <- list(XMLNode("child", attrs=list(atr="test")))
	# set node children
	XMLChildren(sampleXMLTree) <- sampleXMLChild
	sampleXMLTree.children <- XMLChildren(sampleXMLTree)

	expect_that(
		sampleXMLTree.children,
		equals(sampleXMLChild))
})

test_that("set and get XML tree DTD info", {
	load("sample_XML_tree.RData")
	sampleDTD <- list(doctype="html", decl="PUBLIC",
		id="-//W3C//DTD XHTML 1.0 Transitional//EN",
		refer="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd")
	# set missing values
	XMLDTD(sampleXMLTree) <- sampleDTD
	# try to get them back
	sampleXMLTree.DTD <- XMLDTD(sampleXMLTree)

	expect_that(
		sampleXMLTree.DTD,
		equals(sampleDTD))
})

test_that("set and get XML tree decl info", {
	load("sample_XML_tree.RData")
	sampleDecl <- list(version="1.0", encoding="UTF-8")
	# set missing values
	XMLDecl(sampleXMLTree) <- sampleDecl
	# try to get them back
	sampleXMLTree.decl <- XMLDecl(sampleXMLTree)

	expect_that(
		sampleXMLTree.decl,
		equals(sampleDecl))
})

test_that("set and get XML tree file info", {
	load("sample_XML_tree.RData")
	# set missing values
	XMLFile(sampleXMLTree) <- "somefile.xml"
	# try to get them back
	sampleXMLTree.file <- XMLFile(sampleXMLTree)

	expect_that(
		sampleXMLTree.file,
		equals("somefile.xml"))
})


context("getter/setter methods: XMLScan")

test_that("scan XML tree for node names", {
	load("sample_XML_tree.RData")

	# this should return a list of 3
	sampleXMLTree.nodes <- XMLScan(sampleXMLTree, "empty")

	expect_is(
		sampleXMLTree.nodes,
		"list")
	expect_that(
		length(sampleXMLTree.nodes),
		equals(3))
})

test_that("remove XML nodes from tree by name", {
	load("sample_XML_tree.RData")

	# this should remove all nodes,
	# exept the parent "tree" node
	XMLScan(sampleXMLTree, "empty") <- NULL

	expect_identical(
		sampleXMLTree,
		XMLTree(XMLNode("tree")))
})

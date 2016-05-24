get.metacyc.all <-
function() {
	url = 'http://websvc.biocyc.org/xmlquery?[x:x<-meta^^reactions]'

  message('Download MetaCyc reaction list. It takes a while')
	h = basicTextGatherer()
	
  curl.result = tryCatch({
	  curlPerform(url = url, writefunction = h$update)
	}, warning = function(w) {
	  message('WARNING: MetaCyc server is unstable now. We are trying to get alternative server, but we recommend try this function again later')
	}, error = function(e) {
	  message('WARNING: MetaCyc server is unstable now. We are trying to get alternative server, but we recommend try this function again later')
    alter_url = 'http://compbio.korea.ac.kr/rbiorxn/metacyc_all.html'
	  curlPerform(url = alter_url, writefunction = h$update)
	})
	
	xml = h$value()
  if(xml == '') {
    message('get.metacyc.all() function is not properly performed. Please email the author if it keeps happening')
    return(NA)
  }
  
	xml = unlist(strsplit(xml, '\n'))
	
	index = grep('^  <Reaction ID', xml)
	reg_ex = "(.*META:)(.*)(' orgid.+)"
	reactionIds = sub(reg_ex, '\\2', xml[index])
  
  message(sprintf('%s reaction entries are being downloaded', length(reactionIds)))
	metacycDf = .parse.metacyc.biopax(reactionIds)
	return(metacycDf)
}
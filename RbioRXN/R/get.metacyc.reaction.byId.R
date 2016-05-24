get.metacyc.reaction.byId <-
function(metacycId) {
    if(length(metacycId) == 0) {
        message('Please enter more than one MetaCyc ID')
    }
  
    result_df = c()
    for(i in metacycId) {
        cat('processing',i,'\n')
        urlBase = 'http://websvc.biocyc.org/META/pathway-biopax?type=3&object=%s'
        url = sprintf(urlBase, i)
    
        h = basicTextGatherer()
        
        curl.result = tryCatch({
          curlPerform(url = url, writefunction = h$update)
        }, warning = function(w) {
          message('WARNING: MetaCyc server is unstable now. We are trying to get alternative server, but we recommend try this function again later')
        }, error = function(e) {
          message('WARNING: MetaCyc server is unstable now. We are trying to get alternative server, but we recommend try this function again later')
          alter_urlBase = 'http://compbio.korea.ac.kr/rbiorxn/%s.html'
          alter_url = sprintf(alter_urlBase, i)
          curlPerform(url = alter_url, writefunction = h$update)
        })
          
        if(h$value() == '') {
          return(NA)
        } else if(h$value() != '') {
          biopax = h$value()
          biopax = unlist(strsplit(biopax, '\n'))
          
          parsedBiopax = tryCatch({
            .parse.biopax(biopax)
          }, error = function(cond) {
            message(sprintf('WARNING: RbioRXN could not parse reaction %s. It\'s going to be empty data frame', metacycId))
            MetaCyc = i
            result = data.frame(MetaCyc, stringsAsFactors=F)
            return(result)
          }, warning = function(cond) {
            message(sprintf('WARNING: RbioRXN could not parse reaction %s. It\'s going to be empty data frame', metacycId))
            MetaCyc = i
            result = data.frame(MetaCyc, stringsAsFactors=F)
            return(result)
          })
          result_df = rbind.fill(result_df, data.frame(parsedBiopax, stringsAsFactors=F))
        }
    }
    result_df[is.na(result_df)] = ''
    return(result_df)
}

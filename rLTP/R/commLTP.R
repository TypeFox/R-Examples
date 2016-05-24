commLTP = function(input,mission='ws',api_key)
{
    addr = 'http://ltpapi.voicecloud.cn/analysis/'
    if(mission=='ws')
        return_format = 'plain'
    else
        return_format = 'json'
    result = postForm(uri = addr,
                      style = "post",
                      'text' = input,
                      'api_key' = api_key,
                      'pattern' = mission,
                      'format' = return_format)

    Encoding(result)='UTF-8'
    if (return_format!='plain')
    {
        if (Sys.info()['sysname']=='Windows')
            result = iconv(result,'UTF-8','GBK')
        return(result)
    }
    
    result = parseLTP(result)
    result
}


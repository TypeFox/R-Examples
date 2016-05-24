getAccessToken <-
function(client.id, client.secret){
    fields <- list(
        client_id = client.id,
        client_secret = client.secret,
        scope = 'http://api.microsofttranslator.com',
        grant_type = 'client_credentials'
        )

    return(
        fromJSON(postForm('https://datamarket.accesscontrol.windows.net/v2/OAuth2-13',
                          .params = fields,
                          style = 'POST'))[['access_token']]
        )
}

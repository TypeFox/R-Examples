
irtDemo <- function(text){

    demo <- switch(text,
                   dich    = '113ae624d4bb0a9c18d49cb11d516add',
                   eapMap  = 'ac2031709ccfe7281057baf193f7d611',
                   est2PL  = 'b4bacc35ce94136e9a792fb0044cbf16',
                   est3PL  = '70bebb68ad43a05632742cd5dc6e5667',
                   mle     = '49bdd572e430aa7435e0564341e2aa42',
                   mirt    = '2488b4fa4487481d14ae63bae66601a5',
                   gpcm    = 'f92ce43157dd5bb7b5cffe62fc04b149',
                   grm     = 'a818a2af7d5fe8f9ec726fce81677f1f',
                   grsm    = '0038e50f7c924b5d34ac9487c11c46c7',
                   nrm     = 'c2ee5461e11eddebc76bc3a6d37fb3e8',
                   stop(paste("Enter one of the following: \n 'mle' \n 'est2pl' \n 'est3pl' \n 'eapMap' \n 'dich' \n 'gpcm' \n 'grm' \n 'grsm' \n 'nrm' \n 'mirt' "))
    )
    shiny::runGist(paste(demo))
}



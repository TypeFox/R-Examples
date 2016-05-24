library(nws)
host = 'localhost'
port = 8765
wsname = 'hello'

nws = netWorkSpace(wsname, host, port)

count = 10
cat('hello: iterations:', count, '\n')
nwsStore(nws, 'hello example', count)

for (i in 1:count) {
  nwsStore(nws, 'hello', i)
  j = nwsFetch(nws, 'hello')
}

nwsFetch(nws, 'hello example')
cat('Success\n')

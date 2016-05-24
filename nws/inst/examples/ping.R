library(nws)
host = 'localhost'
port = 8765
wsname = 'ping-pong'

nws = netWorkSpace(wsname, host, port)

nwsStore(nws, 'game', 0)

cat('Ping-pong server ', wsname, ' starting\n')
while (TRUE) {
  pong = nwsFetch(nws, 'ping')
  cat('Got a ping from ', pong, '\n')
  nwsStore(nws, pong, 'pong')
}

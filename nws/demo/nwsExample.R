ws = netWorkSpace('r place')
cat('connected, listing contents of netWorkSpace (should be nothing there).\n', nwsListVars(ws), '\n')

nwsStore(ws, 'x', 1)
cat('should now see x.\n', nwsListVars(ws), '\n')

cat('nwsFind (but don\'t consume) x.\n', nwsFind(ws, 'x'), '\n')
cat('check that it is still there.\n', nwsListVars(ws), '\n')

cat('associate another value with x.\n')
nwsStore(ws, 'x', 2)
cat(nwsListVars(ws), '\n')

cat('consume values for x, should see them in order saved.\n',
    nwsFetch(ws, 'x'), '\n', nwsFetch(ws, 'x'), '\n')
cat('no more values for x... .\n', nwsListVars(ws), '\n')

cat('so try to nwsFetch and see what happens... .\n',
    nwsFetchTry(ws, 'x', 'no go'), '\n')

cat('create a single-value variable.\n')
nwsDeclare(ws, 'pi', 'single')
cat(nwsListVars(ws), '\n')

cat('get rid of x.\n')
nwsDeleteVar(ws, 'x')
cat(nwsListVars(ws), '\n')

cat('try to nwsStore two values to pi.\n')
nwsStore(ws, 'pi', 2.171828182)
nwsStore(ws, 'pi', 3.141592654)
cat(nwsListVars(ws), '\n')

cat('check that the right one was kept.\n', nwsFind(ws, 'pi'), '\n')

cat('what about the rest of the world?', nwsListWss(ws@server), '\n')


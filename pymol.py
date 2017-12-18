# import aaindex
# aaindex.init(path='.')
#
# aaindex.grep('volume')
# x = aaindex.get('KRIW790103')
# print x
# print x.get('A')
#
# aaindex.grep('blosum')
# x = aaindex.get('HENS920102')
# print x.get('A', 'K')


import aaindex
aaindex.init(path='.')

aaindex.grep('volume')
x = aaindex.get('KRIW790103')
print(x)
print(x.get('L'))

#aaindex.grep('blosum')
#x = aaindex.get('HENS920102')
#print(x.get('A', 'K'))

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


def getIndexChar(char):
    import aaindex
    aaindex.init(path='.')

    dic = aaindex._aaindex.values()
    listIndices =[]
    for item in dic:
        # print(item.key)
        listIndices.append(item.key)
    listIndices= ['KARS160114','KARS160106','MIYS990101','CASG920101','FAUJ880111','FAUJ880112']

    listComputedValue = []
    # forbiddenList=['KARS160114','TANS760101','MONM990201','NADH010101']

    for indice in listIndices:
        # if indice not in forbiddenList:
        try:

            x = aaindex.get(indice)
            listComputedValue.append(x.get(char))
            print(indice)
        except Exception as e:
            # print(indice)
            # print("")
            # l=2
            print("ECHEC AVEC INDICES : " +indice)
    return listComputedValue

def getIndicesUsed():
    return  ['KARS160106','MIYS990101','CASG920101','FAUJ880111','FAUJ880112']
    # return  ['KARS160114','KARS160106','MIYS990101','CASG920101','FAUJ880111','FAUJ880112']

def getIndexSeq(seq):
    import aaindex
    aaindex.init(path='.')

    dic = aaindex._aaindex.values()
    listIndices= getIndicesUsed()

    listComputedValue = [0]*len(listIndices)

    forbiddenList=['KARS160114','TANS760101','MONM990201','NADH010101']
    for char in seq:
        for i in range(len(listIndices)):
            indice = listIndices[i]
            try:
                x = aaindex.get(indice)
                listComputedValue[i]+=(x.get(char))
            except Exception as e:
                print("ECHEC AVEC INDICES : " +indice)

    return listComputedValue



print(getIndexSeq("MESSKKMDSPGALQTNPPLKLHTDRSAGTPVFVPEQGGYKEKFVKTVEDKYKCEKCHLVLCSPKQTECGHRFCESCMAALLSSSSPKCTACQESIVKDKVFKDNCCKREILALQIYCRNESRGCAEQLMLGHLLVHLKNDCHFEELPCVRPDCKEKVLRKDLRDHVEKACKYREATCSHCKSQVPMIALQKHEDTDCPCVVVSCPHKCSVQTLLRSELSAHLSECVNAPSTCSFKRYGCVFQGTNQQIKAHEASSAVQHVNLLKEWSNSLEKKVSLLQNESVEKNKSIQSLHNQICSFEIEIERQKEMLRNNESKILHLQRVIDSQAEKLKELDKEIRPFRQNWEEADSMKSSVESLQNRVTELESVDKSAGQVARNTGLLESQLSRHDQMLSVHDIRLADMDLRFQVLETASYNGVLIWKIRDYKRRKQEAVMGKTLSLYSQPFYTGYFGYKMCARVYLNGDGMGKGTHLSLFFVIMRGEYDALLPWPFKQKVTLMLMDQGSSRRHLGDAFKPDPNSSSFKKPTGEMNIASGCPVFVAQTVLENGTYIKDDTIFIKVIVDTSDLPDP"))


# print(aaindex._aaindex.values())

# aaindex.grep('volume')
# x = aaindex.get('KRIW790103')
# print(x)
# print(x.get('L'))

#aaindex.grep('blosum')
#x = aaindex.get('HENS920102')
#print(x.get('A', 'K'))

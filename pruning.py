f = open("C:/Users/escroc/Documents/projectBioInformatique/Supp-B.txt", 'r')
#
# while(f.read()):
#     print(text)
data =""
for line in f:
    if not line[0].isdigit():
        data+=line

print(data)

save = f = open("C:/Users/escroc/Documents/projectBioInformatique/Supp-B-prunned.txt", 'w')
save.write(data)

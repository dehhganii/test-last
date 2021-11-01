file1 = open('tmp_output', 'r')
Lines = file1.readlines()
address_BF = "BF_output"
address_KMP = "KMP_output"
address_BM = "BM_output" 
count = 0
for line in Lines:
	if count %3 == 0:
		f = open(address_BF, "a")
		f.write(str(line))
		f.close()
	if count%3 == 1:
		f = open(address_KMP, "a")
		f.write(str(line))
		f.close()
	if count%3 == 2:
		f = open(address_BM, "a")
		f.write(str(line))
		f.close()
	count +=1




def pKvalue (aa):
	if (aa == 'C'):
		return 8.18
	elif (aa == 'D'):
		return 3.65
	elif (aa == 'E'):
		return 4.25
	elif (aa == 'H'):
		return 6
	elif (aa == 'K'):
		return 10.53
	elif (aa == 'R'):
		return 12.48
	elif (aa == 'Y'):
		return 10.07
	elif (aa == 'c'):
		return 2.34
	elif (aa == 'n'):
		return 9.69
	else:
		return 0

#sample usage
#print(calculateCharge(7.0, "FLPVLAGLTPSIVPKLVCLLTKKC"))
#2.87315
def calculateCharge(pH, seq):
	charge = 0.0;

	# n-terminus
	charge += (1 / (1 + pow(10,(1 * (pH - pKvalue('n'))))))

	# Amino acid sequence
	for aminoAcid in seq:

		if (aminoAcid == 'R'):
			charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('R'))))))
		elif (aminoAcid == 'H'):
			charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('H'))))))
		elif (aminoAcid == 'K'):
			charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('K'))))))
		elif (aminoAcid == 'D'):
			charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('D'))))))
		elif (aminoAcid == 'E'):
			charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('E'))))))
		elif (aminoAcid == 'C'):
			charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('C'))))))
		elif (aminoAcid == 'Y'):
			charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('Y'))))))

	# c-terminus
	charge += (-1 / (1 + pow(10,(-1 * (pH - pKvalue('c'))))));

	return charge;

def calculateChargeAminoAcid(pH, aminoAcid):
	charge = 0.0;

	# n-terminus
	charge += (1 / (1 + pow(10,(1 * (pH - pKvalue('n'))))))

	# Amino acid sequence

	if (aminoAcid == 'R'):
		charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('R'))))))
	elif (aminoAcid == 'H'):
		charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('H'))))))
	elif (aminoAcid == 'K'):
		charge += (1 / (1 + pow(10, (1 * (pH - pKvalue('K'))))))
	elif (aminoAcid == 'D'):
		charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('D'))))))
	elif (aminoAcid == 'E'):
		charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('E'))))))
	elif (aminoAcid == 'C'):
		charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('C'))))))
	elif (aminoAcid == 'Y'):
		charge += (-1 / (1 + pow(10, (-1 * (pH - pKvalue('Y'))))))

	# c-terminus
	charge += (-1 / (1 + pow(10,(-1 * (pH - pKvalue('c'))))));

	return charge;

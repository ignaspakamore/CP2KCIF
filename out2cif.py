import sys
import numpy as np

class CP2KCIF():
	def __init__(self, file):
		self.file = file

	def get_coords(slef):
		parsing = False
		data =[]
		elements = [] 
		
		for line in self.file:
			line = line.strip()
			if line.startswith('&COORD'):
				parsing = True
			if line.startswith('&END COORD'):
				parsing = False

			if line.startswith('UNIT'):
				parsing = False

			if parsing and not line.startswith('&COORD'):
				data.append(line)
		coord = np.zeros((len(data), 4))
		for i, j in enumerate(data):
			j = j.split()
			elements.append(j[0])
			coord[i][0] = j[1] 
			coord[i][1] = j[2]
			coord[i][2] = j[3]

		return coord, elements

	def get_uc(self):
		parsing = False
		data = []

		for line in self.file:
			line = line.strip()
			if line.startswith('&CELL'):
				parsing = True
			if line.startswith('&END CELL'):
				parsing = False
			if line.startswith('MULTIPLE_UNIT_CELL'):
				parsing = False
			if parsing and not line.startswith('&CELL'):
				data.append(line)

		matrix = np.zeros((len(data), 3))
		for i, j in enumerate(data):
			j = j.split()
			matrix[i][0] = j[1]
			matrix[i][1] = j[2]
			matrix[i][2] = j[3]
		return matrix


	def frac2cart(self):
		pass
	def cart2frac(self):
		pass
	def gen_cif(self):
		pass
		

if __name__ == '__main__':
	inpt = sys.argv[1]

	f = open(inpt, 'r')
	
	cp2cif = CP2KCIF(f)
	cp2cif.get_uc()

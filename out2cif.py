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
	def det3(mat):
		return ((mat[0][0]*mat[1][1]*mat[2][2]) + (mat[0][1]*mat[1][2]*mat[2][0]) + (mat[0][2]*mat[1][0]*mat[2][1]) - (mat[0][2]*mat[1][1]*mat[2][0]) - (mat[0][1]*mat[1][0]*mat[2][2]) - (mat[0][0]*mat[1][2]*mat[2][1]))


	def frac2cart(cellParam, fracCoords):
		'''
		 a function that takes the cell parameters, in angstrom, and a list of fractional coordinates
		 and returns the structure in Cartesian coordinates
		
		 Assumes the cellParam matrix is of the form:
		 | ax  ay  az |
		 | bx  by  bz |
		 | cx  cy  cz |
		
		'''
	  cartCoords = []
	  for i in fracCoords:
	    xPos = i[1]*cellParam[0][0] + i[2]*cellParam[1][0] + i[3]*cellParam[2][0]
	    yPos = i[1]*cellParam[0][1] + i[2]*cellParam[1][1] + i[3]*cellParam[2][1]
	    zPos = i[1]*cellParam[0][2] + i[2]*cellParam[1][2] + i[3]*cellParam[2][2]
	    cartCoords.append([i[0], xPos, yPos, zPos])
	  return cartCoords


	def cart2frac(cellParam, cartCoords):
		#####################################
		# a function that takes the cell parameters, in angstrom, and a list of Cartesian coordinates
		# and returns the structure in fractional coordinates
		#
		# Uses Cramer's Rule to solve for the fractional coordinates
		#
		# Assumes the cellParam matrix is of the form:
		# | ax  ay  az |
		# | bx  by  bz |
		# | cx  cy  cz |
		#
		# Need to use the transpose of this matrix in calculation, so call transpose function first
		#
		# Assumes cartCoords are of the form:
		# | X  x0  y0  z0 |
		# | X  x1  y1  z1 |
		# | X  x2  y2  z2 |
		# | ............. |
		# where X is the element symbol
		
	  latCnt = [x[:] for x in [[None]*3]*3]
	  for a in range(3):
	    for b in range(3):
		  latCnt[a][b] = cellParam[b][a]

	  fracCoords = []
	  detLatCnt = det3(latCnt)
	  for i in cartCoords:
	    aPos = (det3([[i[1], latCnt[0][1], latCnt[0][2]], [i[2], latCnt[1][1], latCnt[1][2]], [i[3], latCnt[2][1], latCnt[2][2]]])) / detLatCnt
	    bPos = (det3([[latCnt[0][0], i[1], latCnt[0][2]], [latCnt[1][0], i[2], latCnt[1][2]], [latCnt[2][0], i[3], latCnt[2][2]]])) / detLatCnt
	    cPos = (det3([[latCnt[0][0], latCnt[0][1], i[1]], [latCnt[1][0], latCnt[1][1], i[2]], [latCnt[2][0], latCnt[2][1], i[3]]])) / detLatCnt
	    fracCoords.append([i[0], aPos, bPos, cPos])
	  return fracCoords

	def gen_cif(self):
		pass
		

if __name__ == '__main__':
	inpt = sys.argv[1]

	f = open(inpt, 'r')
	
	cp2cif = CP2KCIF(f)
	cp2cif.get_uc()

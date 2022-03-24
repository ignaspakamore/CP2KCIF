import sys
import numpy as np
import math

class CP2KCIF():
	def __init__(self, file):
		self.file = file
		self.unit = ''
		self.scaled = ''
		self.method = ''
		if self.scaled == 'T' or self.scaled == 'TRUE' or self.scaled == 'True':
			self.method = 'frac2cart'
		else:
			self.method = 'cart2frac'
		self.coord = self.get_coords()
		self.uc = self.get_uc()

	def get_coords(self):
		parsing = False
		data =[]
		f = open(self.file, 'r')
		for line in f:
			line = line.strip()
			if line.startswith('&COORD'):
				parsing = True
			if line.startswith('&END COORD'):
				parsing = False
			if line.startswith('SCALED'):
				line.strip()
				self.scaled = line.split()[1]

			if line.startswith('UNIT'):
				line = line.strip()
				self.unit = line.split()[1]
				parsing = False

			if parsing and not line.startswith('&COORD'):
				data.append(line)
		coord = []
		for i, j in enumerate(data):
			j = j.split()
			element = j[0]
			x = float(j[1]) 
			y = float(j[2])
			z = float(j[3])
			coord.append([element, x, y, z])

		return coord

	def get_uc(self):
		cell={'a':0,'b':0,'c':0,'alpha':0,'beta':0,'gamma':0}
		parsing = False
		data = []
		f = open(self.file, 'r')
		for line in f:
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

		'''
		a= |a| =  sgr(d11^2 + d12^2 + d13^2)
		b = |b| = sqr( d21^2 + d22^2 + d23^2)
		c = |c| = sqr( d31^2 + d32^2 + d33^2)
		'''
		a = cell['a'] = math.sqrt(sum(matrix[0]**2))
		b = cell['b'] = math.sqrt(sum(matrix[1]**2))
		c = cell['c'] = math.sqrt(sum(matrix[2]**2))
		'''
		alpha = acos(  b.c/b*c )
		beta  = acos(  c.a/c*a )
		gamma = acos(  a.b/a*b )

		radians to degrees (x*180)/pi
		'''

		cell['alpha'] = ((math.acos(np.dot(matrix[1], matrix[2])/b*c))*180)/math.pi
		cell['beta'] = ((math.acos(np.dot(matrix[2], matrix[0])/c*a))*180)/math.pi
		cell['gamma'] = ((math.acos(np.dot(matrix[0], matrix[1])/a*b))*180)/math.pi

		return matrix, cell

	def det3(self, mat):
		return ((mat[0][0]*mat[1][1]*mat[2][2]) + (mat[0][1]*mat[1][2]*mat[2][0]) + (mat[0][2]*mat[1][0]*mat[2][1]) - (mat[0][2]*mat[1][1]*mat[2][0]) - (mat[0][1]*mat[1][0]*mat[2][2]) - (mat[0][0]*mat[1][2]*mat[2][1]))

	def frac2cart(self, cellParam, fracCoords):
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


	def cart2frac(self, cellParam, cartCoords):
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
		detLatCnt = self.det3(latCnt)
		for i in cartCoords:
		  aPos = (self.det3([[i[1], latCnt[0][1], latCnt[0][2]], [i[2], latCnt[1][1], latCnt[1][2]], [i[3], latCnt[2][1], latCnt[2][2]]])) / detLatCnt
		  bPos = (self.det3([[latCnt[0][0], i[1], latCnt[0][2]], [latCnt[1][0], i[2], latCnt[1][2]], [latCnt[2][0], i[3], latCnt[2][2]]])) / detLatCnt
		  cPos = (self.det3([[latCnt[0][0], latCnt[0][1], i[1]], [latCnt[1][0], latCnt[1][1], i[2]], [latCnt[2][0], latCnt[2][1], i[3]]])) / detLatCnt
		  fracCoords.append([i[0], aPos, bPos, cPos])

		return fracCoords

	def gen_cif(self):
		if self.method == 'cart2frac':
			coord = self.cart2frac(self.uc[0], self.coord)
		else:
			pass

		cif_head   = f"""\
data_CP2K
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                      {self.uc[1]['a']}
_cell_length_b                      {self.uc[1]['b']}
_cell_length_c                      {self.uc[1]['c']}
_cell_angle_alpha                   {self.uc[1]['alpha']}
_cell_angle_beta                    {self.uc[1]['beta']}
_cell_angle_gamma                   {self.uc[1]['gamma']}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z\n"""

		f = open('cp2k.cif', 'w')

		f.write(cif_head)
		kind = {}
		for i in range(len(coord)):
		 kind[coord[i][0]] = 0
		 
		for i in range(len(coord)):
			el = coord[i][0]
			kind[el] = kind[el]+1
			x = coord[i][1]
			y = coord[i][2]
			z = coord[i][3]
			cif_atoms = f"""{el+str(kind[el])}{' '*2}{el}{' '*2}{x}{' '*2}{y}{' '*2}{z}\n"""
			f.write(cif_atoms)
		f.close()

	def gen_xyz(self, coord):
		if self.method == 'frac2cart':
			self.frac2cart(self.uc[0], self.coord)
		else:
			pass
			'''
		 n_atoms = f"""{len(coord)}"""
		 xyz = f"""{x} {y} {z}"""
			'''

if __name__ == '__main__':
	fle_type = sys.argv[1]
	inpt = sys.argv[2]
	cp2cif = CP2KCIF(inpt)
	if fle_type == '-cif':
		cp2cif.gen_cif()
	elif fle_type == '-xyz':
		cp2cif.gen_xyz()
	else:
		print('options: -xyz <file> or -cif <file>')
	




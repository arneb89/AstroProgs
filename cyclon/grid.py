import numpy as np
from scipy.interpolate import interp1d

class abs_coeffs:
	mode = 0
	first_nu_rat = 0
	delta_nu_rat = 0
	number_nu_rat = 0
	tes = []
	thetas = []
	nu_rats = []
	coeffs_o = []
	coeffs_x = []
	
	def openfile(self, path):
		file = open(path, 'r')
		self.mode = int(file.readline())
		# nu_s
		pars = file.readline().split()
		self.first_nu_rat = float(pars[0])
		self.delta_nu_rat = float(pars[1])
		self.number_nu_rat = int(pars[2])
		self.nu_rats = np.zeros(self.number_nu_rat)
		for i in range(0, self.number_nu_rat):
			self.nu_rats[i] = self.first_nu_rat + i * self.delta_nu_rat
		# Te
		pars = file.readline().split()
		self.tes = np.zeros(int(pars[0]))
		for i in range(0, len(self.tes)):
			self.tes[i] = float(pars[i+1])
		# Thetas
		pars = file.readline().split()
		self.thetas = np.zeros(int(pars[0]))
		for i in range(0, len(self.thetas)):
			self.thetas[i] = float(pars[i+1])
		# Coeffs
		self.coeffs_o = np.zeros([len(self.tes), len(self.thetas), self.number_nu_rat])
		self.coeffs_x = np.zeros([len(self.tes), len(self.thetas), self.number_nu_rat])
		for i in range(0, len(self.tes)):
			for j in range(0, self.number_nu_rat):
				pars = file.readline().split()
				for k in range(0, len(self.thetas)):
					self.coeffs_o[i, k, j] = np.log(float(pars[k+1]))
				for k in range(0, len(self.thetas)):
					self.coeffs_x[i, k, j] = np.log(float(pars[k+len(self.thetas)+1]))
	
	def interp(self, mode, temp_num, theta, nu_rat):
		coeffs_thets = np.zeros(len(self.thetas))
		if(mode=='o'):
			for i in range(0, len(coeffs_thets)):
				y = self.coeffs_o[temp_num, i, :]
				f = interp1d(self.nu_rats, y, kind='linear', assume_sorted = True)
				coeffs_thets[i] = f(nu_rat)
		if(mode=='x'):
			for i in range(0, len(coeffs_thets)):
				y = self.coeffs_x[temp_num, i, :]
				f = interp1d(self.nu_rats, y, kind='linear', assume_sorted = True)
				coeffs_thets[i] = f(nu_rat)
		x = self.thetas
		y = coeffs_thets
		f = interp1d(x, y, kind='linear')
		res = f(theta)
		return res
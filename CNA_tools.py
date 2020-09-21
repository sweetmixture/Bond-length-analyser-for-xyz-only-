#!/bin/python

import sys
from math import sqrt

class CNA:

	def __init__(self,*IN):
	
		self.NOA=None					# Number of atoms in the input 'xyz' file
		self.input_xyz=IN[0]				# first argument, takes the main input 'xyz' file
		
		self.if_req_cutoff = False
		if len(IN) == 2:				# second argument, takes bond cutoff distance ... but currently does not have any functions
			self.cutoff = IN[1]
			self.if_req_cutoff = True
		else:
			self.cutoff = 999.

	def load_xyz(self):					# load 'xyz' file ... self.NOA : number of atoms
								# 		      self.config[i][j] : follows 'xyz' file convention
		try:
			with open(self.input_xyz,"r") as f:

				self.NOA = int(f.readline())
				rl	 = f.readline()
				self.config = [[0 for i in range(4)] for j in range(self.NOA)]
				for i in range(self.NOA):
					rl = f.readline()
					spl= rl.split()
					self.config[i][0] = spl[0]
					self.config[i][1] = float(spl[1])
					self.config[i][2] = float(spl[2])
					self.config[i][3] = float(spl[3])

		except FileNotFoundError:
			print("Error, input xyz not found ...")
			sys.exit()

	def get_atom_kind(self):
		elem = [ self.config[i][0] for i in range(self.NOA) ]					# SLICE THE 1st COLUMN ELEMENTS (WHICH SAVING ELEMENT NAMES)
		self.kind_list = []
		for i in elem:
			if i not in self.kind_list:							# IF NOT IN 'kind_list' THEN APPEND THE ELEMENT
				self.kind_list.append(i)
		self.NOE = len(self.kind_list)								# SAVE NUMBER OF UNIQUE ELEMENTS 

	def get_dist(self,A1,A2):				# get atom distance
		kind = [ A1[0], A2[0] ]
		kind.sort()
		dx = A1[1] - A2[1]
		dy = A1[2] - A2[2]
		dz = A1[3] - A2[3]
		return kind, sqrt(dx*dx + dy*dy + dz*dz)	# return 1. between which species 
								# 	 2. the distance

	def run_cna(self):

		self.cna_res = []	# Saves temporal CNA results :self.cna_res[i] = [ [kind , kind], distance ]
		for i in range(self.NOA):
			for j in range(i+1,self.NOA):
				combo, dist = self.get_dist(self.config[i],self.config[j])
				self.cna_res.append([combo,dist])
		
		combo = [ self.cna_res[i][0] for i in range(len(self.cna_res)) ] # Get unique combination of atoms in the system
		self.combo_list = []						 # e.g., self.combo_list = [ [a1,a1], [a1,a2], [a2,a2] ]
		for i in combo:
			if i not in self.combo_list:
				self.combo_list.append(i)

		self.final = [ [] for i in range( len(self.combo_list) ) ]# Get final result ...
									  # e.g., self.final[i] = [ [a_n,a_m] , r1, r2 ,r3 ,r4 ... ] (note that distance 'r' are not ordered)
		final_offset = 0
		for i in range(len(self.kind_list)):
			for j in range(i,len(self.kind_list)):
				X_X = None
				X_X = [self.kind_list[i],self.kind_list[j]]
				X_X.sort()
				self.final[final_offset].append(X_X)

				for k in range(len(self.cna_res)):
					if X_X == self.cna_res[k][0]:
						self.final[final_offset].append(self.cna_res[k][1])
				final_offset += 1
				

		### WRITE RESULT

		for i in range(len(self.combo_list)):

			outname = self.combo_list[i][0] + "_" + self.combo_list[i][1] + ".txt"
			tmp = self.final[i][1:]
			tmp.sort()
				
			with open(outname,"w") as w:

				for j in range(len(tmp)):
					w.write("%4d\t%12.6f\n" % (j+1,tmp[j]))

	def show(self):
		print(" %s%12s" % ("intput file name is :",self.input_xyz))
		print("\n config is ...\n")
		for i in range(self.NOA):
			print("%3s%12.6f%12.6f%12.6f" % (self.config[i][0],self.config[i][1],self.config[i][2],self.config[i][3]))

if __name__ == '__main__':

	if len(sys.argv) == 2:
		inst=CNA(sys.argv[1])
	elif len(sys.argv) == 3:
		inst=CNA(sys.argv[1],sys.argv[2])

	inst.load_xyz()
	inst.get_atom_kind()
	inst.run_cna()
	#inst.show()

# eigenstep.sage
# Copyright (C) 2015 Tim Haga <timhaga@ebene6.org>
# This file is released under the GPLv2

#Example from Fickus et al.
#x = 2/3
#y = 2/3
#lA = matrix([[0,0,0,0,0,0],[0,1,2-y,5/3,5/3,5/3],[0,0,y,4/3-x,5/3,5/3],[0,0,0,x,2/3,5/3]])

#Example from 2015-01-27

#lA = matrix([[0,0,0,0,0],[0,0,1/2,1,2],[0,1,5/2,5,8]])
#lA = matrix([[0,0,0,0,0,0],[0,1,3/2,2,5/2,5/2],[0,0,1/2,1,3/2,5/2]])

#Example from the whiteboard
#x = 1/3
#lA =  matrix([[0,0,0,0,0],[0,1,2-x,2,2],[0,0,x,1,2]])

#Determine N and d

def ConstructFrameFromEigensteps(lA):
	d = lA.dimensions()[0] - 1
	N = lA.dimensions()[1] - 1
	#Choose f_1 as first unit vector
	F = matrix(SR,d,N)
	F[0,0] = 1

	#For sake of simplicity choose U as identity. Every unitary matrix is allowed
	U = matrix.identity(d)

	for n in [1..N-1]:
		# step B.1
		# For sake of simplicity choose V as identity. Unitary block matrices with
		# blocksize according to multiplicity of the entrys of the respective column of lA
		V = matrix.identity(d)

		#step B.2
		I = [1..d]
		J = [1..d]

		for m in range(d,0,-1):
			deleter = false
			for m_1 in J:
				if lA[m][n] == lA[m_1][n+1]:
					deleter = true
			if deleter == true:
				del I[m-1]

				for m_2 in J:
					if lA[m_2][n+1] == lA[m][n]:
						m_3 = m_2
				del J[m_3-1]
		R = len(I)
		pi_I = Permutation(list(I) + list(Set([1..d]).difference(I))).inverse()
		Pi_I = pi_I.to_matrix()
		pi_J = Permutation(list(J) + list(Set([1..d]).difference(J))).inverse()
		Pi_J = pi_J.to_matrix()

		#step B.3
		#v is a size R vector but is filled up with padding zeros. So just set it to size d.
		v = matrix(SR,d,1)
		w = matrix(SR,R,1)
		for m in I:
			I_1 = list(Set(I).difference(Set([m])))
			v[pi_I[m-1]-1,0] = sqrt(-(prod((lA[m][n] - lA[m_2][n+1]) for m_2 in J))/(prod((lA[m][n] - lA[m_2][n]) for m_2 in I_1)))

		for m_1 in J:
			J_1 = list(Set(J).difference(Set([m_1])))
			w[pi_J[m_1-1]-1,0] = sqrt((prod((lA[m_1][n+1] - lA[m_2][n]) for m_2 in I))/(prod((lA[m_1][n+1] - lA[m_2][n+1]) for m_2 in J_1)))

		#step B.4
		#Compute the next frame vector
		F[:,n] = U*V*Pi_I.transpose()*v

		#step B.5
		W = matrix(SR,d,d)
		for m in I:
			for m_1 in J:
				W[pi_I[m-1]-1,pi_J[m_1-1]-1] = ((1/(lA[m_1][n+1] - lA[m][n]))*v[pi_I[m-1]-1,0]*w[pi_J[m_1-1]-1,0]).simplify_full()
		for i in range(R+1,d+1):
			W[i-1,i-1] = 1

		U = U*V*Pi_I.transpose()*W*Pi_J

	# end

	# do some simplification for easy-to-read output
	for i in range(d):
		for j in range(N):
			F[i,j] = F[i,j].simplify_radical()

	return F

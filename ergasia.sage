from itertools import permutations, combinations
from textwrap import wrap
import random

vec2pol = lambda L : sum([b*x^a for a,b in enumerate((L))])
str2vector = lambda s : vector(GF(2),[int(i) for i in s])
vec2bin = lambda c : ''.join([str(x) for x in c])


def FactorPolynomials(P):
	P = list( factor(P) )
	factors = [ x[0] for x in P ]
	return factors


def CyclicCodesList(P): # Calculate all available g(x) and h(x)
	P_origin_factors = FactorPolynomials(P) # [x + 1, x^3 + x + 1, x^3 + x^2 + 1] (all x^7 +1 factors)
	if len(P_origin_factors) > 1:
		tmp1 = list(combinations(P_origin_factors,2)) # [(x + 1, x^3 + x + 1), (x + 1, x^3 + x^2 + 1), (x^3 + x + 1, x^3 + x^2 + 1)]
		for i in range(len(tmp1)): # (x + 1, x^3 + x + 1)
			tmp1[i] = tmp1[i][0] * tmp1[i][1] # [(x + 1)*(x^3 + x + 1), (x + 1)*(x^3 + x^2 + 1), (x^3 + x + 1)*(x^3 + x^2 + 1)]

		tmp2 = []
		for i in tmp1: # Make all possible sets
			for y in P_origin_factors:
				if y not in FactorPolynomials(i):
					tmp2.append((i,y)) # (x + 1), (x^3 + x + 1)*(x^3 + x^2 + 1) | g(x), h(x)
					tmp2.append((y,i)) # (x^3 + x + 1)*(x^3 + x^2 + 1), (x + 1) | g(x), h(x)
		return tmp2

	else:
		raise ValueError(f'P(x) = {P} is not factorial.')


def g_h_chooser():
	g,h = random.choice(C_list)
	C_list.remove((g,h)) # To not be chosen again by C2
	return g,h
	

def Polynomial2Binary(P):
	return ''.join(str(x) for x in list(P))


def RandomSwapMatrixRows(M): # https://trac.sagemath.org/attachment/ticket/10426/trac_10426-matrix-row-column-swapping.patch
	tmpM = M
	while tmpM == M:
		for i in range(random.randint(M.nrows(), 2*M.nrows())):
			a = b = 0
			while a == b:
				a,b = [random.randint(0,M.nrows()-1), random.randint(0,M.nrows()-1)]
			M = M.with_swapped_rows(a,b)
	return M


def Str2Bin(msg):
	return "".join([bin(ord(c))[2:].zfill(8) for c in msg])# ['01010010', '00110100', '01000111', '01001110', '00110100', '01010010'] --> '010100100011010001000111010011100011010001010010'


padded = [False,False]
pad_cache = -1
def PadMsg(msg,n): # 010100100011010001000111010011100011010001010010 --> 010100100011010001000111010011100011010001010010 1011 00000000000000000000000000110000
	global padded
	global pad_cache
	if len(msg) > 256**4:
		raise ValueError(f'Too big msg.')

	pad_cache += 1
	if len(msg) % n != 0:
		padded[pad_cache] = True

		tmp = hex( len(msg) )[2:]
		suffix = bytes.fromhex( tmp.zfill( len(tmp) % 2 + len(tmp) ) ).rjust(4,bytes([0]))

		padded_msg = msg

		for i in range( ( n - ( (len(msg) + 4*8) % n) ) ):
			padded_msg += str(random.randint(0,1))

		padded_msg += "".join([bin(c)[2:].zfill(8) for c in suffix])
		return padded_msg # 010100100011010001000111010011100011010001010010101100000000000000000000000000110000

	else:
		return msg


def C1Encoder(msg):
	msg = wrap(msg, G_.nrows())
	c = ''
	for i in msg:
		c += vec2bin((str2vector(i)*G_))
	return c


def C2Encoder(msg):
	msg = wrap(msg, C2.generator_matrix().nrows())
	c = ''
	for i in msg:
		c += vec2bin( C2.encode( str2vector(i) ) )
	return c


def RandomNoise(msg):
	clear_msg_length = len(msg)
	if padded[1] == True:
		 clear_msg_length -= 4*8
		 
	encoded_msg = msg[:clear_msg_length]
	encoded_msg = wrap(encoded_msg, n)
	noise = ''

	C2_minimum_distance = C2.minimum_distance()

	for i in encoded_msg:
		while 1:
			a = str2vector(i)
			b = str2vector( ''.join([str(random.randint(0,1)) for i in range( len(a)) ]) )

			if vec2bin(a + b).count('1') < C2_minimum_distance:
				noise += vec2bin(b)
				break

	if padded[1] == True:
		noise += '0'*4*8

	return vec2bin( str2vector(c_) + str2vector(noise))


def C2Decoder(c):
	if padded[1] == True:
		s = c[:8*4]
		s = int(s,2)
		c_blocks = wrap(c[:s], n)
	else:
		c_blocks = wrap(c, n)
	c_decoded = ''
	for i in c_blocks:
		c_decoded += vec2bin( C2.decode_to_message(str2vector(i)) )

	return c_decoded


R.<x> = PolynomialRing(GF(2),"x")

n = int(input('n = ')) # length


P_origin = vec2pol( str2vector(input('Polynomial (example: 10000001 = x^7 + 1) = ') ) ) # x^7 + 1
C_list = CyclicCodesList(P_origin) # All possible Cs

# C1
g1,h1 = g_h_chooser()
C1 = codes.CyclicCode(generator_pol = g1, length = n)
G1 = codes.encoders.CyclicCodeVectorEncoder(C1)
H1 = C1.parity_check_matrix()

print(f'C1:\ng1(x) = {g1}\nh1(x) = {h1}\n\nG1 =\n{G1.generator_matrix()}\nH1 =\n{H1}\n\n')




# C2
g2,h2 = g_h_chooser()
C2 = codes.CyclicCode(generator_pol = g2, length = n)
G2 = codes.encoders.CyclicCodeVectorEncoder(C2)
H2 = C2.parity_check_matrix()

print(f'C2:\ng2(x) = {g2}\nh2(x) = {h2}\n\nG2 =\n{G2.generator_matrix()}\nH2 =\n{H2}\n\n')



# G' = S * G1 * P
S = random_matrix(GF(2),C1.generator_matrix().nrows(),C1.generator_matrix().nrows()) 
while not S.is_invertible():
	S = random_matrix(GF(2),C1.generator_matrix().nrows(),C1.generator_matrix().nrows())




G_ = S*G1.generator_matrix()

P = matrix.identity(G_.ncols())
# [1 0 0]
# [0 1 0]
# [0 0 1]

P = RandomSwapMatrixRows(P)
# [1 0 0]
# [0 0 1]
# [0 1 0]

G_ = G_*P

print(f"S =\n{S}\nP =\n{P}\nG' = S*G1*P =\n{G_}")





# MSG INPUT
msg = PadMsg( Str2Bin(input('msg (Example: Hello) = ')), G_.nrows() )
print(f'm = {msg}')




c = C1Encoder(msg)
# print(f'c = {c}')

c = PadMsg(c, C2.generator_matrix().nrows())

c_ = C2Encoder(c)
# print(f"c' (without noise) = {c_}")
c_ = RandomNoise(c_)
print(f"c' (with noise) = {c_}")

print(F"c' (decoded) = {C2Decoder(c_)}")
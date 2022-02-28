import math

###initalize class of elliptic curve by passing in arguments:  Elliptic_curve = Elliptic_curve(p, a, b, P, n, Q)
###BASIC POLLARD RHO: 
    ###call function Elliptic_curve.basic_prho()
    ###c, d, c', d' returned
    ###c, d, c_p, d_p = Elliptic_curve.basic_prho()


##use class of elliptic curve
class Elliptic_curve:

    def __init__(self, p, a, b, P, n, Q):
        ##point given
        self.P = P
        ##used to define curve
        self.a, self.b = a, b
        ##size of P
        self.n = n
        ##multiple of P
        self.Q = Q
        self.p = p
        ##origin 
        self.origin = None
        self.x_index = None
        ##all points on curve
        self.all_points = None

    ##following pseudocode from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    def eea(self, a, b):  
        old_r, r = a, b
        old_s, s = 1, 0
        old_t, t = 0, 1
            
        while r != 0:
            quotient = math.floor(old_r / r)
            (old_r, r) = (r, old_r - (quotient * r))
            (old_s, s) = (s, old_s - (quotient * s))
            (old_t, t) = (t, old_t - (quotient * t))
        
        if b != 0:
            t = math.floor((old_r - (old_s * a)) / b)
        else:
            t = 0
        
        return old_r, old_s, t

    ##inverts denominator of slope 
    def modular_inv(self, a, p):
        gcd, s, t = self.eea(a, p)
        if gcd != 1:
            raise ValueError
        return s % p

    ##negative modular of y value to check if P = -Q
    def neg_modular(self, y, p):
        y = p - (abs(y) % p)
        return y

    ##add points
    def point_add(self, P, Q, p, a):
        ##P, n, p, a
        x1, y1 = P[0], P[1]
        x2, y2 = Q[0], Q[1]

        ##if P or Q are the origin
        if P == self.origin:
            return Q
        elif Q == self.origin:
            return P
        ##3 conditions to consider
        # 1. x1 != x2
        # 2. x1 == x2 and y1 != y2
        # 3. x1 == x2 and y1 == y2   

        elif x1 == x2: 
            ##condition 2
            if y1 == self.neg_modular(y2, p):
                return self.origin
            ##condition 3
            elif y1 == y2:
                slope = ( ((3* pow(x1, 2)) + a) * self.modular_inv(2*y1, p) ) % p
                x3=(pow(slope, 2) - x1 - x2) % p
                y3=(slope*(x1-x3) - y1) % p
                return (x3, y3)
        ##condition 1
        else:
            slope = ((y2-y1) * self.modular_inv(x2-x1, p)) % p
            x3=(pow(slope, 2) - x1 - x2) % p
            y3=(slope*(x1-x3) - y1) % p
            return (x3, y3)

    ##double points
    def point_double(self, P, p, a):
        if P == self.origin:
            return self.origin
        x1, y1 = P[0], P[1]
        ##if  y1 = 0 then point double == inf
        if y1 == 0:
            return self.origin        
        slope = ((3 * pow(x1, 2) + a) * self.modular_inv(2*y1, p)) % p
        x3 = (pow(slope, 2) - (2*x1)) % p
        y3 = (((slope) * (x1-x3)) - y1) % p
        return(x3, y3)

    ##following: https://en.wikipedia.org/wiki/Euler's_criterion and Cryptogoraphy theory and practice: D. Stinson
    ##Euler criterion 
    def quadratic_res(self, n, p):
        if n % p == 0:
            return True
        else:
            power = (p-1) // 2
            z = pow(n, power, p)
            return z == 1
    
    #find all points of curve using tonelli-shanks algorithm, following pseudocode from: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    def tonelli_shank(self, p, n):
        ##y is 0
        if n % p == 0:
            return 0

        ##check if n has square roots
        if not self.quadratic_res(n, p):
            return None

        ##if n is quadratic residue and p == 3mod4
        if p == 3 % 4:
            power = (p + 1) // 4
            return pow(n, power, p)

        ##Find quadratic non-residue 
        if self.not_quadratic_res == None:
            z = 2
            while self.quadratic_res(z, p) == True:
                z += 1
            self.not_quadratic_res  = z
        else:
            z = self.not_quadratic_res 

        ##factor out powers of 2
        i = p - 1
        Q = 0
        S = 0            
        while i % 2 == 0:
            S += 1
            i //= 2        
        Q = i

        M = S
        c = pow(z, Q, p)
        t = pow(n, Q, p)
        power = (Q+1) // 2
        R = pow(n, power, p)

        ##repeated squaring find i
        while t != 1:
            i = 0
            t_p = t 
            if t == 0:
                return 0
            while t_p != 1:
                i += 1
                t_p = pow(t_p, 2, p)

            power = pow(2, (M-i-1))
            b = pow(c, power, p)
            M = i
            c = pow(b, 2, p)
            t = (t * pow(b, 2)) % p
            R = (R*b) % p
        return R
	
    ##if  valid point
    def y_is_point(self, y, x, a, b, p):
        return ( pow(y, 2) - ( pow(x, 3) + (a*x) + b) ) % p == 0

    ##get all the points on the curve
    def get_curve(self, p, a, b):
        all_points = []
        x = 0
        self.not_quadratic_res = None
        i = 0
        
        ##add intial point to list of points
        while len(all_points) < 1:
            n = ( pow(x, 3) + (a * x) + b) % p
            y =  self.tonelli_shank(p, n)
            if y == None:
                x += 1
            elif self.y_is_point(y, x, a, b, p):
                all_points.append((x, y))
                i += 1
                if y != 0:
                    all_points.append((x, self.neg_modular(-y, p)))
                    i += 1
                x += 1

        ##while beginning of list not reached(found all points)
        while (x % p != all_points[0][0]):
            n = (pow(x, 3) + a * x + b) % p
            y =  self.tonelli_shank(p, n)
            if y == None:
                x += 1
            elif self.y_is_point(y, x, a, b, p):
                all_points.append((x, y))
                i += 1
                ##if y not zero, get other y
                if y != 0:
                    all_points.append((x, self.neg_modular(-y, p)))
                    i += 1
                x += 1

        ##add origin to points
        all_points.append( (0, 0) )
        ##save origin
        self.origin = (0, 0)
        self.all_points = all_points
        return all_points

    ##turn number to binary
    def dec_to_bin(self, n):
        b = list()
        while n > 0:
            b.append(n % 2)
            n = n // 2
        return b[::-1]

    ##scalar multiplication - when point multiplied by scalar
    ##double and add algorithm
    def scalar_multiplication(self, n, P, p, a):
        ##get binary of number, n
        binary_n = self.dec_to_bin(n)

        n0 = P 
        ##double and add
        for i in range(0, len(binary_n)):
            if i == 0:
                continue
            ##if 0 only double the point
            elif binary_n[i] == 0:
                n0 = self.point_double(n0, p, a)
            ##if 1, double the point and then add 
            elif binary_n[i] == 1:
                n0 = self.point_double(n0, p, a)
                n0 = self.point_add(P, n0, p, a)
        return n0

    ##iterative function for pollard rho, changes value of x and s and/or t depending on whether x in s1, s2 or s3
    def iterative_f(self, S, x, s, t):
        if S == 1:
            x = self.point_add(self.Q, x, self.p, self.a)
            t = (t + 1) % self.p
            return x, s, t
        elif S == 2:
            x = self.point_double(x, self.p, self.a)
            s = (2 * s) % self.p
            t = (2 * t) % self.p
            return x, s, t
        elif S == 3:
            x = self.point_add(self.P, x, self.p, self.a)
            s = (s + 1) % self.p
            return x, s, t
        else:
            return "error"
    
    ##find if point in S1, S2 or S3
    def find_S(self, x, S1, S2, S3):
        if x in S1:
            return 1
        elif x in S2:
            return 2
        elif x in S3:
            return 3
        else:
            print("ERROR")

    ##basic pollard rho algorithm
    def basic_prho(self, Q = None): 
        ##find all points on curve
        if self.all_points == None:
            self.all_points = self.get_curve(self.p, self.a, self.b)
        

        ##split points into three lists
        L = len(self.all_points) // 3
        L_p = L+L
        
        S1 = self.all_points[:L]
        S2 = self.all_points[L:L_p]
        S3 = self.all_points[L_p:]
        
        ##if no starting point given, start from 0
        if self.x_index == None:
            self.x_index = 1

        ##initlalize s and t to 0
        s, t = 0, 0

        ##find x
        x = self.all_points[self.x_index]
        ##find where x is located
        S = self.find_S(x, S1, S2, S3)
        
        ##if x starting point in s2
        while S == 2:
            self.x_index = (self.x_index + 1) % self.n
            x = self.all_points[self.x_index]
            S = self.find_S(x, S1, S2, S3)

        ##find x', s', t'
        x_p, s_p, t_p = self.iterative_f(S, x, s, t)

        ##loop through until find collision
        while x != x_p:
            ##find where x
            S = self.find_S(x, S1, S2, S3)
            ##change value, x and s and/or t
            x, s, t = self.iterative_f(S, x, s, t)
            ##triple property x = g**s A**t
            x = self.point_add(self.scalar_multiplication(s, self.P, self.p, self.a), self.scalar_multiplication(t, self.Q, self.p, self.a), self.p, self.a)
            
            ##big jump - do twice
            ##find where x'
            S_p = self.find_S(x_p, S1, S2, S3)
            ##change value of x and s and/or t
            x_p, s_p ,t_p = self.iterative_f(S_p, x_p, s_p, t_p)
            ##find where x'
            S_p = self.find_S(x_p, S1, S2, S3)
            ##change value of x and s and/or t
            x_p, s_p ,t_p = self.iterative_f(S_p, x_p, s_p, t_p)
            ##triple property x = g**s A**t
            x_p = self.point_add(self.scalar_multiplication(s_p, self.P, self.p, self.a), self.scalar_multiplication(t_p, self.Q, self.p, self.a), self.p, self.a)
            
        return s, t, s_p, t_p
   

#Input:
p = 16001
a = 10
b = 1
P = (1654, 7208)
n = 8026
Q = (5000, 1283)

Elliptic_curve = Elliptic_curve(p, a, b, P, n, Q)

c, d, c_p, d_p = Elliptic_curve.basic_prho()
print("collision:")
print("c:", c)
print("d:", d)
print("c':", c_p)
print("d':", d_p)

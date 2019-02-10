import numpy as np
import math
import matplotlib.pyplot as plt
pi = math.pi

class Complex:
    def __init__(self, real, imag):
        self.re = real
        self.im = imag
        
    def __str__(self):
        if self.im == 0:
            return str(self.re)
        elif self.im < 0:
            return str(self.re)+' - '+str(-self.im)+' i'
        else:
            return str(self.re)+' + '+str(self.im)+' i'
    
    def __eq__(self, other):
        return self.re == other.re and self.im == other.im

    def __hash__(self):
        return hash((self.re, self.im))
    
    def abs(self):
        return math.hypot(self.re, self.im)
    
    def phase(self):
        return math.atan2(im, re)
        
    def plus(self, c):
        real = c.re + self.re
        imag = c.im + self.im
        return Complex(real, imag)
    
    def minus(self, c):
        real = self.re - c.re
        imag = self.im - c.im
        return Complex(real, imag)
    
    def times(self, c):
        real = self.re * c.re - self.im * c.im
        imag = self.re * c.im + self.im * c.re
        return Complex(real, imag)
    
    def scale(self, alpha):
        return Complex(self.re*alpha, self.im*alpha)
    
    def conjugate(self):
        return Complex(self.re, -self.im)
    
    def reciprocal(self):
        scale = self.re*self.re + self.im*self.im
        return Complex(self.re/scale, -self.im/scale)
    
    def real(self):
        return self.re
    
    def imag(self):
        return self.im
    
    def divides(self, c):
        return self.times(c.reciprocal())
    
    def exp(self):
        re = self.re
        im = self.im
        return Complex(math.exp(re) * math.cos(im), math.exp(re) * math.sin(im))
    
    def sin(self):
        re = self.re
        im = self.im
        return Complex(math.sin(re) * math.cosh(im), math.cos(re) * math.sinh(im))
    
    def cos(self):
        re = self.re
        im = self.im
        return Complex(math.cos(re) * math.cosh(im), -math.sin(re) * math.sinh(im))
    
    def tan(self):
        return self.sin().divides(self.cos())

def fft(x):        
    if x[0] is None:
        print "error"
        return None
    n = len(x)
    if n == 1:
        return [x[0]]
    if n % 2 != 0:
        print "n is not a power of 2"
        return None
    even = [None] * (n/2)
    for k in range(0, n/2):
        even[k] = x[2*k]
    q = fft(even)
    
    odd = even # reuse the array
    for k in range(0, n/2):
        odd[k] = x[2*k + 1]
    r = fft(odd)
    
    # combine even and odd
    y = [None] * n
    for k in range(0, n/2):
        kth = -2*k*pi/n
        wk = Complex(math.cos(kth), math.sin(kth))
        y[k]       = q[k].plus(wk.times(r[k]))
        y[k + n/2] = q[k].minus(wk.times(r[k]))
    return y

# compute the inverse FFT
def ifft(x):
    n = len(x)
    
    y = [None] * n
    # take conjugate
    for i in range(n):
        y[i] = x[i].conjugate()
    
    y = fft(y)
    
    # take conjugate again
    for i in range(n):
        y[i] = y[i].conjugate()
        
    for i in range(n):
        y[i] = y[i].scale(1.0 / n)
        
    return y

# compute the circular convolution of x and y
def cconvolve(x, y):
    if len(x) != len(y):
        print "error, x y length does not match"
        return None
    n = len(x)
    # compute FFT of each sequence
    a = fft(x)
    b = fft(y)
    
    # point-wise multiply
    c = [None] * n
    for i in range(n):
        c[i] = a[i].times(b[i])
    
    # compute inverse FFT
    return ifft(c)

# linear convolution of x and y
def convolve(x, y):
    zero = Complex(0, 0)
    a = [None] * (2*len(x))
    for i in range(len(x)):
        a[i] = x[i]
    for i in range(len(x), 2*len(x)):
        a[i] = zero
        
    b = [None] * (2*len(y))
    for i in range(len(y)):
        b[i] = y[i]
    for i in range(len(y), 2*len(y)):
        b[i] = zero
    
    return cconvolve(a, b)

def show(x):
    for item in x:
        print item

def hammingWindow(n):
    pi = math.pi
    w = [None] * n
    for i in range(n):
        w[i] = 0.54 - 0.46*math.cos(2*pi*i/(n-1))
    return w

file_names = ["p1-1_Female_20-29_170-179cm_Hand_held.txt",
"p1-1_Female_20-29_170-179cm_Trousers_back_pocket.txt",
"p1-2_Female_20-29_170-179cm_Handheld_using.txt",
"p1-3_Female_20-29_170-179cm_Trousers_back_pocket.txt",
"p1-4_Female_20-29_170-179cm_Handbag.txt",
"p10-1_Male_20-29_170-179cm_Hand_held.txt",
"p10-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p10-2_Male_20-29_170-179cm_Handheld_using.txt",
"p10-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p10-3_Male_20-29_170-179cm_Backpack.txt",
"p11-1_Male_20-29_170-179cm_Hand_held.txt",
"p11-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p11-2_Male_20-29_170-179cm_Handheld_using.txt",
"p11-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p11-3_Male_20-29_170-179cm_Backpack.txt",
"p12-2_Male_15-19_180-189cm_Hand_held.txt",
"p12-2_Male_15-19_180-189cm_Trousers_front_pocket.txt",
"p12-3_Male_15-19_180-189cm_Handheld_using.txt",
"p12-3_Male_15-19_180-189cm_Trousers_back_pocket.txt",
"p12-4_Male_15-19_180-189cm_Backpack.txt",
"p13-1_Female_20-29_160-169cm_Hand_held.txt",
"p13-2_Female_20-29_160-169cm_Handheld_using.txt",
"p13-2_Female_20-29_160-169cm_Trousers_front_pocket.txt",
"p13-3_Female_20-29_160-169cm_Trousers_back_pocket.txt",
"p13-4_Female_20-29_160-169cm_Backpack.txt",
"p14-1_Male_20-29_160-169cm_Hand_held.txt",
"p14-1_Male_20-29_160-169cm_Trousers_front_pocket.txt",
"p14-2_Male_20-29_160-169cm_Handheld_using.txt",
"p14-2_Male_20-29_160-169cm_Trousers_back_pocket.txt",
"p14-3_Male_20-29_160-169cm_Backpack.txt",
"p15-1_Male_20-29_180-189cm_Hand_held.txt",
"p15-1_Male_20-29_180-189cm_Trousers_front_pocket.txt",
"p15-2_Male_20-29_180-189cm_Handheld_using.txt",
"p15-2_Male_20-29_180-189cm_Trousers_back_pocket.txt",
"p15-3_Male_20-29_180-189cm_Backpack.txt",
"p16-1_Male_15-19_180-189cm_Hand_held.txt",
"p16-1_Male_15-19_180-189cm_Trousers_front_pocket.txt",
"p16-2_Male_15-19_180-189cm_Handheld_using.txt",
"p16-2_Male_15-19_180-189cm_Trousers_back_pocket.txt",
"p16-3_Male_15-19_180-189cm_Backpack.txt",
"p17-1_Female_20-29_150-159cm_Hand_held.txt",
"p17-1_Female_20-29_150-159cm_Trousers_front_pocket.txt",
"p17-2_Female_20-29_150-159cm_Handheld_using.txt",
"p17-2_Female_20-29_150-159cm_Trousers_back_pocket.txt",
"p17-3_Female_20-29_150-159cm_Handbag.txt",
"p18-1_Male_20-29_180-189cm_Hand_held.txt",
"p18-1_Male_20-29_180-189cm_Trousers_front_pocket.txt",
"p18-2_Male_20-29_180-189cm_Handheld_using.txt",
"p18-2_Male_20-29_180-189cm_Trousers_back_pocket.txt",
"p18-3_Male_20-29_180-189cm_Backpack.txt",
"p19-1_Female_20-29_150-159cm_Hand_held.txt",
"p19-2_Female_20-29_150-159cm_Handheld_using.txt",
"p19-2_Female_20-29_150-159cm_Trousers_back_pocket.txt",
"p19-3_Female_20-29_150-159cm_Handbag.txt",
"p2-1_Male_20-29_180-189cm_Trousers_front_pocket.txt",
"p2-2_Male_20-29_180-189cm_Hand_held.txt",
"p2-2_Male_20-29_180-189cm_Trousers_back_pocket.txt",
"p2-3_Male_20-29_180-189cm_Handheld_using.txt",
"p2-3_Male_20-29_180-189cm_Shirt_pocket.txt",
"p2-4_Male_20-29_180-189cm_Backpack.txt",
"p20-1_Male_20-29_170-179cm_Hand_held.txt",
"p20-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p20-2_Male_20-29_170-179cm_Handheld_using.txt",
"p20-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p20-3_Male_20-29_170-179cm_Backpack.txt",
"p21-1_Male_20-29_180-189cm_Hand_held.txt",
"p21-1_Male_20-29_180-189cm_Trousers_front_pocket.txt",
"p21-2_Male_20-29_180-189cm_Handheld_using.txt",
"p21-2_Male_20-29_180-189cm_Trousers_back_pocket.txt",
"p21-3_Male_20-29_180-189cm_Backpack.txt",
"p22-1_Male_20-29_170-179cm_Hand_held.txt",
"p22-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p22-2_Male_20-29_170-179cm_Handheld_using.txt",
"p22-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p22-3_Male_20-29_170-179cm_Backpack.txt",
"p23-1_Female_20-29_160-169cm_Hand_held.txt",
"p23-1_Female_20-29_160-169cm_Trousers_front_pocket.txt",
"p23-2_Female_20-29_160-169cm_Handheld_using.txt",
"p23-2_Female_20-29_160-169cm_Trousers_back_pocket.txt",
"p23-3_Female_20-29_160-169cm_Handbag.txt",
"p24-1_Male_20-29_170-179cm_Hand_held.txt",
"p24-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p24-2_Male_20-29_170-179cm_Handheld_using.txt",
"p24-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p24-3_Male_20-29_170-179cm_Backpack.txt",
"p25-1_Female_20-29_170-179cm_Hand_held.txt",
"p25-3_Female_20-29_170-179cm_Handheld_using.txt",
"p25-4_Female_20-29_170-179cm_Handbag.txt",
"p26-1_Female_20-29_150-159cm_Backpack.txt",
"p26-2_Female_20-29_150-159cm_Hand_held.txt",
"p26-2_Female_20-29_150-159cm_Trousers_front_pocket.txt",
"p26-3_Female_20-29_150-159cm_Handheld_using.txt",
"p27-1_Male_15-19_170-179cm_Hand_held.txt",
"p27-1_Male_15-19_170-179cm_Trousers_front_pocket.txt",
"p27-2_Male_15-19_170-179cm_Trousers_back_pocket.txt",
"p27-4_Male_15-19_170-179cm_Handheld_using.txt",
"p27-5_Male_15-19_170-179cm_Backpack.txt",
"p3-1_Male_20-29_170-179cm_Hand_held.txt",
"p3-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p3-2_Male_20-29_170-179cm_Handheld_using.txt",
"p3-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p3-3_Male_20-29_170-179cm_Shirt_pocket.txt",
"p3-4_Male_20-29_170-179cm_Backpack.txt",
"p4-1_Male_20-29_170-179cm_Hand_held.txt",
"p4-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p4-2_Male_20-29_170-179cm_Handheld_using.txt",
"p4-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p4-3_Male_20-29_170-179cm_Backpack.txt",
"p5-1_Female_20-29_160-169cm_Hand_held.txt",
"p5-1_Female_20-29_160-169cm_Trousers_front_pocket.txt",
"p5-2_Female_20-29_160-169cm_Handheld_using.txt",
"p5-3_Female_20-29_160-169cm_Handbag.txt",
"p6-1_Male_15-19_180-189cm_Hand_held.txt",
"p6-1_Male_15-19_180-189cm_Trousers_front_pocket.txt",
"p6-2_Male_15-19_180-189cm_Handheld_using.txt",
"p6-2_Male_15-19_180-189cm_Trousers_back_pocket.txt",
"p6-3_Male_15-19_180-189cm_Backpack.txt",
"p7-1_Male_20-29_180-189cm_Hand_held.txt",
"p7-1_Male_20-29_180-189cm_Trousers_front_pocket.txt",
"p7-2_Male_20-29_180-189cm_Handheld_using.txt",
"p7-2_Male_20-29_180-189cm_Trousers_back_pocket.txt",
"p7-3_Male_20-29_180-189cm_Backpack.txt",
"p8-1_Male_20-29_170-179cm_Hand_held.txt",
"p8-1_Male_20-29_170-179cm_Trousers_front_pocket.txt",
"p8-2_Male_20-29_170-179cm_Handheld_using.txt",
"p8-2_Male_20-29_170-179cm_Trousers_back_pocket.txt",
"p8-3_Male_20-29_170-179cm_Backpack.txt",
"p9-1_Female_15-19_160-169cm_Hand_held.txt",
"p9-1_Female_15-19_160-169cm_Trousers_front_pocket.txt",
"p9-2_Female_15-19_160-169cm_Handheld_using.txt",
"p9-2_Female_15-19_160-169cm_Trousers_back_pocket.txt",
"p9-3_Female_15-19_160-169cm_Backpack.txt"]

for i in range(len(file_names)):
    file_names[i] = "dataset/"+file_names[i]

for file_name in file_names:
    accMagn = []
    fr = open(file_name, 'r')
    for line in fr.readlines():
        items = line.split(';')
        x = float(items[1])
        y = float(items[2])
        z = float(items[3])
        magn = 1/3.0 * math.sqrt(x*x + y*y + z*z)
        accMagn.append(magn)
    fr.close()
    print file_name.split('/')[1] + " *read over", len(accMagn)

    # n denotes the length of the DFT window, 400 tuples -> 4 second, 
    # according to the ubicomp paper
    fs = 100 # 100 HZ sampling
    n = 400 
    start = 0
    sliding = 30
    low_fs = 0.6
    high_fs = 2.5
    energy_thred = 1.0
    walking = False
    while(start + n <= len(accMagn)):
        l = 1024
        w = hammingWindow(n)
        # zero padding to 1024
        # for FFT
        fft_input = [Complex(0.0, 0.0)] * l
        for i in range(n):
            val = accMagn[start+i] * w[i]
            fft_input[i] = Complex(val, 0)

        fft_res = fft(fft_input)
        fft_output = [0.0] * l
        for i in range(l):
            fft_output[i] = fft_res[i].abs()

        energy = 0.0
        for i in range(int(low_fs/fs*l), int(high_fs/fs*l)+1):
            energy += fft_output[i] * fft_output[i] * 1.0 / l

        if energy > energy_thred:
            if walking == False:
                print start, "start walking"
                walking = True
        else:
            if walking == True:
                print start, "stop walking"
                walking = False

        start = start + sliding
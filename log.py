from __future__ import division
import numpy as np
from pylab import plot, show, xlabel, ylabel
import random
import web
from Bio import motifs
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import Axes3D
#from bs4 import BeautifulSoup
import time

def gradientDescent(x, y, theta, alpha, m, numIterations):
    J_history = np.zeros(shape=(numIterations, 1))
    xTrans = x.transpose()
    for i in range(0, numIterations):
        hypothesis = np.dot(x, theta)

        loss = hypothesis - y
        # avg cost per example (the 2 in 2*m doesn't really matter here.
        # But to be consistent with the gradient, I include it)
        cost = np.sum(loss ** 2) / (2 * m)
        #print("Iteration %d | Cost: %f" % (i, cost))
        # avg gradient per example
        gradient = np.dot(xTrans, loss) / m
        # update
        theta = theta - alpha * gradient
        J_history[i][0] = cost
    return theta, J_history

def compute_cost(x, y,m, theta):
    '''
    Comput cost for linear regression
    '''
    #Number of training samples
    #m = y.size
    predictions = np.dot(x, theta)
    sqErrors = (predictions - y)
    J = (1.0 / (2 * m)) * sqErrors.T.np.dot(sqErrors)
    return J
def gradient_descent(x, y, theta, alpha, numIterations):
    '''
    Performs gradient descent to learn theta
    by taking num_items gradient steps with learning
    rate alpha
    '''
    #m = y.size
    J_history = np.zeros(shape=(numIterations, 1))
    for i in range(numIterations):
        predictions = np.dot(x, theta)
        theta_size = theta.size
        for it in range(theta_size):
            temp = x[:, it]
            temp.shape = (m, 1)
            errors_x1 = (predictions - y) * temp
            theta[it][0] = theta[it][0] - alpha * (1.0 / m) * errors_x1.sum()
        J_history[i, 0] = compute_cost(x, y, theta)
    return theta, J_history



print "\t\t\t -35"


instances = [Seq("TTGACG"),
Seq("TTTACA"),
Seq("TTGACA"),
Seq("CTGATA"),
Seq("TTGACA"),
Seq("TTTACG"),
Seq("TTTACG"),
Seq("TTTACG"),
Seq("CTGACA"),
Seq("TTTACA"),
Seq("TTTACG"),
Seq("TTGACG"),
Seq("CTGATA"),
Seq("CTGATG"),
Seq("TTTATG"),
Seq("TTTATA"),
Seq("TTGACA"),
Seq("TTGACA"),
Seq("TTGACG"),
]
""
m = motifs.create(instances)
print(m)
print(m.counts);
pwm = m.counts.normalize(pseudocounts= {'A':0.49, 'C': 0.51, 'G': 0.51, 'T': 0.49})
print(pwm)
pssm = pwm.log_odds()
print(pssm)
#return pssm
result = [[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0]
] 

def calculateX(a,b,c,d,e,f,x):
    result[x][0] = pssm[a,0] + pssm[b,1] + pssm[c,2] + pssm[d,3] + pssm[e,4] + pssm[f,5]
    #return result[x][0] 
calculateX('T', 'T', 'G', 'A', 'C', 'G',0)    
calculateX('T', 'T', 'T', 'A', 'C', 'A',1)
calculateX('T', 'T', 'G', 'A', 'C', 'A',2)    
calculateX('C', 'T', 'G', 'A', 'T', 'A',3)
calculateX('T', 'T', 'G', 'A', 'C', 'A',4)    
calculateX('T', 'T', 'T', 'A', 'C', 'G',5)
calculateX('T', 'T', 'T', 'A', 'C', 'G',6)    
calculateX('T', 'T', 'T', 'A', 'C', 'G',7)
calculateX('C', 'T', 'G', 'A', 'C', 'A',8)    
calculateX('T', 'T', 'T', 'A', 'C', 'A',9)
calculateX('T', 'T', 'T', 'A', 'C', 'G',10)
calculateX('T', 'T', 'G', 'A', 'C', 'G',11)
calculateX('C', 'T', 'G', 'A', 'T', 'A',12)    
calculateX('C', 'T', 'G', 'A', 'T', 'A',13)
calculateX('T', 'T', 'T', 'A', 'T', 'G',14)    
calculateX('T', 'T', 'T', 'A', 'T', 'A',15)
calculateX('T', 'T', 'G', 'A', 'C', 'A',16)    
calculateX('T', 'T', 'G', 'A', 'C', 'A',17)
calculateX('T', 'T', 'G', 'A', 'C', 'G',18)   
outputResult = np.log([[[1.0],
[0.7],
[0.86],
[0.01],
[0.72],
[0.24],
[0.47],
[0.36],
[0.51],
[0.04],
[0.33],
[0.58],
[0.01],
[0.01],
[0.1],
[0.15],
[0.16],
[0.06],
[0.56]
]])
print ""
print "X2 Values"
print result 
print ""
print type(result)
print ""
print len(result)
print ""
print "Y Values"
print outputResult


print ""
print "\t\t\t\t -10"

instances2 = [Seq("TACAGT"),
Seq("TATTAT"),
Seq("TACTGT"),
Seq("GATTAT"),
Seq("TATTGT"),
Seq("TACTAT"),
Seq("TATAGT"),
Seq("TATTAT"),
Seq("TATAAT"),
Seq("GACTGT"),
Seq("TACAAT"),
Seq("TATAGT"),
Seq("GATTAT"),
Seq("GATTAT"),
Seq("TACAAT"),
Seq("TACAAT"),
Seq("GACTAT"),
Seq("GATTGT"),
Seq("TATTGT")
]


m2 = motifs.create(instances2)
print(m2.counts);
pwm2 = m2.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G': 0.51, 'T': 0.49})
print(pwm2)
pssm2 = pwm2.log_odds()
print(pssm2)
        
result2 = [[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[0]
] 


def calculateX2(a,b,c,d,e,f,x):
    result2[x][0] = pssm2[a,0] + pssm2[b,1] + pssm2[c,2] + pssm2[d,3] + pssm2[e,4] + pssm2[f,5]
    #return result[x][0] 
calculateX2('T', 'A', 'C', 'A', 'G', 'T',0)    
calculateX2('T', 'A', 'T', 'T', 'A', 'T',1)
calculateX2('T', 'A', 'C', 'T', 'G', 'T',2)    
calculateX2('G', 'A', 'T', 'T', 'A', 'T',3)
calculateX2('T', 'A', 'T', 'T', 'G', 'T',4)    
calculateX2('T', 'A', 'C', 'T', 'A', 'T',5)
calculateX2('T', 'A', 'T', 'A', 'G', 'T',6)    
calculateX2('T', 'A', 'T', 'T', 'A', 'T',7)
calculateX2('T', 'A', 'T', 'A', 'A', 'T',8)    
calculateX2('G', 'A', 'C', 'T', 'G', 'T',9)
calculateX2('T', 'A', 'C', 'A', 'A', 'T',10)
calculateX2('T', 'A', 'T', 'A', 'G', 'T',11)
calculateX2('G', 'A', 'T', 'T', 'A', 'T',12)    
calculateX2('G', 'A', 'T', 'T', 'A', 'T',13)
calculateX2('T', 'A', 'C', 'A', 'A', 'T',14)    
calculateX2('T', 'A', 'C', 'A', 'A', 'T',15)
calculateX2('G', 'A', 'C', 'T', 'A', 'T',16)    
calculateX2('G', 'A', 'T', 'T', 'G', 'T',17)
calculateX2('T', 'A', 'T', 'T', 'G', 'T',18)   
outputResult2 = np.log([[1],
[0.7],
[0.86],
[0.01],
[0.72],
[0.24],
[0.47],
[0.36],
[0.51],
[0.04],
[0.33],
[0.58],
[0.01],
[0.01],
[0.1],
[0.15],
[0.16],
[0.06],
[0.56]
])
print ""
print "X2 Values"
print result2 
print ""
print type(result2)
print ""
print "Y Values"
print outputResult2



a = [
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
[1,0,0],
]

print len(a)

"""
for i in range(0, 6):
    print result[i][0]
    a[i][1] = result[i][0]
    i+=1
    print a[i][1]
"""
print ""
a[0][1] = result[0][0]
a[1][1] = result[1][0]
a[2][1] = result[2][0]
a[3][1] = result[3][0]
a[4][1] = result[4][0]
a[5][1] = result[5][0]
a[6][1] = result[6][0]
a[7][1] = result[7][0]
a[8][1] = result[8][0]
a[9][1] = result[9][0]
a[10][1] = result[10][0]
a[11][1] = result[11][0]
a[12][1] = result[12][0]
a[13][1] = result[13][0]
a[14][1] = result[14][0]
a[15][1] = result[15][0]
a[16][1] = result[16][0]
a[17][1] = result[17][0]
a[18][1] = result[18][0]

a[0][2] = result2[0][0]
a[1][2] = result2[1][0]
a[2][2] = result2[2][0]
a[3][2] = result2[3][0]
a[4][2] = result2[4][0]
a[5][2] = result2[5][0]
a[6][2] = result2[6][0]
a[7][2] = result2[7][0]
a[8][2] = result2[8][0]
a[9][2] = result2[9][0]
a[10][2] = result2[10][0]
a[11][2] = result2[11][0]
a[12][2] = result2[12][0]
a[13][2] = result2[13][0]
a[14][2] = result2[14][0]
a[15][2] = result2[15][0]
a[16][2] = result2[16][0]
a[17][2] = result2[17][0]
a[18][2] = result2[18][0]

print ""
print "Matrix A (Input Matrix)"

for x in a:
    print x
    print ""




#newchanges
b= np.log([ 1,
0.7,
0.86,
0.01,
0.72,
0.24,
0.47,
0.36,
0.51,
0.04,
0.33,
0.58,
0.01,
0.01,
0.1,
0.15,
0.16,
0.06,
0.56
])


print ""
print "Matrix A (Input Matrix)"

for x in b:
    print x
    print ""







#a = [[1,2,3],[1,4,5],[1,5,6]]
#b = [1,2,3]
#print a
print "X Matrix"
x = np.asarray(a)
print x
print ""
print "Y Matrix"
y = np.asarray(b)
print y
print ""
#print type(x1)
#print type(y1)
#print x1
#print y1
#print x1
#print type(x1)
"""x,y = scale(x), y
print(x)
print(y)
"""
m, n = np.shape(x)
numIterations= 200000
alpha = 0.015
theta = np.ones(n)

theta, J_history = gradientDescent(x, y, theta, alpha,m,numIterations)
print "theta"
print(theta)
print ""
print "hx"
hx = x.dot(theta)
print hx
#hx_flatten = hx.flatten()
print ""
print "Difference"
diff = hx - y
print diff
print ""
print "Difference Square"
diff_square = diff*diff
print diff_square 
print ""
print "Sum"
sum = np.sum(diff_square) 
print sum
print ""
print "Cost function"
#cost = 1/(2*m)*sum(diff_square)
temp = 1/(2*m)
cost = temp*sum
print cost
#X_norm, mean_r, std_r = feature_normalize(x)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for c, m in [('r', 'o')]:
    xs = x[:,1]
    #xs = x[i,1]
    print "xs"
    print (xs)
    ys = x[:,2]
    #ys = x[i,2]
    print "ys"
    print(ys)
    zs = y
    print"zs"
    print(zs)
    ax.scatter(xs, ys, zs, c=c, marker=m)

ax.set_xlabel('-10 hexamer')
ax.set_ylabel('-35 hexamer')
ax.set_zlabel('strength of promoter')
plt.show()
#theta = zeros(shape=(3, 1))

print theta, J_history
plot(np.arange(numIterations), J_history)
xlabel('Iterations')
ylabel('Cost Function')
show()

print theta


strength = np.array([1.0,  8.50, 8.65 ]).dot(theta)
l = np.exp(strength)   #newchanges
print 'Predicted strength of promoter : %f' % (l)

meany = np.mean(y)
print meany
print y
sumsqmeany = np.sum((y - meany)**2)
print sumsqmeany

sumsqmeanysum = np.sum((y - hx)**2)/sumsqmeany

R = 1 - sumsqmeanysum
print "The R value is:"
print R 

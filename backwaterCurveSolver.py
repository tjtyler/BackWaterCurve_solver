#---------- FUNCTION DEFINITIONS --------------
def determineNormalDepths(b,n,S0,EPS,Q,c):
  y_try = Q/((b**2)*(1/5))
  Q_t = 0
  deltaQ = Q - Q_t
  while (abs(deltaQ) > (EPS*Q)):
    Q_t = (c/n)*(b*y_try)*((b*y_try)/(b+2*y_try))**(2/3)*(S0)**(1/2)
    deltaQ = Q - Q_t
    y_better = (Q/Q_t)*y_try
    y_best = (y_better + y_try)/2
    y_try = y_best
  return y_try

def determineSequentDepths(y0,q,g):
  ys = (y0/2)*((1+(8*q**2)/(g*y0**3))**(1/2) - 1)
  return ys
  
def jumpOccurs(y0_1,y0_2,y_c):
  if (y0_2 > y_c and y_c > y0_1):
    return True
  else:
    return False

def initFinalBWCDepths(jumpLoc, y0_1, y0_2, ys_1, ys_2):
  if jumpLoc[0] == 'u':
    return (ys_2, y0_2)
  else:
    return (y0_1, ys_1)

def findBWCLength(initDepth, FinlDepth, EPS, Q, b, g, n, c, S0):
  totalL1 = 0
  totalL2 = 0
  its = 2
  deltaTotalL = totalL2
  while deltaTotalL >= EPS * totalL2:
    y1 = initDepth
    totalL1 = totalL2
    totalL2 = 0
    deltaY = (FinlDepth - initDepth) / its
    for i in range(2,its+2):
      y2 = y1 + deltaY
      V1 = Q/(b*y1)
      V2 = Q/(b*y2)
      R1 = (b*y1)/(b+2*y1)
      R2 = (b*y2)/(b+2*y2)
      E1 = y1 + (V1**2)/(2*g)
      E2 = y2 + (V2**2)/(2*g)
      V_avg = (V1 + V2)/2
      R_avg = (R1 + R2)/2
      AS_E = ((V_avg*n)/(c*R_avg**(2/3)))**2
      deltaL = (E2-E1)/(S0 - AS_E)
      y1 = y2
      totalL2 = totalL2 + deltaL
    deltaTotalL = abs(totalL2 - totalL1)
    print(f'iteration: {its}; cumulative length: {totalL2}; length added or subtracted: {deltaTotalL}')
    its +=1
#---------- END FUNCTION DEFINITIONS ---------------


#---------- SETTING VARIABLES -----------------------
#---------- user should only change these values----
b = 15
n = .022
S0_1 = .015
S0_2 = .0005
EPS = .001
Q = 215
c = 1.49
g = 32.2
q = Q/b
y_c = ((q**2)/g)**(1/3)
#---------- END SETTING VARIABLES ----------------------


#---------- SOLVING THE PROBLEM ---------------------
y0_1 = determineNormalDepths(b,n,S0_1,EPS,Q,c)
ys_2 = determineSequentDepths(y0_1,q,g)
y0_2 = determineNormalDepths(b,n,S0_2,EPS,Q,c)
ys_1 = determineSequentDepths(y0_2,q,g)
jump = jumpOccurs(y0_1,y0_2,y_c)

jumpStrg = ''
if jump:
  jumpStrg ='Yes, there is a jump'
else:
  jumpStrg = 'No, there is not a jump'

jumpLocStrng = ''
if jump:
  if y0_2 > ys_2 and ys_1 < y0_1:
    jumpLocStrng = 'upstream'
  else:
    jumpLocStrng = 'downstream'

initFinal = initFinalBWCDepths(jumpLocStrng, y0_1,y0_2,ys_1,ys_2)

S0 = float
if jumpLocStrng[0] == 'u':
  S0 = S0_1
else:
  S0 = S0_2

print(f'\ny0_1 = {y0_1} ft and ys_2 = {ys_2} ft')
print(f'y0_2 = {y0_2} ft and ys_1 ={ys_1} ft\n')
print(jumpStrg)
if jump:
  print('The jump is ' + jumpLocStrng + '\n')
  print(f'initial depth is {initFinal[0]} and final depth is + {initFinal[1]} \n')
  findBWCLength(initFinal[0], initFinal[1], EPS, Q, b, g, n, c, S0)

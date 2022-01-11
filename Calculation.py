from icecream import ic
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
# !Define (x,y) phase A,B,C
xa = 0  # phase A
ya = 15.1
xb = 4  # phase B
yb = 15.1
xc = 4  # phase C
yc = 12.6
r_a = 115
theta_a = 0
r_b = 115
theta_b = -120
r_c = 115
theta_c = 120
r = 12.825  # radius of cable or GMR
line_current = 1606.539
EPSILON_0 = 8.854*pow(10, -12)

#define function polar form <--> cartesian


def pol2cart(r, theta):
    z = r * np.exp(1j * theta)
    x, y = z.real, z.imag
    return x, y


def cart2pol(x, y):
    z = x + y * 1j
    r, theta = np.abs(z), np.angle(z)*57.2958
    # 1rad = 57.2958 deg
    return r, theta

#*define function find Distance between 2 point


def distance(Xa, Xb, Ya, Yb):
  distance = np.sqrt((Xa-Xb)**2+(Ya-Yb)**2)
  return distance

def calculate(xp, yp, size):
    xp=-xp
    xv_a, yv_a = pol2cart(r_a, np.radians(theta_a))
    xv_b, yv_b = pol2cart(r_b, np.radians(theta_b))
    xv_c, yv_c = pol2cart(r_c, np.radians(theta_c))
    v_complex = np.array([[complex(xv_a, yv_a)/np.sqrt(3)], [complex(xv_b, yv_b)/np.sqrt(3)], [complex(xv_c, yv_c)/np.sqrt(3)]])
    dab = distance(xa, xb, ya, yb)
    dac = distance(xa, xc, ya, yc)
    dbc = distance(xb, xc, yb, yc)
    dab_p = distance(xa, xb, ya, -yb)
    dac_p = distance(xa, xc, ya, -yc)
    dbc_p = distance(xb, xc, yb, -yc)
    Matrix = np.array([[np.log(2*ya/r), np.log(dab_p/dab), np.log(dac_p/dac)],
                        [np.log(dab_p/dab), np.log(2*yb/r), np.log(dbc_p/dbc)],
                        [np.log(dac_p/dac), np.log(dbc_p/dbc), np.log(2*yc/r)]])
    q_cart = 2*np.pi*EPSILON_0*np.matmul(np.linalg.inv(Matrix), v_complex)
    dpa = distance(xp, xa, yp, ya)
    dpb = distance(xp, xb, yp, yb)
    dpc = distance(xp, xc, yp, yc)
    dpa_p = distance(xp, xa, yp, -ya)
    dpb_p = distance(xp, xb, yp, -yb)
    dpc_p = distance(xp, xc, yp, -yc)
    vpe = ((q_cart[0]*np.log(dpa_p/dpa))+(q_cart[1]*np.log(dpb_p/dpb))+(q_cart[2]*np.log(dpc_p/dpc)))/(2*np.pi*EPSILON_0)

    #====== MAGNETIC FIELD ======
    Dpa = distance(xp, xa, 0, ya)
    Dpb = distance(xp, xb, 0, yb)
    Dpc = distance(xp, xc, 0, yc)

    Ia, Ima = pol2cart(line_current, np.radians(theta_a))
    Ib, Imb = pol2cart(line_current, np.radians(theta_b))
    Ic, Imc = pol2cart(line_current, np.radians(theta_c))

    Iphase = np.array([[complex(Ia, Ima)], [complex(Ib, Imb)], [complex(Ic, Imc)]])

    SuperPosition = np.array([[np.log(Dpa/dpa), np.log(Dpb/dpb), np.log(Dpc/dpc)]])

    matrix2 = np.matmul(SuperPosition, Iphase)

    ep = 2*(10**-7)*100*np.pi*matrix2
    vpm = ep*size
    vp = vpm+vpe
    (VP, VI) = cart2pol(np.real(vp), np.imag(vp))
    #print(f"Induce Voltage is {VP[0][0]:.2f}∠{VI[0][0]:.2f}°")
    return round(VP[0][0], 2)

m_xp = []
v_induce = []
for x in np.arange(0.1, 10, 0.1):
  x = round(x, 1)
  m_xp.append(x)
  v_induce.append(calculate(x, 15.1, 40))
data = {'distance': m_xp, 'Induced Voltage': v_induce}
df = pd.DataFrame.from_dict(data)
df.set_index('distance', inplace=True)
ic(df)
df.plot(figsize=(10, 8), ylabel='Induce Voltage(V)', xlabel='Distance(m)', title='Induce Voltage & Distance', grid=True, xlim=0, ylim=0)
plt.show()

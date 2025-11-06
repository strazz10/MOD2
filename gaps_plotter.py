import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
from matplotlib.gridspec import GridSpec

########################################################################
###################Prepare yourself for some############################
########################shitty ass code#################################
########################################################################

data_dir = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/DatiPBC/"  #trova e carica file
file_pattern = "*.txt" 

xx = np.linspace(3, 22, 50)
xx2 = np.linspace(10, 22, 15)

files = sorted(glob.glob(data_dir + file_pattern))
if not files:
    raise FileNotFoundError("No .txt files found in directory.")

print(f"Found {len(files)} files.")

fig = plt.figure(figsize=(14, 10))
label = [ 5, 11, 15, 19]   #misure della catena
gs = GridSpec(2, 2, figure=fig)
axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#ff7f0e', '#2ca02c', '#1f77b4']

myrand1 = np.random.randint(0, len(colors))
myrand2 = np.random.randint(0, len(colors))
myrand3 = np.random.randint(0, len(colors))

for index,filename in enumerate(files):
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]    						   #campo g
        data[:, 1:4] = np.sort(data[:, 1:4], axis=1) #sistemo lo sbocco del davidson
        #print(data[:, 1:])
        mask1 = (x >= 1.0)
        mask2 = (x <= 1.0)
        
        e0 = data[:, 1]
        e1 = data[:, 2]
        e2 = data[:, 3]
        #e3 = data[:, 4]
        
        d0 = e1-e0
        d1 = e2-e0 #???
        #d2 = e3-e0
        
        axs[index].plot(x[mask1], 2*(x[mask1]-1), ls='-', color='b', label=r'$2(g-1)$')
        axs[index].plot(x[mask2], -4*(x[mask2]-1), ls='-', color='g', label=r'$-4(g-1)$')
        
        axs[index].plot(
                x, d0, 
                marker='o',
                color='k',
                ls='none',
                label= r'$\Delta_0$'
                )
        axs[index].plot(
                x, d1, 
                marker='^',
                color='r',
                ls='none',
                label= r'$\Delta_1$'  
                )
        #axs[index].plot(
                #x, d2, 
                #marker='h',
                #color='b',
                #ls='none',
                #label= 'r$\Delta_2$'  
                #)
        
    except Exception as e:
        print(f"Error reading {filename}: {e}")

xlabels = [r"g"]

for i, ax in enumerate(axs):
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(fontsize=9, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.set_title(f'L= {label[i]} ')
    
axs[len(axs)-1].set_xlabel(xlabels[0], fontsize=10)
axs[len(axs)-2].set_xlabel(xlabels[0], fontsize=10)
fig.suptitle('PBC Energy Gaps')

plt.tight_layout()
plt.show()

##ripeto per obc

data_dir = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/DatiOBC/"  #trova e carica file
file_pattern = "*.txt" 

files = sorted(glob.glob(data_dir + file_pattern))
if not files:
    raise FileNotFoundError("No .txt files found in directory.")

print(f"Found {len(files)} files.")

fig = plt.figure(figsize=(14, 10))
label = [ 6, 10, 14, 18]   #misure della catena
gs = GridSpec(2, 2, figure=fig)
axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#ff7f0e', '#2ca02c', '#1f77b4']

myrand1 = np.random.randint(0, len(colors))
myrand2 = np.random.randint(0, len(colors))
myrand3 = np.random.randint(0, len(colors))

for index,filename in enumerate(files):
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]    						   #campo g
        data[:, 1:4] = np.sort(data[:, 1:4], axis=1) #sistemo lo sbocco del davidson
        #print(data[:, 1:])
        mask1 = (x >= 1.0)
        mask2 = (x <= 1.0)
        
        e0 = data[:, 1]
        e1 = data[:, 2]
        e2 = data[:, 3]
        #e3 = data[:, 4]
        
        d0 = e1-e0
        d1 = e2-e0 #???
        #d2 = e3-e0
        
        axs[index].plot(x[mask1], 2*(x[mask1]-1), ls='-', color='b', label=r'$2(g-1)$')
        axs[index].plot(x[mask2], -2*(x[mask2]-1), ls='-', color='g', label=r'$-4(g-1)$')
        
        axs[index].plot(
                x, d0, 
                marker='o',
                color='k',
                ls='none',
                label= r'$\Delta_0$'
                )
        axs[index].plot(
                x, d1, 
                marker='^',
                color='r',
                ls='none',
                label= r'$\Delta_1$'  
                )
        #axs[index].plot(
                #x, d2, 
                #marker='h',
                #color='b',
                #ls='none',
                #label= 'r$\Delta_2$'  
                #)
        
    except Exception as e:
        print(f"Error reading {filename}: {e}")

xlabels = [r"g"]

for i, ax in enumerate(axs):
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(fontsize=9, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.set_title(f'L= {label[i]} ')
    
axs[len(axs)-1].set_xlabel(xlabels[0], fontsize=10)
axs[len(axs)-2].set_xlabel(xlabels[0], fontsize=10)
fig.suptitle('OBC Energy Gaps')

plt.tight_layout()
plt.show()

##roba a g fissato

fig = plt.figure(figsize=(14, 10))
gs = GridSpec(1, 2, figure=fig)
axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(2)]

try:
	filename1 = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/DatiPBC/gaps_0.75_d0_pbc.dat"  #trova e carica file (ripetere 3 volte)
except Exception as e:
	print(f"Error reading {filename}: {e}")
data = np.loadtxt(filename1)
x = data[:, 0]    						             #lunghezza L
data[:, 1:4] = np.sort(data[:, 1:4], axis=1)         #sistemo lo sbocco del davidson

e0 = data[:, 1]
e1 = data[:, 2]
d0 = e1-e0
axs[0].plot(x, d0, ls='none', color='b', marker='v',label = r'g=0.75')
axs[1].plot(x, d0, ls='none', color='b', marker='v',label = r'g=0.75')

def delta_exp(x,b,c):
	return b*np.exp(c*x)
maskk=(x>14)
popt1, pcov1 = curve_fit(delta_exp, x[maskk], d0[maskk],p0=(0.2, -0.5), maxfev=5000)
print(f"\n g<1 Exponential parameters: {popt1} ± {np.sqrt(np.diag(pcov1))}")

try:
	filename2 = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/DatiPBC/gaps_1.00_d0_pbc.dat" 
except Exception as e:
	print(f"Error reading {filename2}: {e}") 
data = np.loadtxt(filename2)
x = data[:, 0]    						             
data[:, 1:4] = np.sort(data[:, 1:4], axis=1)        

def power_law(x, a, b, c):
	return a+b*x**c
	
maskk=(x>14)

e0 = data[:, 1]
e1 = data[:, 2]
d0 = e1-e0
axs[0].plot(x, d0, ls='none', color='r', marker='o',label = r'g=1.00')
axs[1].plot(x, d0, ls='none', color='r', marker='o',label = r'g=1.00')

popt2, pcov2 = curve_fit(power_law, x[maskk], d0[maskk], p0=(10,1,-1), maxfev=50000)
print(f"\nCritical exponent for g=1: {popt2[2]} ± {np.sqrt(pcov2[2,2])}")

try:
	filename3 = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/DatiPBC/gaps_1.25_d0_pbc.dat"  
except Exception as e:
	print(f"Error reading {filename3}: {e}")
data = np.loadtxt(filename3)
x = data[:, 0]    						             
data[:, 1:4] = np.sort(data[:, 1:4], axis=1)   

e0 = data[:, 1]
e1 = data[:, 2]
d0 = e1-e0
axs[0].plot(x, d0, ls='none', color='k', marker='h',label = r'g=1.25')
axs[1].plot(x, d0, ls='none', color='k', marker='h',label = r'g=1.25')

maskk=(x>13)     
def power_exp(x,a,c):
	return a+np.exp(c*x)
popt3, pcov3 = curve_fit(power_exp, x[maskk], d0[maskk], p0=(0.1, -0.1), maxfev=50000)
print(f"\n: g>1 Exponential parameters: {popt3} ± {np.sqrt(np.diag(pcov3))}")   

axs[0].set_ylabel(r'$\Delta_0$', fontsize=14)
axs[1].set_xscale('log')
axs[1].minorticks_on()
axs[1].grid(which='minor', linestyle='--')

axs[0].plot(xx2, power_exp(xx2, *popt3), color='y', label='Paramagnetic Zone (Exponential+Corrections)')
axs[0].plot(xx2, delta_exp(xx2, *popt1), color='m', label='Ferromagnetic Zone (Exponential)')
axs[1].plot(xx, power_law(xx, *popt2), color='c', label=' Criticality (Power Law)')

for i, ax in enumerate(axs):
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(fontsize=12, frameon=False)
    ax.tick_params(axis='both',which='minor', labelsize=9)
    ax.set_xlabel(r'$L$', fontsize=14)
    ax.set_yscale('log')
    
plt.tight_layout()
plt.show()

##cpu times

try:
	filenamecpu = "/home/strazz/Magistrale/NumMeth/MOD2/Progetto/cputimes.txt"  #trova e carica file (ripetere 3 volte)
except Exception as e:
	print(f"Error reading {filenamecpu}: {e}")

fig = plt.figure(figsize=(10, 7))
data = np.loadtxt(filenamecpu)

cputimes = data[:,0]
lenght = data[:,1]
mask=(lenght>15)

cputimes_fit = cputimes[mask]
lenght_fit = lenght[mask]
def test_exp(x,b,c):
	return -1+b*np.exp(x*c)

plt.plot(lenght, cputimes, marker='^', color='g', ls='none', label='Real Times')
popt, pcov = curve_fit(test_exp, lenght, cputimes, maxfev=5000)
plt.plot(xx, test_exp(xx, *popt), color='#1f77b4', label='Best Fit') 

print(f"\nCpu fit parameters: {popt} ± {np.sqrt(np.diag(pcov))}")

plt.grid(True, linestyle="--", alpha=0.7)
plt.tick_params(axis='both', which='major', labelsize=9)
plt.xlabel(r"$L$")
plt.ylabel("CPU Time (s)")
plt.tight_layout(pad=3.0)
plt.legend()
plt.title(r'CPU Times, PBC, $g=0.75$')

plt.show()

##krylov
#dati (possibilmente da importare dal file)
values = np.array([
    -0.026189639210315363, -0.04851139626316581, -0.03209264025372249, -0.026345120925725496, -0.05746394117613818,
    -0.0009419793505003327, -0.0008774298485150211, -0.0005409015611803625, -0.0005604971647699131, -0.0004417837863002205,
    -1.6782766579126474e-05, -1.3396970643952955e-05, -7.736945008218754e-06, -6.474335350503679e-06, -3.648809979495127e-05,
    -5.273113856674172e-07, -1.284097379539162e-07, -2.892556949518621e-07, -2.6422912924317643e-07, -1.3131466403137892e-07,
    -3.3887772588059306e-09, -3.4178810892626643e-09, -2.954948286060244e-09,
    -3.3651303965598345e-11, -5.729816621169448e-11, -3.637978807091713e-11,
    5.4569682106375694e-12, 1.8189894035458565e-12, 1.6370904631912708e-11
])
scales = np.array([5]*5 + [6]*5 + [7]*5 + [8]*5 + [9]*3 + [10]*3 + [11]*3)

plt.figure(figsize=(8,5))
plt.scatter(scales, np.abs(values), color='r', marker='h')
plt.yscale('log')
plt.xlabel('Krylov Space Dimension', fontsize=14)
plt.ylabel('Absolute Error (Log)', fontsize=14)
plt.grid(True, which='both', ls='--', lw=0.5)

plt.show()

########################################################################
#############################suffering##################################

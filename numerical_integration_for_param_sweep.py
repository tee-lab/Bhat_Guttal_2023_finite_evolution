'''
Author: Shikhara Bhat
Email ID: shikharabhat@gmail.com
Date created: 2023-07-04 10:56:35
Date last modified: 2023-08-17 11:05:07
Purpose: Numerically solve for stationary state of the resource competition model 
that we introduce as an example in Bhat and Guttal 2023

Note: The parameter sweep is SLOW. If you want to directly import the results,
run this script from line 76 onwards
'''
import numpy as np
import pandas as pd
from scipy import integrate as int
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

#fixed parameters
mu = 0.001 #mutation rate


#define some useful functions
def f(p,K,eps): #drift
    return p*((1-p)**2)*((2-eps*(K-1))/K)+mu*(1-2*p)*(1-(1/K))

def g(p,K,eps): #square of diffusion
    return (p*(1-p)*(2+(2+eps)*(-p**2 + 2*p - 1))+mu*(1-3*p*(1-p)))/K

def dU(p,K,eps): #this is the function we need to integrate
    return f(p,K,eps)/g(p,K,eps)

def U(p, K, eps, p0=0): #potential function

    return int.quad(dU, p0, p, args=(K,eps))[0] #output of quad is (result, estimated error). We only want the result.

def m(p, K, eps): #speed density

    return np.exp(U(p,K,eps))/g(p,K,eps)

#################################################
#Computing the probabilities we want

#parameter space to scan over
#Klist = np.arange(300,1301,25)
#epslist = np.arange(0,0.0041,0.0001)
Klist = np.arange(300,1301,1)
epslist = np.arange(0,0.0041,0.0000041)

#this will store the values for us to plot
plot_matrix = np.zeros((len(Klist),len(epslist)))

#This is SLOW. Running time approx 2 days on a Ryzen 7 6800H with 8 cores
#I'm sure I could have made the code faster with parallelization, but brute
#force was quicker to just implement and run
for i in range(len(Klist)):
    for j in range(len(epslist)):

        #probability that freq of type 1 individuals is >= 0.5
        plot_matrix[i,j] = 1-((int.quad(m, 0, 0.5, args=(Klist[i],epslist[j]))[0])/(int.quad(m, 0, 1, args=(Klist[i],epslist[j]))[0])) 

########################################
#Convert to dataframe
df = pd.DataFrame(plot_matrix)
df.index = Klist #rows index K values
df.columns = np.round(epslist,7) #columns index epsilon values
#df.columns = np.round(epslist,5) #columns index epsilon values


df.to_csv('data/numerical_eps_K_param_sweep.csv')

########################################
#Plotting
#import data (in case you don't want to run the above code each time, just
# run the script from here on. Comment out the line below if you want to simulate each time)
df = pd.read_csv('data/numerical_eps_K_param_sweep.csv',index_col=[0],header=[0])
df.index = Klist #rows index K values
df.columns = np.round(epslist,7) #columns index epsilon values

#df = gaussian_filter(df,sigma=3)

#Make all labels bold
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.titleweight'] = "bold"


#plot        
plt.figure(figsize=[12, 6])

#twoslopenorm makes sure there are two different scales for above vs below 0.5
#integer values of xticklabels and yticklabels samples every nth element
hmap = sns.heatmap(df.T,cmap='icefire',norm=colors.TwoSlopeNorm(vcenter=0.5,vmin=0.3,vmax=0.6),xticklabels=200,yticklabels=200)
#Add a white line seperating <0.5 from >0.5. We can do this using a mask
hmap = sns.heatmap(df.T, mask=(np.abs(df.T - 0.5) > 0.0002), cmap=colors.ListedColormap(['white']), cbar=False, ax=hmap,xticklabels=200,yticklabels=200)

sampled_yticks = epslist[0::200]

#aesthetics
hmap.invert_yaxis()
hmap.figure.axes[0].set_xlabel('K', size=16)
hmap.figure.axes[0].set_ylabel('eps', size=22)
hmap.figure.axes[-1].set_ylabel("integral", size=16)
hmap.set_xticklabels(hmap.get_xmajorticklabels(), fontsize = 16)
#hmap.set_yticklabels(hmap.get_ymajorticklabels(), fontsize = 16)
hmap.set_yticklabels([f'{x:.1e}' for x in sampled_yticks], fontsize = 16)
plt.yticks(rotation=0) 


save_loc = 'plots/'
plt.savefig(save_loc+'K_epsilon_numerical_param_sweep.png',dpi=1000)

plt.show()
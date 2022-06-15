#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nik_Sharma, June 2020
"""

import seaborn as sns
import scipy.stats as scs
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
style.use('seaborn-talk')
import decimal
from astropy.visualization import hist

# =============================================================================
# Importing data from excel file into Python as a dataframe
# =============================================================================

'''
Enter file path of the excel file in which data is stored 
'''

# for windows #

data = pd.read_excel(r'FILE - PATH') #file path of .xls to be entered

#data = pd.read_excel(r'')

WC = pd.DataFrame(data, columns = ['Length'])

# =============================================================================
# Calculating iterative median
# =============================================================================

d50_values = WC['Length'].values
iterative_median_photo_wc = np.zeros((200,1)) #'200' is the total number of measurements; it will change accordingly

for i in range(200):
    iterative_median_photo_wc[i] = np.median(d50_values[0:i+1])
    
WC['Iterative Median'] = iterative_median_photo_wc


# =============================================================================
# Calculating cumulative percentage of wolman count data
# =============================================================================

WC['cumulative_percentage'] = 100*WC.Length.cumsum()/WC.Length.sum()

# =============================================================================
# Normalising wolman count data according to the psi scale (psi scale = positive of phi scale)
# =============================================================================

WC['normalised_wolman($\psi$ scale)'] = np.log2(WC.Length)

# =============================================================================
# D in mm from psi scale (psi scale = positive of phi scale, D = 2 ** psi scale)
# =============================================================================

WC['wolman_mm_from_psi'] = 2**WC['normalised_wolman($\psi$ scale)']

# =============================================================================
# Calculation of D50 and other important variables
# =============================================================================

# Shapiro-Wilk test calculation #

Shapiro_Wilk = scs.shapiro(WC['normalised_wolman($\psi$ scale)'])
Shapiro_Wilk_Rounded = round(Shapiro_Wilk[0],2),round(Shapiro_Wilk[1],2)

# D50 (psi) and D50 (mm) calculation #

D_Fifty_psi = np.around(np.percentile(WC['normalised_wolman($\psi$ scale)'],50),2)

D_Fifty_mm = np.around(2**D_Fifty_psi,2)

# D84 (psi) and D84 (mm) calculation #

D_EightyFour_psi = np.around(np.percentile(WC['normalised_wolman($\psi$ scale)'],84),2)

D_EightyFour_mm = np.around(2**D_EightyFour_psi,2)

# standard error calculation psi and mm calculation  #

Standard_Error_psi = round((decimal.Decimal(scs.sem(WC['normalised_wolman($\psi$ scale)']))),2)

Standard_Error_mm = round(decimal.Decimal(scs.sem(WC['wolman_mm_from_psi'])),2)


# standard deviation (sigma) calculation psi and mm calculation  #

Standard_Deviation_psi = round((decimal.Decimal(np.std(WC['normalised_wolman($\psi$ scale)']))),2)

Standard_Deviation_mm = round(decimal.Decimal(np.std(WC['wolman_mm_from_psi'])),2)

# =============================================================================
# Skewness and Kurtosis calculation
# =============================================================================

Skewness = np.around(scs.skew(WC['normalised_wolman($\psi$ scale)']),2)


Kurtosis = np.around(scs.kurtosis(WC['normalised_wolman($\psi$ scale)']),2)


# =============================================================================
# Plotting of Histogram and various other subplots 
# =============================================================================

figsize = (14,14)

fig = plt.figure(figsize=figsize)
fig.subplots_adjust(hspace=0.8, wspace=0.7)
plt.suptitle('Location: Olson, Sample  - Wolman Count (100% Representative)', fontsize='xx-large')

plt.subplot(3,3,1)
plt.grid(False)
plt.axis('off')

plt.subplot(3,3,2)
hist(WC['Length'],bins="scott",alpha=0.4,rwidth=0.6, color='b')
#plt.legend(loc='best')
#plt.xlim(0,100)
#plt.ylim(0,70)
plt.xlabel('Clast Size [D (mm)]', fontsize='large')
plt.ylabel('Frequency', fontsize='large')
plt.title('Raw Data', fontsize='x-large')

plt.subplot(3,3,3)
plt.grid(False)
plt.axis('off')

plt.subplot(3,3,4)
sns.distplot(WC['normalised_wolman($\psi$ scale)'],kde=True, rug=True, kde_kws={"lw":2},color='b')
plt.xlim(0,8)
plt.xlabel('$\psi$ Scale', fontsize='large')
plt.ylabel('Desity', fontsize='large')
plt.title('Normalised Data Distribution', fontsize='x-large')
mm = lambda psi: 2**psi
xmin, xmax = plt.gca().get_xlim()
ax1 = plt.twiny()
ax1.set_xlabel("$D \: [log(mm)]$")
ax1.set_xlim((mm(xmin),mm(xmax)))
ax1.set_xscale('log')

plt.subplot(3,3,5)
sns.distplot(WC['normalised_wolman($\psi$ scale)'],hist_kws=dict(cumulative=True, alpha=0),kde_kws=dict(cumulative=True,color='b'))
plt.xlim(0,8)
plt.xlabel('$\psi$ Scale', fontsize='large')
plt.ylabel('CDF', fontsize='large')
#plt.xlabel('$D \: [log(mm)]$', fontsize='large')
plt.title('Cumulative Distribution Function (CDF)', fontsize='x-large')
mm = lambda psi: 2**psi
xmin, xmax = plt.gca().get_xlim()
ax1 = plt.twiny()
ax1.set_xlabel("$D \: [log(mm)]$")
ax1.set_xlim((mm(xmin),mm(xmax)))
ax1.set_xscale('log')
plt.annotate('$D_{50}$',xy=(D_Fifty_psi,0.50),xytext=(D_Fifty_mm+0.6,0.03),color='r')
plt.axhline(0.50,ls='-.',lw=0.5,color='b')
plt.axvline(D_Fifty_mm,ls='-.',lw=0.5,color='b')
percent = lambda percent: percent*100
ymin, ymax = plt.gca().get_ylim()
ax2 = plt.twinx()
ax2.set_ylabel("$\%$")
ax2.set_ylim((percent(ymin),percent(ymax)))


ax1 = plt.subplot(3,3,6)
probplot = sm.ProbPlot(WC['normalised_wolman($\psi$ scale)'],scs.t,fit=True)
probplot.qqplot(ax=ax1,line='45')
plt.title('Normal Q-Q Plot')
#plt.xlim(-4,4)
#plt.ylim(-4,4)


# =============================================================================
# Creating a table
# =============================================================================

Table_Data = [[Shapiro_Wilk_Rounded,Skewness,Kurtosis,D_Fifty_psi,D_Fifty_mm,D_EightyFour_psi,D_EightyFour_mm,Standard_Error_mm,Standard_Deviation_mm]]
        
#Table_Data = [[Shapiro_Wilk_Rounded,Skewness,Kurtosis,D_Fifty_psi,D_Fifty_mm,Standard_Error_mm,Standard_Deviation_mm]]

rows = ['']

columns = ['Shapiro-Wilk Test','Skewness','Kurtosis', '$D_{50} \: (\psi)$', '$D_{50} \: (mm)$', '$D_{84} \: (\psi)$', '$D_{84} \: (mm)$','$\pm \Delta D_{50} \: (mm)$','$\sigma \: (mm)$']

#columns = ['Shapiro-Wilk Test','Skewness','Kurtosis', '$D_{50} \: (\psi)$', '$D_{50} \: (mm)$','$\pm \Delta D_{50} \: (mm)$','$\sigma \: (mm)$']

# =============================================================================
# Displaying the table as subplot
# =============================================================================


plt.subplot(3,3,8)
plt.axis('off')
table = plt.table(cellText=Table_Data,rowLabels=rows,colLabels=columns,colWidths=[0.40]*9,cellLoc='center',loc='top')
table.auto_set_font_size(False)
table.set_fontsize(14)
table.scale(2,2)

# =============================================================================
# Exporting the final figure as a pdf
# =============================================================================


plt.savefig("_WC.pdf",bbox_inches='tight')
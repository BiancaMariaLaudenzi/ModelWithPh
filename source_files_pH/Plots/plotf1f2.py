import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# # Leggi il file CSV
data = pd.read_csv('f1f2_vs_pO2pCO2_wider.csv')

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
heatmap_data = data.pivot_table(index='pO2', columns='pCO2', values='f1') 
heatmap_data = heatmap_data.sort_index(axis=0).sort_index(axis=1)
pO2 = heatmap_data.index.values   # asse Y
pCO2 = heatmap_data.columns.values  # asse X
Z = heatmap_data.values
X, Y = np.meshgrid(pCO2, pO2)
ax = sns.heatmap(heatmap_data, cmap='coolwarm', annot=False, cbar_kws={'label': 'f1'})
plt.title('f1')
plt.xlabel('pCO2 [kPa]')  
plt.ylabel('pO2 [kPa]')
contours = plt.contour(X, Y, Z, levels=[0], colors='black', linewidths=2)
plt.clabel(contours, fmt='f1=0', fontsize=10)
plt.grid(True)

plt.subplot(1,2,2)
heatmap_data2 = data.pivot_table(index='pO2', columns='pCO2', values='f2') 
heatmap_data2 = heatmap_data2.sort_index(axis=0).sort_index(axis=1)
pO2 = heatmap_data2.index.values   # asse Y
pCO2 = heatmap_data2.columns.values  # asse X
Z2 = heatmap_data2.values
X2, Y2 = np.meshgrid(pCO2, pO2)
ax = sns.heatmap(heatmap_data2, cmap='coolwarm', annot=False, cbar_kws={'label': 'f2'})
plt.title('f2')
plt.xlabel('pCO2 [kPa]')  
plt.ylabel('pO2 [kPa]')
contours2 = plt.contour(X2, Y2, Z2, levels=[0], colors='black', linewidths=2)
plt.clabel(contours2, fmt='f2=0', fontsize=10)
plt.grid(True)

plt.tight_layout()
plt.savefig('f1f2_vs_pO2pCO2_old.png', dpi=300)
plt.show()
exit(-1)

# pO2_vals = np.sort(data['pO2'].unique())
# pCO2_vals = np.sort(data['pCO2'].unique())

# # Creo griglia
# X, Y = np.meshgrid(pCO2_vals, pO2_vals)
# f1 = data['f1']
# f2 = data['f2']
# Z1 = f1.values.reshape(len(pO2_vals), len(pO2_vals))
# Z2 = f2.values.reshape(len(pO2_vals), len(pO2_vals))

# # --- Plot f1 ---
# plt.figure(figsize=(12,5))
# plt.subplot(1,2,1)
# contourf1 = plt.contourf(X, Y, Z1, levels=50, cmap='coolwarm')
# plt.colorbar(contourf1, label='f1 value')
# contour0 = plt.contour(X, Y, Z1, levels=[0], colors='black', linewidths=2)
# plt.clabel(contour0, fmt='f1=0', colors='black')
# plt.title('Heatmap of f1 (cO2 error)')
# plt.xlabel('pCO2 [kPa]')
# plt.ylabel('pO2 [kPa]')
# plt.grid(True)

# # --- Plot f2 ---
# plt.subplot(1,2,2)
# contourf2 = plt.contourf(X, Y, Z2, levels=50, cmap='coolwarm')
# plt.colorbar(contourf2, label='f2 value')
# contour0 = plt.contour(X, Y, Z2, levels=[0], colors='black', linewidths=2)
# plt.clabel(contour0, fmt='f2=0', colors='black')
# plt.title('Heatmap of f2 (cCO2 error)')
# plt.xlabel('pCO2 [kPa]')
# plt.ylabel('pO2 [kPa]')
# plt.grid(True)

# plt.tight_layout()
# plt.savefig('f1f2_vs_pO2pCO2_wider.png', dpi=300)
# plt.show()

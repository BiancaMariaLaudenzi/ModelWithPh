import csv
import matplotlib.pyplot as plt

cCO2_vals = []
f_vals = []
with open('f_vs_cCO2.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cCO2_vals.append(float(row['cCO2']))
        f_vals.append(float(row['f']))

cCO2_vals_NEWTON1 = []
f_vals_NEWTON1 = []
with open('f_vs_cCO2_NEWTON1.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cCO2_vals_NEWTON1.append(float(row['cCO2']))
        f_vals_NEWTON1.append(float(row['f']))

cCO2_vals_NEWTON2 = []
f_vals_NEWTON2 = []
with open('f_vs_cCO2_NEWTON2.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cCO2_vals_NEWTON2.append(float(row['cCO2']))
        f_vals_NEWTON2.append(float(row['f']))

cCO2_vals_NEWTON3 = []
f_vals_NEWTON3 = []
with open('f_vs_cCO2_NEWTON3.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cCO2_vals_NEWTON3.append(float(row['cCO2']))
        f_vals_NEWTON3.append(float(row['f']))

cCO2 = []
pH = []
with open('cCO2vsPH.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        cCO2.append(float(row['cCO2'])*1000)
        pH.append(float(row['pH']))

x = []
y = []
with open('pO2vsSaO2_pH7punto4.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        x.append(float(row['pO2'])/ 0.133322)
        y.append(float(row['SaO2']))

x1 = []
y1 = []
with open('pCO2vscCO2_pO213.3.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        x1.append(float(row['pCO2'])/ 0.133322)
        y1.append(float(row['cCO2'])*1000)

x2 = []
y2 = []
with open('pCO2vscCO2_pO25.3.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        x2.append(float(row['pCO2'])/ 0.133322)
        y2.append(float(row['cCO2'])*1000)

plt.figure()
# plt.plot(x, y, 'g-', label= "pH=7.4")
plt.plot(x1, y1, 'r-', label= "pO2 = 5.3 kPa")
plt.plot(x2, y2, 'b-', label= "pO2 = 13.3 kPa")
plt.xlabel('pCO2 [mmHg]')
plt.ylabel('cCO2 [mmol/L]')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('cCO2_vs_pCO2.png')



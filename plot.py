import matplotlib.pyplot as plt
import math as m

with open('./SimulationsResults/out.txt') as f:
    data = f.readlines()[2:]

data = [row.strip().split(' ') for row in data]
x = [abs(float (row[0])) for row in data]
y = [abs( float (row[1])) for row in data]

plt.scatter(x,y,label='Results from BP')
plt.title('Results for')
plt.xlabel('Temperature')
plt.ylabel('Magnitude of magnetization')
plt.legend()
plt.grid(True)
plt.show()
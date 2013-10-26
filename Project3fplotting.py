import numpy
from matplotlib.pylab import *

RunID = '0'

#reading file
Planets = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
figure()
hold('on')

for planet in Planets:
    f = open(''.join([RunID, '_', planet, '.dat']))
    data = [[], []]
    for line in f:
        elem = line.strip().split(' ')
        data[0].append(float(elem[0]))
        data[1].append(float(elem[-1]))
    f.close()
    plot(data[0], data[1])


title(r"Motion of the planets in the solar system (and Pluto)")
legend(Planets)
xlabel('x[AU]')
ylabel('y[AU]')
axis('equal')
show()

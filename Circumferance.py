import math
from math import pi
import numpy as np
import matplotlib.pyplot as plt

#rysuje 200 detektorów z rozpietoscią 90 stopni

rozpietosc_katowa_punktow = 90  #w stopniach, w radianach: alfa*math.pi/180
liczba_emiterow = 200

img_size = 256


def points_on_circumference(center, r, n=round(liczba_emiterow*360/rozpietosc_katowa_punktow)):
    return [
        (
            center[0] + (math.cos(2*pi/n * x) * r),  # x
            center[1] + (math.sin(2*pi/n * x) * r)  # y

        ) for x in range(math.floor(-liczba_emiterow/2), math.ceil(liczba_emiterow/2))]


image = np.zeros((img_size+1, img_size+1))
for x,y in points_on_circumference(center=(img_size/2, img_size/2), r=img_size/2):
    image[round(x)][round(y)] = 255
plt.subplot(1, 1, 1)
plt.imshow(image, cmap='gray')
plt.show()

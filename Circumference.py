import math
from math import pi
import numpy as np
import matplotlib.pyplot as plt

#rysuje 200 detektorów z rozpietoscią 90 stopni

rozpietosc_katowa_punktow = math.radians(90)  #w stopniach, w radianach: alfa*math.pi/180
liczba_emiterow = 200

img_size = 256


def points_on_circumference(center, r, n=round(liczba_emiterow*360/rozpietosc_katowa_punktow)):
    return [
        (
            center[0] + (math.cos(2*pi/n * x) * r),  # x
            center[1] + (math.sin(2*pi/n * x) * r)  # y

        ) for x in range(math.floor(-liczba_emiterow/2), math.ceil(liczba_emiterow/2))]

#def emiter(center, r):
#    x,y = center
#    return x-r, y

def odbiorniki(r, nr_odbiornika, alfa=0):
    x = r * math.cos(alfa + math.pi - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    y = r * math.sin(alfa + math.pi - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    return x+r, y+r     #srodek nie w (0,0) a (r,r)

def emitery(r, nr_odbiornika, alfa=0):
    x = r * math.cos(alfa - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))
    y = r * math.sin(alfa - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))

    return x+r, y+r     #srodek nie w (0,0) a (r,r)


img = np.zeros((img_size+1, img_size+1))
for alpha in range(0, 360, 10):
    image = np.copy(img)
    for i in range(liczba_emiterow):
        x,y = odbiorniki(r=img_size/2, nr_odbiornika=i, alfa=math.radians(alpha))
        image[round(x)][round(y)] = 255
    for i in range(liczba_emiterow):
        x, y = emitery(r=img_size / 2, nr_odbiornika=i, alfa=math.radians(alpha))
        image[round(x)][round(y)] = 255
    image[round(x)][round(y)] = 255
    plt.subplot(1, 1, 1)
    plt.imshow(image, cmap='gray')
    plt.show()

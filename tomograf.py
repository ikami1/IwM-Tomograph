from skimage import io
from matplotlib import pyplot as plt
import numpy as np
import math
from PIL import Image

#name = "example_image.png"
name = "phantom256.gif"
ile_emiterow = 200
#rozpietosc_katowa_punktow = 90  #w stopniach, w radianach: alfa*math.pi/180
kat_kroku = 2


def plotLineLow(x0, y0, x1, y1, img, img2, img_size):
    suma=0
    licznik=1
    dx = x1 - x0
    dy = y1 - y0
    yi = 1
    if dy < 0:
        yi = -1
        dy = -dy
    D = 2*dy - dx
    y = y0

    for x in range(x0,x1):
        if (x>0) and (x<img_size) and (y>0) and (y<img_size):
            suma=suma+img[x,y]
            licznik=licznik+1
        if D > 0:
           y = y + yi
           D = D - 2*dx
        D = D + 2*dy
    return suma/licznik

def plotLineHigh(x0,y0, x1,y1, img, img2, img_size):
    suma=0
    licznik=1
    dx = x1 - x0
    dy = y1 - y0
    xi = 1
    if dx < 0:
        xi = -1
        dx = -dx
    D = 2*dx - dy
    x = x0

    for y in range(y0,y1):
        if (x>0) and (x<img_size) and (y>0) and (y<img_size):
            suma=suma+img[x,y]
            licznik=licznik+1
        if D > 0:
            x = x + xi
            D = D - 2*dy
        D = D + 2*dx
    return suma/licznik

def rotuj_o_kat(x,y,kat,sx,sy):
    for i in range(x.size):
        tx=x[i]
        tx=tx-sx
        ty=y[i]
        ty=ty-sy
        xp=tx*math.cos(kat)-ty*math.sin(kat)
        yp=tx*math.sin(kat)+ty*math.cos(kat)
        x[i]=xp+sx
        y[i]=yp+sy
    return x, y

def init_y(x, sx, sy, rad):
    y_em=[]
    y_czuj=[]
    for i in x:
        y1=sy-rad
        y2=sy+rad

        if y1<y2:
            y_em.append(y1)
            y_czuj.append(y2)
        else:
            y_em.append(y2)
            y_czuj.append(y1)

    return y_em, y_czuj

def bresenham(x0,y0,x1,y1, img, img2, img_size):
    suma=0
    if abs(y1 - y0) < abs(x1 - x0):
        if x0 > x1:
            suma=suma+plotLineLow(x1, y1, x0, y0, img, img2, img_size)
        else:
            suma=suma+plotLineLow(x0, y0, x1, y1, img, img2, img_size)
    else:
        if y0 > y1:
            suma=suma+plotLineHigh(x1, y1, x0, y0, img, img2, img_size)
        else:
            suma=suma+plotLineHigh(x0, y0, x1, y1, img, img2, img_size)
    return suma

def bresenham2TEST(xx1,yy1,xx2,yy2,img,img2,img_size):
    x1 = int(xx1)
    y1 = int(yy1)
    x2 = int(xx2)
    y2 = int(yy2)
    suma=0
    licznik=1
    d=0
    dx=0
    dy=0
    ai=0
    bi=0
    xi=0
    yi=0
    x = x1
    y = y1
    if (x1 < x2):
         xi = 1
         dx = x2 - x1
    else:
         xi = -1
         dx = x1 - x2
    if (y1 < y2):
         yi = 1
         dy = y2 - y1
    else:
         yi = -1
         dy = y1 - y2
    if (x>0) and (x<img_size) and (y>0) and (y<img_size):
        suma=suma+img[x,y]
        licznik=licznik+1
    if (dx > dy):
         ai = (dy - dx) * 2
         bi = dy * 2
         d = bi - dx
         while (x != x2):
             if (d >= 0):
                 x += xi
                 y += yi
                 d += ai
             else:
                 d += bi
                 x += xi
             if (x>0) and (x<img_size) and (y>0) and (y<img_size):
                 suma=suma+img[x,y]
                 #img2[x,y]=127
                 licznik=licznik+1
    else:
         ai = ( dx - dy ) * 2
         bi = dx * 2
         d = bi - dy
         while (y != y2):
             if (d >= 0):
                 x += xi
                 y += yi
                 d += ai
             else:
                 d += bi
                 y += yi
             if (x>0) and (x<img_size) and (y>0) and (y<img_size):
                 suma=suma+img[x,y]
                 #img2[x,y]=127
                 licznik=licznik+1
    return suma/licznik

def get_sinogram(img, img_size, D):
    ile_krokow = round(180/kat_kroku)
    sx = round(D+(img_size/2))
    sy = round(D+(img_size/2))
    r = img_size/2
    x_start=np.linspace(D,D+img_size,ile_emiterow, dtype=int, endpoint=False)
    #emit_y, odbiornik_y=init_y(x_start, sx, sy, r)
    emit_y, odbiornik_y=init_y(x_start, sx, sy, int(img_size/2))
    emit_x=x_start
    odbiornik_x=np.copy(x_start)
    sin_img=np.zeros((ile_krokow,ile_emiterow))
    for krok in range(ile_krokow):
        #sin_img[krok]=krok
        img2=np.copy(img)
        for emit in range(ile_emiterow):
            sin_img[krok][emit]=bresenham2TEST(emit_x[emit],emit_y[emit],odbiornik_x[emit],odbiornik_y[emit],img, img2, img_size)
        #io.imsave("kat"+str(krok)+".jpg",img2)
        emit_x, emit_y=rotuj_o_kat(emit_x,emit_y,kat_kroku*math.pi/180,sx,sy)
        odbiornik_x, odbiornik_y = rotuj_o_kat(odbiornik_x, odbiornik_y, kat_kroku*math.pi/180, sx, sy)

    return sin_img


def inv_radon(img_sin):
    inv_sin = img_sin
    return inv_sin


def main():
    img = io.imread(name, as_grey=True)
    width, height = img.shape
    radius = width/2
    Delt=0
    '''nimg=np.zeros((100+R*2,100+R*2))
    S=int(X/2)
    Delt=50+R-S
    for eh in range(X):
        for yh in range(Y):
            nimg[Delt+eh,Delt+yh]=img[eh,yh]
    img=nimg'''

    if width == height:
        plt.subplot(2, 2, 1)
        plt.title('Obraz wejsciowy')
        plt.axis('off')
        io.imshow(img, cmap="gray")

        plt.subplot(2, 2, 2)
        sinogram = get_sinogram(img, width, Delt)
        #sinogram = sinogram.resize((sinogram.size[0] + 100, sinogram.size[1]), Image.ANTIALIAS)
        #io.imsave("sinogram.jpg",sinogram)
        plt.title('Sinogram')
        io.imshow(sinogram, cmap="gray")

        plt.subplot(2, 2, 4)
        re_img = inv_radon(sinogram)
        plt.title('Obraz wyjsciowy')
        io.imshow(re_img, cmap='gray')

        io.show()
    else:
        print("Obrazek nie jest kwadratowy")


if __name__ == "__main__":
    main()

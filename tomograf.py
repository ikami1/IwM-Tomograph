from skimage import io
from skimage.transform import  resize
from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.fftpack import fft, ifft, fftfreq
from scipy.interpolate import interp1d

#name = "example_image.png"
name = "phantom256.gif"
liczba_emiterow = 600
rozpietosc_katowa_punktow = math.radians(150)  #w stopniach, w radianach: alfa*math.pi/180
kat_kroku = 1


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
            img2[x, y] = 127
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
            img2[x, y] = 127
            licznik=licznik+1
        if D > 0:
            x = x + xi
            D = D - 2*dy
        D = D + 2*dx
    return suma/licznik

def bresenham(x0,y0,x1,y1, img, img2, img_size):
    x0 = round(x0)
    y0 = round(y0)
    x1 = round(x1)
    y1 = round(y1)
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

def bresenham2TEST(x1,y1,x2,y2,img,img2,img_size):
    x1 = round(x1)
    y1 = round(y1)
    x2 = round(x2)
    y2 = round(y2)
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
                 img2[x,y]=127
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
                 img2[x,y]=127
                 licznik=licznik+1
    return suma/licznik

'''
def get_sinogram(img, img_size, D):
    ile_krokow = round(180/kat_kroku)
    sx = round(D+(img_size/2))
    sy = round(D+(img_size/2))
    r = img_size/2
    x_start=np.linspace(D,D+img_size,liczba_emiterow, dtype=int, endpoint=False)
    #emit_y, odbiornik_y=init_y(x_start, sx, sy, r)
    emit_y, odbiornik_y=init_y(x_start, sx, sy, int(img_size/2))
    emit_x=x_start
    odbiornik_x=np.copy(x_start)
    sin_img=np.zeros((ile_krokow,liczba_emiterow))
    for krok in range(ile_krokow):
        #sin_img[krok]=krok
        img2=np.copy(img)
        for emit in range(liczba_emiterow):
            sin_img[krok][emit]=bresenham2TEST(emit_x[emit],emit_y[emit],odbiornik_x[emit],odbiornik_y[emit],img, img2, img_size)
        #io.imsave("kat"+str(krok)+".jpg",img2)
        emit_x, emit_y=rotuj_o_kat(emit_x,emit_y,kat_kroku*math.pi/180,sx,sy)
        odbiornik_x, odbiornik_y = rotuj_o_kat(odbiornik_x, odbiornik_y, kat_kroku*math.pi/180, sx, sy)

    return sin_img
'''

#Polozenie detektorow
def odbiornik(r, nr_odbiornika, alfa=0.0):
    x = r * math.cos(alfa + math.pi - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    y = r * math.sin(alfa + math.pi - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    return x+r, y+r     #srodek nie w (0,0) a (r,r)

# Pojedynczy emiter
def emiter(r, alfa=0.0):
    x = r * math.cos(alfa)
    y = r * math.sin(alfa)
    return x+r, y+r     #srodek nie w (0,0) a (r,r)

# Tyle samo emiterow co detektorow
def emitery(r, nr_odbiornika, alfa=0.0):
    x = r * math.cos(alfa - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))
    y = r * math.sin(alfa - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))

    return x+r, y+r     #srodek nie w (0,0) a (r,r)

def get_sinogram(img, img_size, D):
    ile_krokow = math.floor(180/kat_kroku)
    r = img_size/2

    sin_img = np.zeros((ile_krokow, liczba_emiterow))
    for krok in range(ile_krokow):

        #emiter_x, emiter_y = emiter(r=r, alfa=math.radians(krok * kat_kroku))      #Pojedynczy emiter, promieniowanie stożkowe
        emiters = []
        for i in range(liczba_emiterow):
            emiters.append(emitery(r=r, nr_odbiornika=i, alfa=math.radians(krok * kat_kroku)))

        odbiorniki = []
        for i in range(liczba_emiterow):
            odbiorniki.append(odbiornik(r=r, nr_odbiornika=i, alfa=math.radians(krok * kat_kroku)))

        img2 = np.copy(img)
        for i in range(liczba_emiterow):
            odbiornik_x, odbiornik_y = odbiorniki[i]
            emiter_x, emiter_y = emiters[liczba_emiterow - i - 1]
            sin_img[krok][i] = bresenham(emiter_x, emiter_y, odbiornik_x, odbiornik_y, img, img2, img_size)
        #if krok%10 == 0:
        #    io.imsave("./bresenham/kat"+str(krok*kat_kroku)+".jpg", img2)
    return sin_img

def inv_radon(img_sin):
    inv_sin = img_sin
    return inv_sin

def sinogram_circle_to_square(sinogram):
    diagonal = int(np.ceil(np.sqrt(2) * sinogram.shape[0]))
    pad = diagonal - sinogram.shape[0]
    old_center = sinogram.shape[0] // 2
    new_center = diagonal // 2
    pad_before = new_center - old_center
    pad_width = ((pad_before, pad - pad_before), (0, 0))
    return np.pad(sinogram, pad_width, mode='constant', constant_values=0)

def iradon_transform(radon_image, theta=None,interpolation='linear'):
    output_size = radon_image.shape[0]
    radon_image = sinogram_circle_to_square(radon_image)
    th = (np.pi / 180.0) * np.asarray(theta)
    # resize image to next power of two (but no less than 64) for
    # Fourier analysis; speeds up Fourier and lessens artifacts
    projection_size_padded = \
        max(64, int(2 ** np.ceil(np.log2(2 * radon_image.shape[0]))))
    pad_width = ((0, projection_size_padded - radon_image.shape[0]), (0, 0))
    img = np.pad(radon_image, pad_width, mode='constant', constant_values=0)
    f = fftfreq(projection_size_padded).reshape(-1, 1)   # digital frequency
    omega = 2 * np.pi * f                                # angular frequency
    fourier_filter = 2 * np.abs(f)                       # ramp filter
    projection = fft(img, axis=0) * fourier_filter
    radon_filtered = np.real(ifft(projection, axis=0))
    radon_filtered = radon_filtered[:radon_image.shape[0], :]
    reconstructed = np.zeros((output_size, output_size))
    mid_index = radon_image.shape[0] // 2
    [X, Y] = np.mgrid[0:output_size, 0:output_size]
    xpr = X - int(output_size) // 2
    ypr = Y - int(output_size) // 2
    # Reconstruct image by interpolation
    interpolation_types = ('linear', 'nearest', 'cubic')
    if interpolation not in interpolation_types:
        raise ValueError("Unknown interpolation: %s" % interpolation)
    for i in range(len(theta)):
        t = ypr * np.cos(th[i]) - xpr * np.sin(th[i])
        x = np.arange(radon_filtered.shape[0]) - mid_index
        if interpolation == 'linear':
            backprojected = np.interp(t, x, radon_filtered[:, i],
                                      left=0, right=0)
        else:
            interpolant = interp1d(x, radon_filtered[:, i], kind=interpolation,
                                   bounds_error=False, fill_value=0)
            backprojected = interpolant(t)
        reconstructed += backprojected
    radius = output_size // 2
    reconstruction_circle = (xpr ** 2 + ypr ** 2) <= radius ** 2
    reconstructed[~reconstruction_circle] = 0.

    return reconstructed * np.pi / (2 * len(th))

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
        sinogram = resize(sinogram, (256, 256))
        #io.imsave("sinogram.jpg",sinogram)
        plt.title('Sinogram')
        io.imshow(sinogram, cmap="gray")

        plt.subplot(2, 2, 4)
        re_img = iradon_transform(sinogram, theta=[kat_kroku * q for q in range(math.floor(180/kat_kroku))])
        plt.title('Obraz wyjsciowy')
        io.imshow(re_img, cmap='gray')

        io.show()
    else:
        print("Obrazek nie jest kwadratowy")


if __name__ == "__main__":
    main()

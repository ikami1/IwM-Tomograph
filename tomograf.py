from __future__ import division
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq
from scipy.interpolate import interp1d
import math
import time

imagename='phantom256.gif'

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


def plotLineHigh(x0, y0, x1,y1, img, img2, img_size):
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


def bresenham(x0, y0, x1, y1, img, img2, img_size):
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


def inversed_bresenham(x1,y1,x2,y2,sin,img2, img2_licznik,img_size):
    x1 = round(x1)
    y1 = round(y1)
    x2 = round(x2)
    y2 = round(y2)
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
                 img2[x, y] = img2[x, y] + sin
                 img2_licznik[x, y] += 1
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
                 img2[x,y]=img2[x,y]+sin
                 img2_licznik[x,y] += 1


#Polozenie detektorow
def odbiornik(r, nr_odbiornika, alfa=0.0):
    x = r * math.cos(alfa - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    y = r * math.sin(alfa - rozpietosc_katowa_punktow/2 + nr_odbiornika*rozpietosc_katowa_punktow/(liczba_emiterow-1))
    return x+r, y+r     #srodek nie w (0,0) a (r,r)


# Pojedynczy emiter
def emiter(r, alfa=0.0):
    x = r * math.cos(alfa)
    y = r * math.sin(alfa)
    return x+r, y+r     #srodek nie w (0,0) a (r,r)


# Tyle samo emiterow co detektorow
def emitery(r, nr_odbiornika, alfa=0.0):
    x = r * math.cos(alfa + math.pi - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))
    y = r * math.sin(alfa + math.pi - rozpietosc_katowa_punktow / 2 + nr_odbiornika * rozpietosc_katowa_punktow / (liczba_emiterow - 1))

    return x+r, y+r     #srodek nie w (0,0) a (r,r)


def get_sinogram(img, img_size):
    theta = np.linspace(0., 180., math.floor(180/kat_kroku), endpoint=False)
    r = img_size/2

    sin_img = np.zeros((liczba_emiterow, len(theta)))
    index = 0
    for step in theta:
        emiters = []
        for i in range(liczba_emiterow):
            emiters.append(emitery(r=r, nr_odbiornika=i, alfa=math.radians(step)))

        odbiorniki = []
        for i in range(liczba_emiterow):
            odbiorniki.append(odbiornik(r=r, nr_odbiornika=i, alfa=math.radians(step)))

        img2 = np.copy(img)
        for i in range(liczba_emiterow):
            odbiornik_x, odbiornik_y = odbiorniki[i]
            emiter_x, emiter_y = emiters[liczba_emiterow - i - 1]
            sin_img[i][index] = bresenham(emiter_x, emiter_y, odbiornik_x, odbiornik_y, img, img2, img_size)
        index += 1
        #if krok%10 == 0:
        #    io.imsave("./bresenham/kat"+str(krok*kat_kroku)+".jpg", img2)
    return sin_img


def inversed_Radon(radon_image, theta, ax3, ax4):
    output_size = radon_image.shape[0]
    th = (np.pi / 180.0) * theta    # radians

    # uzupelnij zerami zeby rozmiar byl: potega dwojki <= 64
    projection_size_padded = max(64, int(2 ** np.ceil(np.log2(2 * radon_image.shape[0]))))
    pad_width = ((0, projection_size_padded - radon_image.shape[0]), (0, 0))
    img = np.pad(radon_image, pad_width, mode='constant', constant_values=0)

    f = fftfreq(projection_size_padded).reshape(-1, 1)   # digital frequency, discrete Fourier transform
    fourier_filter = 2 * np.abs(f)                       # ramp filter

    projection = fft(img, axis=0) * fourier_filter
    radon_filtered = np.real(ifft(projection, axis=0))
    radon_filtered = radon_filtered[:radon_image.shape[0], :]

    ax3.set_title("Filtered Radon")
    ax3.imshow(radon_filtered, cmap="gray")
    plt.show(block=False)
    plt.pause(0.01)

    reconstructed = np.zeros((output_size, output_size))
    mid_index = radon_image.shape[0] // 2
    [X, Y] = np.mgrid[0:output_size, 0:output_size]
    xpr = X - int(output_size) // 2
    ypr = Y - int(output_size) // 2

    # Reconstruct image by interpolation
    for i in range(len(theta)):
        t = ypr * np.cos(th[i]) - xpr * np.sin(th[i])
        x = np.arange(radon_filtered.shape[0]) - mid_index
        interpolant = interp1d(x, radon_filtered[:, i], kind='cubic',
                               bounds_error=False, fill_value=0)
        backprojected = interpolant(t)
        reconstructed += backprojected

        if (i%5 == 0):
            ax4.imshow(reconstructed, cmap="gray")
            plt.show(block=False)
            plt.pause(0.01)

    radius = output_size // 2
    reconstruction_circle = (xpr ** 2 + ypr ** 2) <= radius ** 2
    reconstructed[~reconstruction_circle] = 0.

    return reconstructed * np.pi / (2 * len(th))


image = misc.imread(imagename, flatten=True).astype('float64')
img_size = max(image.shape)
image = misc.imresize(image, (img_size, img_size))
sinogram = get_sinogram(image, img_size)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 6))

ax1.set_title("Original image")
ax1.imshow(image, cmap="gray")
ax2.set_title("Radon transform\n(Sinogram)")
ax2.imshow(sinogram, cmap="gray", aspect='auto')
ax4.set_title("Reconstruction")

reconstruction_fbp = inversed_Radon(sinogram, np.linspace(0., 180., math.floor(180/kat_kroku), endpoint=False), ax3, ax4)
reconstruction_fbp = misc.imresize(reconstruction_fbp, (img_size,img_size))

error = reconstruction_fbp - image
print('FBP rms reconstruction error: %.3g' % np.sqrt(np.mean(error**2)))

ax3.set_title("Filtered Reconstruction")
ax3.imshow(reconstruction_fbp, cmap="gray")
plt.show()
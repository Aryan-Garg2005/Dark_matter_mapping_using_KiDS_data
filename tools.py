import numpy as np

def apodize(map2d, frac=0.15):
    ny, nx = map2d.shape
    wx = np.ones(nx)
    wy = np.ones(ny)

    nedge_x = int(frac * nx)
    nedge_y = int(frac * ny)

    x = np.linspace(0, np.pi/2, nedge_x)
    taper_x = np.sin(x)**2
    wx[:nedge_x] = taper_x
    wx[-nedge_x:] = taper_x[::-1]

    y = np.linspace(0, np.pi/2, nedge_y)
    taper_y = np.sin(y)**2
    wy[:nedge_y] = taper_y
    wy[-nedge_y:] = taper_y[::-1]

    window = np.outer(wy, wx)
    return map2d * window

def kaiser_squires(g1, g2):
    g1_hat = np.fft.fft2(g1)
    g2_hat = np.fft.fft2(g2)

    ky = np.fft.fftfreq(g1.shape[0])
    kx = np.fft.fftfreq(g1.shape[1])
    kx, ky = np.meshgrid(kx, ky)

    k2 = kx**2 + ky**2
    k2[0,0] = 1

    D1 = (kx**2 - ky**2)/k2
    D2 = (2*kx*ky)/k2

    kappa_hat = D1*g1_hat + D2*g2_hat
    return np.real(np.fft.ifft2(kappa_hat))

def kaiser_squires_EB(g1, g2):
    g1_hat = np.fft.fft2(g1)
    g2_hat = np.fft.fft2(g2)

    ky = np.fft.fftfreq(g1.shape[0])
    kx = np.fft.fftfreq(g1.shape[1])
    kx, ky = np.meshgrid(kx, ky)

    k2 = kx**2 + ky**2
    k2[0, 0] = 1

    D1 = (kx**2 - ky**2) / k2
    D2 = (2 * kx * ky) / k2

    kappa_E_hat = D1 * g1_hat + D2 * g2_hat
    kappa_B_hat = -D2 * g1_hat + D1 * g2_hat

    kappa_E = np.real(np.fft.ifft2(kappa_E_hat))
    kappa_B = np.real(np.fft.ifft2(kappa_B_hat))

    return kappa_E, kappa_B


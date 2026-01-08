from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from tools import apodize, kaiser_squires, kaiser_squires_EB

hdul = fits.open(r"weak_lensing_data.fits",memmap = True)

data = hdul[1].data

ra = data['RAJ2000']
dec = data['DECJ2000']
g1 = data['e1']
g2 = data['e2']
w = data['weight']

patch = (
    (ra > 215.0) & (ra < 217.0) &
    (dec > -2.0) & (dec < 0)
)

ra, dec = ra[patch], dec[patch]
g1, g2 = g1[patch], g2[patch]
w = w[patch]

good = (w > 0) & np.isfinite(g1) & np.isfinite(g2)

ra, dec = ra[good], dec[good]
g1, g2 = g1[good], g2[good]
w = w[good]

print(f"Selected galaxies: {len(ra)}")

#Conversion to cartesian coordinates
ra0  = np.mean(ra)
dec0 = np.mean(dec)

x = (ra - ra0) * np.cos(np.deg2rad(dec0))
y = (dec - dec0)

N = 128

xbins = np.linspace(x.min(), x.max(), N+1)
ybins = np.linspace(y.min(), y.max(), N+1)

gamma1 = np.zeros((N, N))
gamma2 = np.zeros((N, N))
weight = np.zeros((N, N))

ix = np.digitize(x, xbins) - 1
iy = np.digitize(y, ybins) - 1

valid = (ix >= 0) & (ix < N) & (iy >= 0) & (iy < N)

for i in np.where(valid)[0]:
    gamma1[iy[i], ix[i]] += w[i] * g1[i]
    gamma2[iy[i], ix[i]] += w[i] * g2[i]
    weight[iy[i], ix[i]] += w[i]

mask = weight > 0
gamma1[mask] /= weight[mask]
gamma2[mask] /= weight[mask]
gamma1[~mask] = np.nan
gamma2[~mask] = np.nan

gamma1 = np.nan_to_num(gamma1)
gamma2 = np.nan_to_num(gamma2)

gamma1 -= np.mean(gamma1)
gamma2 -= np.mean(gamma2)

kappa = kaiser_squires(gamma1, gamma2)

kappa_smooth = gaussian_filter(kappa, sigma=6.0)

gamma1_ap = apodize(gamma1)
gamma2_ap = apodize(gamma2)

kE, kB = kaiser_squires_EB(gamma1_ap, gamma2_ap)

kE_vis = gaussian_filter(kE, sigma=6.0)
kB_vis = gaussian_filter(kB, sigma=6.0)
print("E-mode std:", np.std(kE_vis))
print("B-mode std:", np.std(kB_vis))
print("B/E ratio:", np.std(kB_vis) / np.std(kE_vis))

vmin, vmax = np.percentile(kappa_smooth, [5, 95])
plt.figure(figsize=(5,5))
plt.imshow(kappa_smooth, origin='lower', cmap='inferno', vmin = vmin, vmax = vmax)
plt.colorbar(label=r'$\kappa$')
plt.title("KiDS-450 G15 Weak-Lensing Convergence Map")
plt.axis('off')
plt.show()


